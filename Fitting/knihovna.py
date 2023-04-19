import os
from tkinter import filedialog
from tkinter import *
from tkinter.messagebox import *
from tkinter import messagebox
import tkinter.simpledialog as sm
import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

# Upozorneni !!!
# Momentalni verze funguje pouze pro pole v x-ovem smeru !!!

class spectra:
    def __init__(self,spectra,path,RT,repeat,time,angle,field,temperature,fc):
        self.spec = spectra # dataframe spekter
        self.path = path # cesta k souboru
        self.RT = RT     # Real Time, skutecny cas mereni
        self.repeat = repeat
        self.time = time
        self.angle = angle
        self.field = field
        self.temp = temperature
        self.fc = fc     # field current

    def __str__(self):
        return 'Real time      '+ str(self.RT) + '\n Repeat        ' + str(self.repeat)+ '\n Time          ' + str(self.time)+ '\n Angle         ' + str(self.angle)+ '\n Field         ' + str(self.field) + '\n Temperature   ' + str(self.temp)+ '\n Field current '+str(self.fc)

class measurement:
    def __init__(self,spec,repeat,time,field,temp,findmin,sample,folder_path):
        self.list = spec              # List spekter v tomto mereni
        self.repeat = repeat          # Repeat
        self.time = time              # Time (repetice ne real time)
        self.field = field            # Pole
        self.temp = temp              # Teplota
        self.findmin = findmin        # jsou to findmin spektra?
        self.sample = sample          # sample name
        self.folder = folder_path     # path do slozky s merenimi
        self.fitted = 'Not defined'   # dataframe nafitovanych dat
        self.fittable = 'Not defined' # bulik jestli jsou to RPE spektra

    def __str__(self):
        return '# Vzorek       ' + str(self.sample)+'\n# Repeat       ' + str(self.repeat)+ '\n# Time         ' + str(self.time)+ '\n# Findmin      ' + str(self.findmin)+ '\n# Field        ' + str(self.field) + '\n# Temperature  ' + str(self.temp)+ '\n# Fittable     '+str(self.fittable)+ '\n# Folder path  '+str(self.folder)
    
    def standart_moke_str(self):
        return '# Vzorek:       ' + str(self.sample)+'\n# Repeat:       ' + str(self.repeat)+ '\n# Time:         ' + str(self.time)+ '\n# Field:        ' + str(abs(self.field)) + '\n# Temperature:  ' + str(self.temp)+ '\n# Folder path:  '+str(self.folder)

    def check_duplicates(self):
        angles = []
        for i in range(0,len(self.list)):
            if self.list[i].angle in angles:
                return True
            else:
                angles.append(self.list[i].angle)
        return False

    def deduplicate(self):
        angles = []
        measurements = []
        for i in range(0,len(self.list)):
            if self.list[i].angle in angles:
                measurements.append(self.list[0:i])
                angles = [self.list[i].angle]
            else:
                angles.append(self.list[i].angle)
        measurements.append(self.list[-len(angles):])
        return measurements

    def get_spectra(self):
        # Pokud se meni uhel, pak nazev sloupcu jsou uhly
        # Pokud se nehybe polarizatorem, pak nazvy sloupcu jsou RT
        if self.list[0].angle != 'Not defined':
            data = self.list[0].spec
            data = data.rename(columns = {1:self.list[0].angle})
            for spec in self.list[1:]:
                data = pd.concat([data,spec.spec.rename(columns = {1:spec.angle})],axis=1, join='inner')
        else:
            data = self.list[0].spec
            data = data.rename(columns = {1:self.list[0].RT})
            for spec in self.list[1:]:
                data = pd.concat([data,spec.spec.rename(columns = {1:spec.RT})],axis=1, join='inner')
        
        data = data.sort_index()
        data = data.interpolate(method = 'index')
        data = data[~data.index.duplicated(keep='first')] # drop duplicate indexes
        data = data[~data.isin([np.nan]).any(axis=1)] 

        return data
    
    def is_fitable(self):
        is_fit = True
        self.fittable = True
        for meas in self.list:
            if meas.angle == 'Not defined':
                is_fit = False
                self.fittable = False
        if self.field == 'Not defined':
            is_fit = False
            self.fittable = False
        if self.findmin:
            is_fit = False
            self.fittable = False
        
        return is_fit
    
    def fit(self): # TODO jak se pocita chyba ??
        if self.is_fitable():
            data = self.get_spectra()
            angles = list(data.columns)
            A = []
            for angle in np.deg2rad(angles):
                A.append([np.power(np.sin(angle),2),np.sin(2*angle),1])
            A = np.linalg.pinv(A)
            for ind,row in data.iterrows():
                data.loc[ind,'Fit'] = np.rad2deg(A.dot(row.values)[1]/A.dot(row.values)[0])
            data['eV'] = 1239.8/data.index.to_series()
            data = data.set_index('eV')
        else:
            data = 'Not fittable'

        self.fitted = data
        return data

class measurement_result:
    def __init__(self,data,repeat,time,field,temp,sample,folder_path,date,type_meas,path):
        self.repeat = repeat          # Repeat
        self.time = time              # Time (repetice ne real time)
        self.field = field            # Pole
        self.temp = temp              # Teplota
        self.sample = sample          # sample name
        self.folder = folder_path     # path do slozky s merenimi
        self.data = data              # dataframe dat
        self.date = date              # Date measured
        self.type = type_meas         # Measurement type eg loop, moke,one polarity of field
        self.path = path              # Path souboru
        if 'No BG' in type_meas:
            self.bg_sub = True
        else:
            self.bg_sub = False

    def __str__(self):
        if self.bg_sub:
            return '# Vzorek:       ' + str(self.sample)+'\n# Repeat:       ' + str(self.repeat)+ '\n# Time:         ' + str(self.time)+ '\n# Field:        ' + str(self.field) + '\n# Temperature:  ' + str(self.temp)+  '\n# Folder path:  '+str(self.folder)+'\n# Background:   Subtracted'
        else:
            return '# Vzorek:       ' + str(self.sample)+'\n# Repeat:       ' + str(self.repeat)+ '\n# Time:         ' + str(self.time)+ '\n# Field:        ' + str(self.field) + '\n# Temperature:  ' + str(self.temp)+  '\n# Folder path:  '+str(self.folder)
class Full_moke:
    def __init__(self,data,field,temp,sample,path,date,type_meas):
        self.field = field            # Pole
        self.temp = temp              # Teplota
        self.sample = sample          # sample name
        self.path = path              # pathy na soubory jednotlivych zpracovanych mereni
        self.data = data              # dataframe dat
        self.date = date              # Dates measured, list
        self.type = type_meas         # Measurement type eg loop, moke,one polarity of field

    def print_fm(self): # print full moke
        dates = ''
        for i in self.date:
            if i not in dates:
                dates = dates +'\n                      '
        rot = ''
        for i in self.path:
            rot = rot+'\n#               '+i
        return '# Type of measurement: '+self.type + '# Vzorek:       ' + str(self.sample)+'\n# Field:        ' + str(abs(self.field)) + '\n# Temperature:  ' + str(self.temp)+'\n# Dates measured: '+dates+ '\n# Raw data:      '+rot+'\n'

    def print_av(self): # Print average
        dates = ''
        for i in self.date:
            if i not in dates:
                dates = dates +'\n                      '
        pathy = ''
        for i in self.path:
            pathy = pathy+'\n#               '+i
        return '# Type of measurement: '+self.type + '# Vzorek:       ' + str(self.sample)+'\n# Field:        ' + str(abs(self.field)) + '\n# Temperature:  ' + str(self.temp)+'\n# Dates measured: '+dates+ '\n# Paths:    '+pathy+'\n'


def help(funkce = None):  # TODO prumer, bg
    fitni_help = '''- fitni(show = True,average = False, show_fit = 4.0)
          - Fitne data ve vybrané složce a to i v podsložkách
            - Parameter \'show\' udává, zda chcete data rovnou ukázat
            - Parameter \'average\' říká, zda chcete data i zprůměrovat
              - MOMENTÁLNĚ NEFUNKČNÍ
            - Parameter \'show_fit\' udává, zda chcete ukázat konkrétní fity
              - Přijímá jako hodnotu energii na které vám ukáže fit 
              - MOMENTÁLNĚ NEFUNKČNÍ
          '''
    uka_help = '''- uka(path = \'Not defined\',num = 1,sep = 1,pm = False, loop_energy = \'Not defined\')
          - Ukáže MOKE/loop vybranych souboru
            - Parameter \'sep\' rozhoduje, jak se data budou shlukovat do grafů:
              - 0: Všechna data v jednom grafu
              - 1: Všechny MOKE data v jednom grafu a všechny loop data v jednom grafu a stejně one polarity of field
              - 2: V jednom grafu vždy všechna měření stejného vzorku a stejného typu měření
              - 3: Každé měření ve svém grafu 
            - Parameter \'num\' udává počet složek, ze kterých chcete měření brát
            - Parameter \'pm\' udává, zda chceme vidět i plus/mínus pole u MOKE fitů
            - Parameter \'path\' je defaultně prázdný, pokud ho určíte pak se Tkinter okno otevře v zadané složce
            - Parameter \'loop_energy\' je pole energií na kterých udělat řez smyčky
              Pokud je elementem pole další pole o dvou prvcích, pak dostanete integrální řez mezi těmito energiemi
    ''' 
    uka_vse_help = '''- uka_vse(path = \'Not defined\',num = 1,pm = False, loop_energy = \'Not defined\')
          - Ukáže všechna MOKE/loop měření ve vybraném adresáři a podadresářích
            - Parameter \'num\' udává počet složek, ze kterých chcete měření brát
            - Parameter \'path\' je volitelný, když zadáte, pak se fce vyhodnotí v zadané složce
              pokud jej nezadáte, otevře se vám Tkinter okno
            - Parameter \'pm\' udává, zda chceme vidět i plus/mínus pole u MOKE fitů
            - Parameter \'loop_energy\' je pole energií na kterých udělat řez smyčky
              Pokud je elementem pole další pole o dvou prvcích, pak dostanete integrální řez mezi těmito energiemi
    '''
    elip_help = '''- elip(path = 'Not defined',show = True, save = True)
          - Spočte elipticitu z měření rotace a desek
            - Parameter \'path\' udává výchozí složku ve které se otevře Tkinter okno
            - Parameter \'show\' udává, zda chcete zobrazit výsledný graf
            - Parameter \'save\' udává, zda chcete výsledná data uložit
    '''
    
    fce = {
        "fitni": fitni_help,
        "uka": uka_help,
        "uka_vse": uka_vse_help,
        "elip": elip_help
        }
    if funkce == None:
        here_you_go = '\nV tomto programu existují následující funkce:\n'
        options = '''fitni\nuka\nuka_vse\nelip\n
- Pokud chcete info o libovolné z funkcí zadejte:
        help('jméno_funkce')
        '''
        print(here_you_go)
        print(options)
    else:
        try:
            print(fce[funkce])
        except:
            print('No such function in this module !!!')

def LoadData(): 
    # Load all txt files in chosen directory
    root = Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    initial_dir = os.getcwd()
    folder = filedialog.askdirectory(initialdir = initial_dir,title = "Select measurement folders")
    root.destroy()
    Path_folder = pathlib.Path(folder)
    files = list(Path_folder.rglob("*.txt"))
    files.sort(key = lambda x: str(x).rsplit('\\',1)[1][0:18])
    # Initialize a list of all spectra and of all measurement configurations
    # spectra_list is a 2D list, for each measurement it contains a list of all spectra
    spectra_list = []
    measurements = []
    last_meas = []

    # for each path to a txt measurement file disect its name into components
    # and load it into the spectra class
    print('\nLoading files')
    for file in tqdm(files):
        name = str(file).rsplit('\\',1)[1]
        path = str(file).rsplit('\\',1)[0]
        try:
            data = pd.read_csv(file, comment='#', index_col=0, sep='\t', header=None)
        except:
            data = 'Not defined'
        RT = name[0:18]
        name = name.split('_')
        name[-1] = name[-1].split('.txt')[0]
        things = ['repeat','time','polarizer','fieldx','temperature','fieldcurrentx']
        properties = []
        for i in things:
            if i in name:
                index = name.index(i)
                properties.append(float(name[index+1]))
            else:
                properties.append('Not defined')
        if 'findmin' in name:
            findmin = True
        else:
            findmin = False
        meas = [properties[0],properties[1],properties[3],properties[4],findmin,name[2],path]
        if meas == last_meas:
            spectra_list[-1].append(spectra(data,file,RT,properties[0],properties[1],properties[2],properties[3],properties[4],properties[5]))
        else:
            measurements.append(meas)
            last_meas = meas
            spectra_list.append([spectra(data,file,RT,properties[0],properties[1],properties[2],properties[3],properties[4],properties[5])])

    # Check if there are duplicate measurements in one element (e.g. -1 T during loop measurement)
    measurement_list = []
    for i in range(0,len(measurements)):
        current = measurement(spectra_list[i],measurements[i][0],measurements[i][1],measurements[i][2],measurements[i][3],measurements[i][4],measurements[i][5],measurements[i][6])
        if current.check_duplicates():
            duplicates = current.deduplicate()
            for j in duplicates:
                deduplicate = measurement(j,measurements[i][0],measurements[i][1],measurements[i][2],measurements[i][3],measurements[i][4],measurements[i][5],measurements[i][6])
                measurement_list.append(deduplicate)
        else:
            measurement_list.append(current)
            
    # Return a list of measurement class objects and opened folder path
    return measurement_list, folder

def fitni(show = True): # TODO path, show_fit
    data,path = LoadData()
    print('\nFitting')
    for i in tqdm(range(0,len(data))):
        data[i].fit()
    splacnuto = []
    done = []
    for i in range(0,len(data)):
        if data[i].fittable and i not in done:
            splacnuto.append([data[i]])
            for j in range(i+1,len(data)):
                if data[i].repeat == data[j].repeat and data[i].time == data[j].time and data[i].temp == data[j].temp and data[i].sample == data[j].sample and data[i].folder == data[j].folder:
                    splacnuto[-1].append(data[j])
                    done.append(j)
    for mereni in splacnuto:
        if len(mereni) == 1:
            save_one_polarity(mereni,path)
        elif len(mereni) == 2:
            save_kerr(mereni,path)
        else:
            save_smycka(mereni,path)
    if show:
        uka_vse(path)

def save_one_polarity(data,path):
    path = path + '\\' + data[0].sample
    if data[0].temp != 'Not defined':
        path = path + '_Temperature_' + str(data[0].temp) + 'K'
    if data[0].repeat != 'Not defined':
        path = path + '_Repeat_' + str(data[0].repeat)
    if data[0].time != 'Not defined':
        path = path + '_Time_' + str(data[0].time)
    if data[0].field != 'Not defined':
        path = path + '_Field_' + str(data[0].field) + 'T'
    path = path +'.dat'

    if os.path.exists(path):
        if askyesno(title='File exists', message='File already exists. Do you want to overwrite it?\nOne polarity measurement\n'+str(data[0])):
            os.remove(path)
        else:
            return
    file = open(path,'a')
    if data[0].list[0].field != 'Not defined':
        file.write('# Type of measurement: One polarity of field at '+ str(data[0].field) + 'T\n' )
    else:
        file.write('# Type of measurement: One polarity of field at ??T\n')
    file.write('# Date measured:       ' + data[0].list[0].RT[6:8] + '.' + data[0].list[0].RT[4:6]+'.'+data[0].list[0].RT[:4]+'\n')
    file.write(str(data[0].standart_moke_str()))
    file.write('\n')
    data[0].fitted = data[0].fitted.rename(columns = {'Fit':'MOKE'})
    data[0].fitted.to_csv(file)
    file.close()

def save_kerr(data,path):
    path = path + '\\' + data[0].sample
    if data[0].temp != 'Not defined':
        path = path + '_Temperature_' + str(data[0].temp) + 'K'
    if data[0].repeat != 'Not defined':
        path = path + '_Repeat_' + str(data[0].repeat)
    if data[0].time != 'Not defined':
        path = path + '_Time_' + str(data[0].time)
    path = path +'.dat'

    if os.path.exists(path):
        if askyesno(title='File exists', message='File already exists. Do you want to overwrite it?\nStandart MOKE measurement\n'+data[0].standart_moke_str()):
            os.remove(path)
        else:
            return
    file = open(path,'a')
    file.write('# Type of measurement: Standart MOKE measurement at +- ' + str(abs(data[0].field)) + 'T\n' )
    file.write('# Date measured:       ' + data[0].list[0].RT[6:8] + '.' + data[0].list[0].RT[4:6]+'.'+data[0].list[0].RT[:4]+'\n')
    file.write(data[0].standart_moke_str())
    file.write('\n')
    file.close()
    data[0].fitted = data[0].fitted.rename(columns = {'Fit':'MOKE '+str(data[0].field)+'T'})
    data[1].fitted = data[1].fitted.rename(columns = {'Fit':'MOKE '+str(data[1].field)+'T'})
    data_save = pd.concat([data[0].fitted[['MOKE '+str(data[0].field)+'T']],data[1].fitted[['MOKE '+str(data[1].field)+'T']]],axis=1, join='inner')
    data_save = data_save.sort_index()
    data_save = data_save.interpolate(method = 'index')
    data_save = data_save[~data_save.index.duplicated(keep='first')] # drop duplicate indexes
    data_save = data_save[~data_save.isin([np.nan]).any(axis=1)] 
    data_save['MOKE'] = (data[0].fitted['MOKE '+str(data[0].field)+'T'] - data[1].fitted['MOKE '+str(data[1].field)+'T'])/2
    data_save.to_csv(path,mode='a')

def save_smycka(data,path):
    path = path + '\\' + data[0].sample + '_smycka'
    if data[0].temp != 'Not defined':
        path = path + '_Temperature_' + str(data[0].temp) + 'K'
    if data[0].repeat != 'Not defined':
        path = path + '_Repeat_' + str(data[0].repeat)
    if data[0].time != 'Not defined':
        path = path + '_Time_' + str(data[0].time)
    path = path +'.dat'

    if os.path.exists(path):
        if askyesno(title='File exists', message='File already exists. Do you want to overwrite it?\Loop MOKE measurement\n'+data[0].standart_moke_str()):
            os.remove(path)
        else:
            return
    file = open(path,'a')
    file.write('# Type of measurement: Loop MOKE measurement\n' )
    file.write('# Date measured:       ' + data[0].list[0].RT[6:8] + '.' + data[0].list[0].RT[4:6]+'.'+data[0].list[0].RT[:4]+'\n')
    max_field = sorted(data,key = lambda x: x.field)
    file.write(max_field[0].standart_moke_str())
    file.write('\n')
    file.close()

    # Prejmenovavani sloupcu. Sloupce maji stejna jmena pri stejnem poli. 
    # Pri pd_readcsv se nactou jako 1T a 1T.1 atd.
    data.sort(key = lambda x: int(x.list[0].RT.replace('_','')))
    data_save = data[0].fitted[['Fit']]
    data_save = data_save.rename(columns = {'Fit': str(data[0].field) + 'T'})
    # name_colum = [str(data[0].field) + 'T']

    # Zakomentovana cast kdyztak prejnemovava sloupce aby se nejmenovali stejne
    for i in range(1,len(data)):
        name_col_single = str(data[i].field) + 'T'
        # j = 2     
        # while name_col_single in name_colum:
        #     name_col_single = str(data[i].field) + 'T ' + str(j)
        #     j = j + 1   
        # name_colum.append(name_col_single)
        data[i].fitted = data[i].fitted.rename(columns = {'Fit': name_col_single})
        data_save = pd.concat([data_save,data[i].fitted[[name_col_single]]],axis=1, join='inner')

    data_save = data_save.sort_index()
    data_save = data_save.interpolate(method = 'index')
    data_save = data_save[~data_save.index.duplicated(keep='first')] # drop duplicate indexes
    data_save = data_save[~data_save.isin([np.nan]).any(axis=1)] 
    data_save.to_csv(path,mode='a')
    
def ukaz_LoadData(filenames): # TODO savitzky golay, background odecet
    # Nacte soubory jiz vyfitovanych dat
    one_pol = []
    moke = []
    loop = []
    fm = []
    pm = []
    for file in filenames:
        try:
            path = file
            data = pd.read_csv(file, comment='#', index_col=0)
            f = open(file,'r')
            type_meas = f.readline().rsplit(': ')[1].strip()
            if 'Full MOKE' in type_meas: 
                vzorek = f.readline().rsplit(': ')[1].strip()
                field = f.readline().rsplit(': ')[1].strip()
                temp = f.readline().rsplit(': ')[1].strip()
                dates = []
                dates.append(f.readline().rsplit(': ')[1].strip())
                lajna = f.readline()
                while 'data' not in lajna:
                    dates.append(lajna.rsplit('# ')[1].strip())
                    lajna = f.readline()
                path = []
                path.append(lajna.rsplit(': ')[1].strip())
                lajna = f.readline()
                while '#' in lajna:
                    path.append(lajna.rsplit('# ')[1].strip())
                    lajna = f.readline()
            elif 'Prumer' in type_meas:
                vzorek = f.readline().rsplit(': ')[1].strip()
                field = f.readline().rsplit(': ')[1].strip()
                temp = f.readline().rsplit(': ')[1].strip()
                dates = []
                dates.append(f.readline().rsplit(': ')[1].strip())
                lajna = f.readline()
                while 'Path' not in lajna:
                    dates.append(lajna.rsplit('# ')[1].strip())
                    lajna = f.readline()
                path = []
                path.append(lajna.rsplit(': ')[1].strip())
                lajna = f.readline()
                while '#' in lajna:
                    path.append(lajna.rsplit('# ')[1].strip())
                    lajna = f.readline()
            else:
                date = f.readline().rsplit(': ')[1].strip()
                vzorek = f.readline().rsplit(': ')[1].strip()
                repeat = f.readline().rsplit(': ')[1].strip()
                time = f.readline().rsplit(': ')[1].strip()
                field = f.readline().rsplit(': ')[1].strip()
                temp = f.readline().rsplit(': ')[1].strip()
                folder = f.readline().rsplit(':  ')[1].strip()
                f.close()
        except:
            data = 'Not defined'
            f = 'Not defined'
            type_meas = 'Not defined'
            date = 'Not defined'
            vzorek = 'Not defined'
            repeat = 'Not defined'
            time = 'Not defined'
            field = 'Not defined'
            temp = 'Not defined'
            folder = 'Not defined'
        if 'One polarity' in type_meas:
            one_pol.append(measurement_result(data,repeat,time,field,temp,vzorek,folder,date,type_meas,path))
        elif 'Standart MOKE' in type_meas:
            moke.append(measurement_result(data,repeat,time,field,temp,vzorek,folder,date,type_meas,path))
        elif 'Loop MOKE' in type_meas:
            loop.append(measurement_result(data,repeat,time,field,temp,vzorek,folder,date,type_meas,path))
        elif 'Full MOKE' in type_meas:
            fm.append(Full_moke(data,field,temp,vzorek,path,dates,type_meas))
        elif 'Prumer' in type_meas:
            pm.append(Full_moke(data,field,temp,vzorek,path,dates,type_meas))
    return one_pol, moke, loop,fm,pm

def uka(path = 'Not defined',num = 1,sep = 1,pm = False, loop_energy = 'Not defined'):
    # Tato funkce ukaze konkretni mereni v adresari, ktere uzivatel vybere
    if path == 'Not defined':
        initial_dir = os.getcwd()
    else:
        initial_dir = path
    filenames = []
    for i in range(num):
        root = Tk()
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        filenames = filenames + list(filedialog.askopenfilename(initialdir = initial_dir,title = "Select files",filetypes = (("txt files","*.dat .txt .KNT"),("all files","*.*")),multiple=True))
    root.destroy()
    onepol,moke,loop,fm,pm = ukaz_LoadData(filenames)
    color_list = [plt.cm.coolwarm,plt.cm.BrBG,plt.cm.PiYG,plt.cm.PuOr,plt.cm.PRGn,plt.cm.twilight_shifted]
    color_list = color_list + max(0,len(color_list)-len(loop))*[plt.cm.coolwarm]
    if sep == 0:
        fig, ax = plt.subplots()
        for meas in onepol:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            if meas.field != 'Not defined':
                legenda = legenda + ' Field ' + str(meas.field)
            ax.plot(meas.data['MOKE'], label = legenda)
        for meas in moke:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            if pm:
                names = meas.data.columns
                for name in names:
                    ax.plot(meas.data[name], label = legenda+' '+str(name))
            else:
                ax.plot(meas.data['MOKE'], label = legenda)
        for meas in loop:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            colors = color_list[loop.index(meas)](np.linspace(0, 1, len(meas.data.columns)))
            for i, col in enumerate(meas.data.columns):
                ax.plot(meas.data.index, meas.data[col], color=colors[i])
            cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=color_list[loop.index(meas)]), ax=ax, ticks=[0, 0.5, 1])
            cbar.ax.invert_yaxis()
            cbar.ax.set_yticklabels([meas.data.columns[0].rsplit('T')[0]+' T', meas.data.columns[len(meas.data.columns)//2].rsplit('T')[0]+' T', meas.data.columns[-1].rsplit('T')[0]+' T'])
            ax.lines[-1].set_label(legenda)
        for meas in fm:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            names = meas.data.columns
            for name in names:
                ax.plot(meas.data[name], label = legenda+' '+str(name))
        for meas in pm:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            names = meas.data.columns
            for name in names:
                ax.plot(meas.data[name], label = legenda+' '+str(name))
        plt.title('All measurements')
    if sep == 1:
        for meas in onepol:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            if meas.field != 'Not defined':
                legenda = legenda + ' Field ' + str(meas.field)
            plt.figure(1)
            plt.plot(meas.data['MOKE'], label = legenda)
            plt.title('All one polarity of field measurements')
        for meas in moke:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            plt.figure(2)
            if pm:
                names = meas.data.columns
                for name in names:
                    plt.plot(meas.data[name], label = legenda+' '+str(name))
            else:
                plt.plot(meas.data['MOKE'], label = legenda)
            plt.title('All MOKE measurements')
        for meas in fm:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            plt.figure(3)
            names = meas.data.columns
            for name in names:
                plt.plot(meas.data[name], label = legenda+' '+str(name))
            plt.title('All Full MOKE measurements')
        for meas in pm:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            plt.figure(4)
            names = meas.data.columns
            for name in names:
                plt.plot(meas.data[name], label = legenda+' '+str(name))
            plt.title('All avereged MOKE measurements')
        if len(loop)>0:    
            fig, ax = plt.subplots() 
        for meas in loop:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                print(meas.temp)
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            colors = color_list[loop.index(meas)](np.linspace(0, 1, len(meas.data.columns)))
            for i, col in enumerate(meas.data.columns):
                ax.plot(meas.data.index, meas.data[col], color=colors[i])
            cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=color_list[loop.index(meas)]), ax=ax, ticks=[0, 0.5, 1])
            cbar.ax.invert_yaxis()
            cbar.ax.set_yticklabels([meas.data.columns[0].rsplit('T')[0]+' T', meas.data.columns[len(meas.data.columns)//2].rsplit('T')[0]+' T', meas.data.columns[-1].rsplit('T')[0]+' T'])
            ax.lines[-1].set_label(legenda)
            plt.title('All loop measurements')  
    if sep == 2:
        samples = []
        for meas in onepol:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            if meas.field != 'Not defined':
                legenda = legenda + ' Field ' + str(meas.field)
            if meas.sample not in samples:
                samples.append(meas.sample)
            plt.figure(samples.index(meas.sample))
            plt.plot(meas.data['MOKE'], label = legenda)
            plt.title('One polarity of field of '+ str(meas.sample))
        add = len(samples)
        samples = []
        for meas in moke:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            if meas.sample not in samples:
                samples.append(meas.sample)
            plt.figure(add+samples.index(meas.sample))
            plt.title('Standart MOKE measurement of '+ str(meas.sample))
            if pm:
                names = meas.data.columns
                for name in names:
                    plt.plot(meas.data[name], label = legenda+' '+str(name))
            else:
                plt.plot(meas.data['MOKE'], label = legenda)
        add = add + len(samples)
        samples = []
        for meas in fm:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.sample not in samples:
                samples.append(meas.sample)
            plt.figure(add+samples.index(meas.sample))
            plt.title('Full MOKE measurement of '+ str(meas.sample))
            names = meas.data.columns
            for name in names:
                plt.plot(meas.data[name], label = legenda+' '+str(name))
        add = add + len(samples)
        samples = []
        for meas in pm:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.sample not in samples:
                samples.append(meas.sample)
            plt.figure(add+samples.index(meas.sample))
            plt.title('Averaged MOKE measurement of '+ str(meas.sample))
            names = meas.data.columns
            for name in names:
                plt.plot(meas.data[name], label = legenda+' '+str(name))
        add = add + len(samples)
        samples = []
        figures = []

        loop.sort(key = lambda x: x.sample)
        iterab = 0
        for meas in loop:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            if meas.sample not in samples:
                samples.append(meas.sample)
                fig, ax = plt.subplots()
                figures.append(ax)
                iterab = 0
            colors = color_list[iterab](np.linspace(0, 1, len(meas.data.columns)))
            for i, col in enumerate(meas.data.columns):
                figures[len(samples)-1].plot(meas.data.index, meas.data[col], color=colors[i])
            cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=color_list[iterab]), ax=ax, ticks=[0, 0.5, 1])
            cbar.ax.invert_yaxis()
            cbar.ax.set_yticklabels([meas.data.columns[0].rsplit('T')[0]+' T', meas.data.columns[len(meas.data.columns)//2].rsplit('T')[0]+' T', meas.data.columns[-1].rsplit('T')[0]+' T'])
            ax.lines[-1].set_label(legenda)
            ax.set_title('Loop measurements of ' + meas.sample)
            iterab += 1
    if sep == 3:
        for meas in onepol:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            if meas.field != 'Not defined':
                legenda = legenda + ' Field ' + str(meas.field)
            plt.figure()
            plt.plot(meas.data['MOKE'], label = legenda)
            plt.title('One polarity of field measurement of '+ legenda)
        for meas in moke:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            plt.figure()
            if pm:
                names = meas.data.columns
                for name in names:
                    plt.plot(meas.data[name], label = legenda+' '+str(name))
            else:
                plt.plot(meas.data['MOKE'], label = legenda)
            plt.title('MOKE measurements of '+legenda) 
        for meas in fm:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            plt.figure()
            names = meas.data.columns
            for name in names:
                plt.plot(meas.data[name], label = legenda+' '+str(name))
            plt.title('Full MOKE measurements of '+legenda) 
        for meas in pm:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            plt.figure()
            names = meas.data.columns
            for name in names:
                plt.plot(meas.data[name], label = legenda+' '+str(name))
            plt.title('Averaged MOKE measurements of '+legenda) 
        for meas in loop:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                print(meas.temp)
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            fig, ax = plt.subplots()
            colors = color_list[0](np.linspace(0, 1, len(meas.data.columns)))
            for i, col in enumerate(meas.data.columns):
                ax.plot(meas.data.index, meas.data[col], color=colors[i])
            cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=color_list[0]), ax=ax, ticks=[0, 0.5, 1])
            cbar.ax.invert_yaxis()
            cbar.ax.set_yticklabels([meas.data.columns[0].rsplit('T')[0]+' T', meas.data.columns[len(meas.data.columns)//2].rsplit('T')[0]+' T', meas.data.columns[-1].rsplit('T')[0]+' T'])
            ax.lines[-1].set_label(legenda)
            plt.title('Loop measurement of '+legenda)  
    
    for i in plt.get_fignums():
        plt.figure(i)
        plt.xlabel('Energy [eV]')
        plt.ylabel('MOKE [deg]')
        plt.legend()
    if len(plt.get_fignums())>= 1:
        last = plt.get_fignums()[-1]
    samples = []
    if isinstance(loop_energy, list):
        for meas in loop:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            if sep == 0 or sep == 1:
                plt.figure(last+1)
            elif sep == 2:
                if meas.sample not in samples:
                    samples.append(meas.sample)
                plt.figure(last+samples.index(meas.sample))
            for energy in loop_energy:
                if isinstance(energy, list):
                    leg = legenda +'\nIntegralni smycka s energiemi '+str(round(energy[0],2))+' - '+str(round(energy[1],2)) + ' eV'
                    meas.data = meas.data[meas.data.index>energy[0]]
                    meas.data = meas.data[meas.data.index<energy[1]]
                    loop_moke = []
                    for column in meas.data.columns:
                        loop_moke.append(meas.data[column].mean())
                    pole = list(meas.data.columns)
                    pole = [float(i.split('T')[0]) for i in pole]
                    plt.plot(pole,loop_moke,label = leg)
                else:
                    rowenergie = meas.data.index.get_indexer([energy], 'nearest')
                    pole = list(meas.data.columns)
                    pole = [float(i.split('T')[0]) for i in pole]
                    leg = legenda +'\nEnergie '+str(round(meas.data.index[rowenergie[0]],2)) + ' eV'
                    plt.plot(pole,list(meas.data.iloc[rowenergie[0]]),label = leg)
    try:
        for i in plt.get_fignums()[plt.get_fignums().index(last)+1:]:
            plt.figure(i)
            plt.xlabel('Field [T]')
            plt.ylabel('MOKE [deg]')
            plt.legend()
    except:
        pass

def uka_vse(path = 'Not defined',num = 1,pm = False, loop_energy = 'Not defined'): 
    # Tato funkce ukaze vsechny mereni v adresari
    if path == 'Not defined':
        root = Tk()
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        initial_dir = os.getcwd()
        path = filedialog.askdirectory(initialdir = initial_dir,title = "Select measurement folders")
        root.destroy()
    Path_folder = pathlib.Path(path)
    filenames = list(Path_folder.rglob("*.dat")) 
    if num>1:
        for i in range(num-1):
            root = Tk()
            root.withdraw()
            root.call('wm', 'attributes', '.', '-topmost', True)
            initial_dir = os.getcwd()
            path = filedialog.askdirectory(initialdir = initial_dir,title = "Select measurement folders")
            root.destroy()
            Path_folder = pathlib.Path(path)
            filenames = filenames + list(Path_folder.rglob("*.dat")) 
    onepol,moke,loop,fm,pm = ukaz_LoadData(filenames)
    samples = []
    for meas in onepol:
        legenda = meas.sample
        if meas.temp != 'Not defined':
            legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
        if meas.repeat != 'Not defined':
            legenda = legenda + ' Repeat ' + str(meas.repeat)
        if meas.time != 'Not defined':
            legenda = legenda + ' Time ' + str(meas.time)
        if meas.field != 'Not defined':
            legenda = legenda + ' Field ' + str(meas.field)
        if meas.sample not in samples:
            samples.append(meas.sample)
        plt.figure(samples.index(meas.sample))
        plt.plot(meas.data['MOKE'], label = legenda)
        plt.title('One poalrity of field of '+ str(meas.sample))
    add = len(samples)
    samples = []
    for meas in moke:
        legenda = meas.sample
        if meas.temp != 'Not defined':
            legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
        if meas.repeat != 'Not defined':
            legenda = legenda + ' Repeat ' + str(meas.repeat)
        if meas.time != 'Not defined':
            legenda = legenda + ' Time ' + str(meas.time)
        if meas.sample not in samples:
            samples.append(meas.sample)
        plt.figure(add+samples.index(meas.sample))
        if pm:
            names = meas.data.columns
            for name in names:
                plt.plot(meas.data[name], label = legenda+' '+str(name))
        else:
            plt.plot(meas.data['MOKE'], label = legenda)
        plt.title('Standart MOKE measurement of '+ str(meas.sample))
    add = add + len(samples)
    samples = []
    for meas in fm:
        legenda = meas.sample
        if meas.temp != 'Not defined':
            legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
        if meas.sample not in samples:
            samples.append(meas.sample)
        plt.figure(add+samples.index(meas.sample))
        plt.title('Full MOKE measurement of '+ str(meas.sample))
        names = meas.data.columns
        for name in names:
            plt.plot(meas.data[name], label = legenda+' '+str(name))
    add = add + len(samples)
    samples = []
    for meas in pm:
        legenda = meas.sample
        if meas.temp != 'Not defined':
            legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
        if meas.sample not in samples:
            samples.append(meas.sample)
        plt.figure(add+samples.index(meas.sample))
        plt.title('Averaged MOKE measurement of '+ str(meas.sample))
        names = meas.data.columns
        for name in names:
            plt.plot(meas.data[name], label = legenda+' '+str(name))
    for meas in loop:
        legenda = meas.sample
        if meas.temp != 'Not defined':
            legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
        if meas.repeat != 'Not defined':
            legenda = legenda + ' Repeat ' + str(meas.repeat)
        if meas.time != 'Not defined':
            legenda = legenda + ' Time ' + str(meas.time)
        fig, ax = plt.subplots()

        colors = plt.cm.coolwarm(np.linspace(0, 1, len(meas.data.columns)))
        for i, col in enumerate(meas.data.columns):
            ax.plot(meas.data.index, meas.data[col], color=colors[i])
        cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.coolwarm), ax=ax, ticks=[0, 0.5, 1])
        cbar.ax.invert_yaxis()
        cbar.ax.set_yticklabels([meas.data.columns[0].rsplit('T')[0]+' T', meas.data.columns[len(meas.data.columns)//2].rsplit('T')[0]+' T', meas.data.columns[-1].rsplit('T')[0]+' T'])
        ax.lines[-1].set_label(legenda)
        ax.lines[-1].set_label(legenda)
        ax.set_title('Loop measurement of ' + meas.sample)
    for i in plt.get_fignums():
        plt.figure(i)
        plt.xlabel('Energy [eV]')
        plt.ylabel('MOKE [deg]')
        plt.legend()
    if len(plt.get_fignums())>= 1:
        last = plt.get_fignums()[-1]
    if isinstance(loop_energy, list):
        for meas in loop:
            legenda = meas.sample
            if meas.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
            if meas.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(meas.repeat)
            if meas.time != 'Not defined':
                legenda = legenda + ' Time ' + str(meas.time)
            plt.figure()
            for energy in loop_energy:
                if isinstance(energy, list):
                    leg = legenda +'\nIntegralni smycka s energiemi '+str(round(energy[0],2))+' - '+str(round(energy[1],2)) + ' eV'
                    meas.data = meas.data[meas.data.index>energy[0]]
                    meas.data = meas.data[meas.data.index<energy[1]]
                    loop_moke = []
                    for column in meas.data.columns:
                        loop_moke.append(meas.data[column].mean())
                    pole = list(meas.data.columns)
                    pole = [float(i.split('T')[0]) for i in pole]
                    plt.plot(pole,loop_moke,label = leg)
                else:
                    rowenergie = meas.data.index.get_indexer([energy], 'nearest')
                    pole = list(meas.data.columns)
                    pole = [float(i.split('T')[0]) for i in pole]
                    leg = legenda +'\nEnergie '+str(round(meas.data.index[rowenergie[0]],2)) + ' eV'
                    plt.plot(pole,list(meas.data.iloc[rowenergie[0]]),label = leg)
    try:
        for i in plt.get_fignums()[plt.get_fignums().index(last)+1:]:
            plt.figure(i)
            plt.xlabel('Field [T]')
            plt.ylabel('MOKE [deg]')
            plt.legend()
    except:
        pass

def elip(path = 'Not defined',show = True, save = True):
    if path == 'Not defined':
        initial_dir = os.getcwd()
    else:
        initial_dir = path
    root = Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    rotation = filedialog.askopenfilename(initialdir = initial_dir,title = "Select Kerr rotation",filetypes = (("txt files","*.dat .txt .KNT"),("all files","*.*")),multiple=False)
    root.destroy()
    root = Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    initial_dir = rotation.rsplit('/',1)[0]
    waveplates = filedialog.askopenfilename(initialdir = initial_dir,title = "Select waveplates",filetypes = (("txt files","*.dat .txt .KNT"),("all files","*.*")),multiple=False)
    root.destroy()
    fuck,rot_moke,this,fm,rot_pm = ukaz_LoadData([rotation])
    fuck,wav_moke,this,fm,wav_pm = ukaz_LoadData([waveplates])
    data =[]
    if len(rot_moke) == 1:
        data.append(rot_moke[0])
    elif len(rot_pm) == 1:
        data.append(rot_pm[0])
        data[0].data = data[0].data.rename(columns={'Average':'MOKE'})
    if len(wav_moke) == 1:
        data.append(wav_moke[0])
    elif len(wav_pm) == 1:
        data.append(wav_pm[0])
        data[1].data = data[1].data.rename(columns={'Average':'MOKE'})
    result = pd.concat([data[0].data['MOKE'],data[1].data['MOKE']],keys=['Rotation', 'Plates'],axis=1, join='inner')
    result = result.sort_index()
    result = result.interpolate(method = 'index')
    result = result[~result.index.duplicated(keep='first')] # drop duplicate indexes
    result = result[~result.isin([np.nan]).any(axis=1)] 
    energy = np.array(result.index.astype(float))
    lmbda = 1239.8/energy
    K = np.power(lmbda,2)
    W = 1- 93.0665**2/K
    Q = 1/(np.power(W,2) * lmbda * np.sqrt(1+136.24/W))
    delta = (Q * 1.69508759865 * 100000 + 2.884488929) * np.pi/180
    elipticita = - (np.array(result['Plates']) - np.array(result['Rotation'])*np.cos(delta))/np.sin(delta)
    result['Ellipticity'] = elipticita.tolist()
    result = result.drop(columns=['Plates'])
    if show:
        ax = plt.subplot()
        ax.plot(result['Rotation'], label = 'Rotace')
        ax.plot(result['Ellipticity'], label = 'Elipticita')
        plt.legend()
    if save:
        if data[0].sample == data[1].sample:
            sample = data[0].sample
        else:
            sample = sm.askstring("Confilcting sample names!!!","Enter a sample name:\t\t\t\t\t\nNames:  "+data[0].sample+'\n'+'           '+data[1].sample)
        path_save = rotation.rsplit('/',1)[0]+'\\'+str(sample)+'_Full_MOKE.dat'
        if os.path.exists(path_save):
            if askyesno(title='File exists', message='File already exists. Do you want to overwrite it?\nOne polarity measurement\n'+str(data[0])):
                os.remove(path_save)
            else:
                return
        file = open(path_save,'a')
        if data[0].field == data[1].field:
            file.write('# Type of measurement: Full MOKE +- ' + str(data[0].field) + 'T\n' )
        else:
            file.write('# Type of measurement: Full MOKE +- ?T\n' )
        
        daticka = ''
        if isinstance(data[0].date,list):
            for i in data[0].date:
                daticka = daticka + i + '\n#                 '
        else:
            daticka = data[0].date +    '\n#                 ' 
        if isinstance(data[1].date,list):
            for i in data[1].date:
                daticka = daticka + i  +    '\n#                 '
        else:
            daticka = daticka+data[1].date+ '\n#                 '
        daticka = daticka.rsplit('\n',1)[0]

        raw = ''
        if isinstance(data[0].path,list):
            for i in data[0].path:
                raw = raw + i + '\n#                 '
        else:
            raw = data[0].path +    '\n#                 ' 
        if isinstance(data[1].path,list):
            for i in data[1].path:
                raw = raw + i  +    '\n#                 '
        else:
            raw = raw+data[1].path+ '\n#                 '
        raw = raw.rsplit('\n',1)[0]

        if data[0].field == data[1].field:
            field = str(data[0].field)
        else:
            field = 'Not defined'
        if data[0].temp == data[1].temp:
            temp = str(data[0].temp)
        else:
            temp = 'Not defined'
        file.write('# Vzorek:         '+str(sample)+'\n# Field:          '+field+'\n# Temperature:    '+temp+'\n# Dates measured: '+daticka+'\n# OG data:        '+raw) 
        file.write('\n')                        
        file.close()
        result.to_csv(path_save,mode='a')

def prumer(path = 'Not defined', num = 1, show = True, save = True): 
    # Prumeruje mereni 
    if path == 'Not defined':
        initial_dir = os.getcwd()
    else:
        initial_dir = path
    filenames = []
    for i in range(num):
        root = Tk()
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        filenames = filenames + list(filedialog.askopenfilename(initialdir = initial_dir,title = "Select files",filetypes = (("txt files","*.dat .txt .KNT"),("all files","*.*")),multiple=True))
    root.destroy()
    onepol,moke,loop,fm,pm = ukaz_LoadData(filenames)
    if len(moke) > 1:
        data_to_ave = moke
        print('Averaging MOKE')
    elif len(onepol) > 1:
        data_to_ave = onepol
        print('Avereging One polarity of field')
    else:
        print('Error no MOKE or one polarity of field data to average')
        return
    fields = set([i.field for i in data_to_ave])
    temps = set([i.temp for i in data_to_ave])
    samples = set([i.sample for i in data_to_ave])
    times = set([i.time for i in data_to_ave])
    repeats = set([i.repeat for i in data_to_ave])
    if len(times) == len(data_to_ave):
        indexer = ['Time '+ str(i) for i in sorted(times)]
        data_to_ave.sort(key = lambda x: x.time)
    elif len(repeats) == len(data_to_ave):
        indexer = ['Repeat '+ str(i) for i in sorted(repeats)]
        data_to_ave.sort(key = lambda x: x.repeat)
    else:
        indexer = list(range(0,len(data_to_ave)))

    if len(fields) == 1:
        field = data_to_ave[0].field
    else:
        field = 'Not defined'
    if len(temps) == 1:
        temp = data_to_ave[0].temp
    else:
        temp = 'Not defined'
    if len(samples) == 1:
        sample = data_to_ave[0].sample
    else:
        sample='Not defined'
    if len(data_to_ave) != len(filenames):
        filenamess = []
        for cesta in data_to_ave:
            filenamess.append(cesta.path)
    else:
        filenamess = filenames
    raw = ''
    if isinstance(filenamess,list):
        for i in filenamess:
            raw = raw + i + '\n#                 '
    else:
        raw = ave.path +    '\n#                 ' 
    raw = raw.rsplit('\n',1)[0]
    if len(data_to_ave) != len(filenames):
        win = Tk()
        win.geometry("900x650")
        win.withdraw()
        not_filenames = ''
        for dataa in filenames:
            if dataa not in filenamess:
                not_filenames= not_filenames + dataa +'\n#                 ' 
        not_filenames = not_filenames.rsplit('\n',1)[0]
        messagebox.showerror('Filenames Error', 'Error: Not all data can be averaged!\n\nAveraged data: '+raw+'\n\n\nData not averaged: '+not_filenames)
    ave = Full_moke(data_to_ave[0].data,field,temp,sample,filenamess,[data_to_ave[0].date],'Prumer '+ data_to_ave[0].type.rsplit('+')[0].strip())
    ave.data = ave.data.rename(columns = {'MOKE':indexer[0]})
    ave.data = ave.data[[indexer[0]]]
    for file in data_to_ave[1:]:
        ave.data = pd.concat([ave.data,file.data['MOKE']],axis=1, join='inner')
        ave.data = ave.data.rename(columns = {'MOKE':indexer[data_to_ave.index(file)]})
        ave.date.append(file.date)
    ave.data = ave.data.sort_index()
    ave.data = ave.data.interpolate(method = 'index')
    ave.data = ave.data[~ave.data.index.duplicated(keep='first')] # drop duplicate indexes
    ave.data = ave.data[~ave.data.isin([np.nan]).any(axis=1)] 
    ave.data['Average'] = ave.data.mean(axis=1)
    if show:
        legenda = sample
        names = ave.data.columns
        for name in names:
            if name == 'Average':
                plt.plot(ave.data[name],'k--',label = legenda+' '+str(name))
            else:
                plt.plot(ave.data[name], label = legenda+' '+str(name))
        plt.title('Avereged MOKE measurements of '+str(sample))
        plt.legend()
        plt.show()
    if save:
        if sample == 'Not defined':
            sample = sm.askstring("Confilcting sample names!!!","Enter a sample name:\t\t\t")
        path_save = filenamess[0].rsplit('/',1)[0]+'\\'+str(sample)+'_Prumer_MOKE.dat'
        if os.path.exists(path_save):
            if askyesno(title='File exists', message='File already exists. Do you want to overwrite it?'):
                os.remove(path_save)
            else:
                return
        file = open(path_save,'a')
        if ave.field != 'Not defined':
            file.write('# Type of measurement: Prumer MOKE +- ' + str(ave.field) + 'T\n' )
        else:
            file.write('# Type of measurement: Prumer MOKE')
        daticka = ''
        if isinstance(ave.date,list):
            for i in ave.date:
                daticka = daticka + i + '\n#                 '
        else:
            daticka = ave.date +    '\n#                 ' 
        daticka = daticka.rsplit('\n',1)[0]
        file.write('# Vzorek:         '+sample+'\n# Field:          '+str(field)+'\n# Temperature:    '+str(temp)+'\n# Dates measured: '+daticka+'\n# Path OG data:   '+raw) 
        file.write('\n')                        
        file.close()
        ave.data.to_csv(path_save,mode='a')

def bg(path = 'Not defined', show = True, save = True):
    # odecte background smycky
    if path == 'Not defined':
        initial_dir = os.getcwd()
    else:
        initial_dir = path
    root = Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    data = filedialog.askopenfilename(initialdir = initial_dir,title = "Select Loop data",filetypes = (("txt files","*.dat"),("all files","*.*")),multiple=True)
    root.destroy()
    one_pol,moke,loop,fm,pm = ukaz_LoadData(data)
    for loo in loop:
        maxx = [i for i in loo.data.columns if abs(float(i.rsplit('T',1)[0])) == float(loo.field)]
        control = sum(list(np.sign([float(i.rsplit('T',1)[0]) for i in maxx])))
        if control != 0:
            win = Tk()
            win.geometry("900x650")
            win.withdraw()
            messagebox.showerror('Bg control error. Control = '+str(control), 'Error: Bg control error!\n\nControl = '+str(control)+'\n\nPay close attention to results!!')
        bg = loo.data[maxx].sum(axis=1)/len(maxx)
        loo.data = loo.data.subtract(bg, axis=0)
        loo.bg_sub = True
        if show:
            legenda = loo.sample+' Subtracted BG'
            if loo.temp != 'Not defined':
                legenda = legenda + ' Temperature ' + str(loo.temp) + 'K'
            if loo.repeat != 'Not defined':
                legenda = legenda + ' Repeat ' + str(loo.repeat)
            if loo.time != 'Not defined':
                legenda = legenda + ' Time ' + str(loo.time)
            fig, ax = plt.subplots()
    
            colors = plt.cm.coolwarm(np.linspace(0, 1, len(loo.data.columns)))
            for i, col in enumerate(loo.data.columns):
                ax.plot(loo.data.index, loo.data[col], color=colors[i])
            cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.coolwarm), ax=ax, ticks=[0, 0.5, 1])
            cbar.ax.invert_yaxis()
            cbar.ax.set_yticklabels([loo.data.columns[0].rsplit('T')[0]+' T', loo.data.columns[len(loo.data.columns)//2].rsplit('T')[0]+' T', loo.data.columns[-1].rsplit('T')[0]+' T'])
            ax.lines[-1].set_label(legenda)
            ax.lines[-1].set_label(legenda)
            ax.set_title('Loop measurement of ' + loo.sample + ' with BG subtraction')
        if save:
            path = loo.path.rsplit('.dat')[0]+'_No-BG.dat'
            if os.path.exists(path):
                if askyesno(title='File exists', message='File already exists. Do you want to overwrite it?\Loop MOKE measurement\n'+loo.__str__()):
                    os.remove(path)
                else:
                    return
            file = open(path,'a')
            file.write('# Type of measurement: Loop MOKE measurement - No BG\n' )
            file.write('# Date measured:       ' + loo.date +'\n')
            file.write(loo.__str__())
            file.write('\n')
            file.close()
            loo.data.to_csv(path,mode='a')
    if show:
        plt.legend()

def sava(path = 'Not defined',win = 51,pol=3, show = True, save = True): #TODO
    # Funkce na Savitzky-Golay filter
    pass

##########################################################################################################################
##########################################################################################################################
start = time.time()

cesta = r'C:\Users\tmale\OneDrive\Documents\Data\LSMO\Francie 2022\MOKE\PLD4150\loop' 


bg(cesta)




end = time.time()
print('\nExecution time: ',round(end-start,2),' seconds')
plt.show()



# Funkce ukafit: ukaze fit na jedne energii 
#              - energie : default 4 eV
#              - field : default 1
#   - Bude fungovat tak, ze si z nafitovaneho souboru precte path a parametry
#     a nasledne je nalezne ve slozce (snad) a nafituje na te konkretni energii
#
#
# Funkce savitzky: proste filtr
#                 - parameter save : default no
#
# Funkce vyplot findmin
#
# Funkce na view random Tim spektra - to by melo byt v uka!
#
# Funkce nacti data for advanced purposes - jen tkinter okno co ti dovoly vybrat data a ty to nacte
#
# Funkce tkinter select jenom proste