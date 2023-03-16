import os
from tkinter import filedialog
from tkinter import *
from tkinter.messagebox import *
import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

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
    def __init__(self,data,repeat,time,field,temp,sample,folder_path,date,type_meas):
        self.repeat = repeat          # Repeat
        self.time = time              # Time (repetice ne real time)
        self.field = field            # Pole
        self.temp = temp              # Teplota
        self.sample = sample          # sample name
        self.folder = folder_path     # path do slozky s merenimi
        self.data = data              # dataframe dat
        self.date = date              # Date measured
        self.type = type_meas         # Measurement type eg loop, moke,one polarity of field

    def __str__(self):
        return '# Vzorek       ' + str(self.sample)+'\n# Repeat       ' + str(self.repeat)+ '\n# Time         ' + str(self.time)+ '\n# Field        ' + str(self.field) + '\n# Temperature  ' + str(self.temp)+  '\n# Folder path  '+str(self.folder)
    

def help(funkce = None):  # TODO now fine, dont forget to update tho
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
    
    fce = {
        "fitni": fitni_help,
        "uka": uka_help,
        "uka_vse": uka_vse_help
        }
    if funkce == None:
        here_you_go = '\nV tomto programu existují následující funkce:\n'
        options = '''fitni\nuka\nuka_vse\n
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
    for file in files:
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

def fitni(show = True): # TODO average, show_fit
    data,path = LoadData()
    for i in range(0,len(data)):
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
    print('\nDone fitting!\n')
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
    data[0].fitted.rename(columns = {'Fit':'MOKE'})
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
    file.write(data[0].standart_moke_str())
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
    

def ukaz_LoadData(filenames): # TODO full MOKE
    # Nacte soubory jiz vyfitovanych dat
    one_pol = []
    moke = []
    loop = []
    for file in filenames:
        try:
            data = pd.read_csv(file, comment='#', index_col=0)
            f = open(file,'r')
            type_meas = f.readline().rsplit(': ')[1].strip()
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
            one_pol.append(measurement_result(data,repeat,time,field,temp,vzorek,folder,date,type_meas))
        elif 'Standart MOKE' in type_meas:
            moke.append(measurement_result(data,repeat,time,field,temp,vzorek,folder,date,type_meas))
        elif 'Loop MOKE' in type_meas:
            loop.append(measurement_result(data,repeat,time,field,temp,vzorek,folder,date,type_meas))
    return one_pol, moke, loop

def uka(path = 'Not defined',num = 1,sep = 1,pm = False, loop_energy = 'Not defined'): # TODO full MOKE # TODO hnusny barvy smycka
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
    onepol,moke,loop = ukaz_LoadData(filenames)
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
            ax.plot(meas.data)
            ax.lines[-1].set_label(legenda)
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
            ax.plot(meas.data)
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
        figures = []
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
            figures[samples.index(meas.sample)].plot(meas.data)
            ax.lines[-1].set_label(legenda)
            ax.set_title('Loop measurements of ' + meas.sample)
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
            ax.plot(meas.data)
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

def uka_vse(path = 'Not defined',num = 1,pm = False, loop_energy = 'Not defined'): # TODO full MOKE # TODO barvicky loops
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
    onepol,moke,loop = ukaz_LoadData(filenames)
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
    for meas in loop:
        legenda = meas.sample
        if meas.temp != 'Not defined':
            legenda = legenda + ' Temperature ' + str(meas.temp) + 'K'
        if meas.repeat != 'Not defined':
            legenda = legenda + ' Repeat ' + str(meas.repeat)
        if meas.time != 'Not defined':
            legenda = legenda + ' Time ' + str(meas.time)
        fig, ax = plt.subplots()
        ax.plot(meas.data)
        ax.lines[-1].set_label(legenda)
        ax.set_title('Loop measurements of ' + meas.sample)
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

def elip(path = 'Not defined',show = True, save = True): #TODO info hlavicka, TODO prumerovani kdyz vyberes vic souboru
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
    fuck,data,this = ukaz_LoadData([rotation,waveplates])
    result = pd.concat([data[0].data,data[1].data.rename(columns = {'MOKE':'Plates'})],axis=1, join='inner')
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
    elipticita = - (np.array(result['Plates']) - np.array(result['MOKE'])*np.cos(delta))/np.sin(delta)
    result['Elipticita'] = elipticita.tolist()
    result.drop(columns=['Plates'])
    if save:
        path_save = rotation.rsplit('/',1)[0]+'\\'+str(data[0].sample)+'_Full_MOKE.dat'
        if os.path.exists(path_save):
            if askyesno(title='File exists', message='File already exists. Do you want to overwrite it?\nOne polarity measurement\n'+str(data[0])):
                os.remove(path_save)
            else:
                return
        file = open(path_save,'a')
        file.write('# Type of measurement: Full MOKE +- ' + str(abs(data[0].field)) + 'T\n' )
        file.write(data[0].standart_moke_str()) # TODO tohle je potreba zadefinovat v classe, nezapomen na datum
        file.write('\n')
        file.close()
        result.to_csv(path_save,mode='a')
    if show:
        ax = plt.subplot()
        ax.plot(result['MOKE'], label = 'Elipticita')
        ax.plot(result['Elipticita'], label = 'Rotace')

##########################################################################################################################
##########################################################################################################################

cesta = r'C:\Users\tmale\OneDrive\Documents\Data\LSMO\Francie 2022\MOKE\PLD4150\loop\raw_data' 

start = time.time()

uka()

end = time.time()
print('Execution time: ',end-start,' seconds')

plt.show()

# Funkce prumeruj - zprumeruje vybrana mereni
#                  - ukaz : default false 

# Funkce ukafit: ukaze fit na jedne energii 
#              - energie : default 4 eV
#              - field : default 1
#   - Bude fungovat tak, ze si z nafitovaneho souboru precte path a parametry
#     a nasledne je nalezne ve slozce (snad) a nafituje na te konkretni energii
#
# Funkce elip: spocte elipticitu a vytvori jeden soubor rotace a elipticity
#              - soubor bude obsahovat data kdy bylo co zmerene, jmeno vzorku
#
# Funkce savitzky: proste filtr
#                 - parameter save : default no
#
# Funkce odecti BG smycka
#
# Funkce vyplot findmin
#
# Help