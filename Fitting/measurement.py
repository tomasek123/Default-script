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

class measurement: #TODO pridej parameter sample name a measurement folder
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
        return '# Vzorek       ' + str(self.sample)+'\n# Repeat       ' + str(self.repeat)+ '\n# Time         ' + str(self.time)+ '\n# Field        ' + str(abs(self.field)) + '\n# Temperature  ' + str(self.temp)+ '\n# Folder path  '+str(self.folder)

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


def LoadData(): 
    # Load all txt files in chosen directory
    root = Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    initial_dir = os.getcwd()
    initial_dir = r'C:\Users\tmale\OneDrive\Documents\Data\LSMO\Francie 2022\MOKE\PLD4150\loop\raw_data'  #TODO delete
    folder = filedialog.askdirectory(initialdir = initial_dir,title = "Select measurement folders")
    root.destroy()
    Path_folder = pathlib.Path(folder)
    files = list(Path_folder.rglob("*.txt"))

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

def fitni():
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

    for mereni in splacnuto:
        if len(mereni) == 1:
            save_one_polarity(mereni,path)
        elif len(mereni) == 2:
            save_kerr(mereni,path)
        else:
            save_smycka(mereni,path)

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
    file.write(str(data[0]))
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

def save_smycka(data,path):# TODO
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

    # Jmena sloupcu, musi byt jina?? Vyzkousej... TODO
    data.sort(key = lambda x: int(x.list[0].RT.replace('_','')))
    data_save = data[0].fitted[['Fit']]
    data_save = data_save.rename(columns = {'Fit': str(data[0].field) + 'T'})
    name_colum = [str(data[0].field) + 'T']
    for i in range(1,len(data)):
        name_col_single = str(data[i].field) + 'T'
        j = 2
        while name_col_single in name_colum:
            name_col_single = str(data[i].field) + 'T ' + str(j)
            j = j + 1
        name_colum.append(name_col_single)
        data[i].fitted = data[i].fitted.rename(columns = {'Fit': name_col_single})
        data_save = pd.concat([data_save,data[i].fitted[[name_col_single]]],axis=1, join='inner')

    data_save = data_save.sort_index()
    data_save = data_save.interpolate(method = 'index')
    data_save = data_save[~data_save.index.duplicated(keep='first')] # drop duplicate indexes
    data_save = data_save[~data_save.isin([np.nan]).any(axis=1)] 
    data_save.to_csv(path,mode='a')
    

start = time.time()

fitni()

end = time.time()

print('Execution time: ',end-start,' seconds')


# Funkce fitni
#               - nasledne ukaz : false
#               - nasledne prumeruj : false     

# Funkce prumeruj

# Funkce ukaz: - pm: default false
#              - energie konkretni na smycku : default None

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