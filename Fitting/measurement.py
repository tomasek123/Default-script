import os
from tkinter import filedialog
from tkinter import *
import pathlib
import pandas as pd
import numpy as np

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
    def __init__(self,spec,repeat,time,field,temp,findmin):
        self.list = spec      # List spekter v tomto mereni
        self.repeat = repeat
        self.time = time
        self.field = field
        self.temp = temp
        self.findmin = findmin

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
        data = data[~data.isin([np.nan]).any(1)] 

        return data
    
    def is_fitable(self):
        is_fit = True
        for meas in self.list:
            if meas.angle == 'Not defined':
                is_fit = False
        if self.field == 'Not defined':
            is_fit = False
        if self.findmin:
            is_fit = False
        
        return is_fit
    
    def fit(self): # todo dodelat
        if self.is_fitable():
            lel = True
        else:
            lel = 'This is not a MOKE measurement!'
        return lel

def LoadData():
    # Load all txt files in chosen directory
    root = Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    initial_dir = os.getcwd()
    initial_dir = r'C:\Users\tmale\Downloads\Test'  #todo delete
    folder = filedialog.askdirectory(initialdir = initial_dir,title = "Select measurement folders")
    root.destroy()
    Path_folder = pathlib.Path(folder)
    files = list(Path_folder.rglob("*.txt"))

    # Initialize a list of all spectra and of all measurement configurations
    # spectra_list is a 2D list, for each measurement it contains a list of all spectra
    spectra_list = []
    measurements = []

    # for each path to a txt measurement file disect its name into components
    # and load it into the spectra class
    for file in files:
        name = str(file).rsplit('\\',1)[1]
        data = pd.read_csv(file, comment='#', index_col=0, sep='\t', header=None)
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
        meas = [properties[0],properties[1],properties[3],properties[4],findmin]
        if meas in measurements:
            indexx = measurements.index(meas)
            spectra_list[indexx].append(spectra(data,file,RT,properties[0],properties[1],properties[2],properties[3],properties[4],properties[5]))
        else:
            measurements.append(meas)
            spectra_list.append([spectra(data,file,RT,properties[0],properties[1],properties[2],properties[3],properties[4],properties[5])])

    for i in range(0,len(measurements)):
        measurements[i] = measurement(spectra_list[i],measurements[i][0],measurements[i][1],measurements[i][2],measurements[i][3],measurements[i][4])
 
    # Return a list of measurement class objects 
    return measurements


result = LoadData()

for i in result:
    print(i.fit())

