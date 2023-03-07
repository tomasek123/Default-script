import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


data = pd.read_csv(r'C:\Users\tmale\OneDrive\Documents\Data\LSMO\Francie 2021\Elipsometrie\PLD3978\PLD3978_19_09_22_opticke_konst.txt',index_col=0)

data['n']=(data['e1r']+1j*data['e1i'])**(1/2)

theta = np.deg2rad(20)

data['thetaT'] = np.arcsin(np.sin(theta)/data['n'])

data['rs'] = np.abs(np.real((np.cos(theta)-data['n']*np.cos(data['thetaT']))/(np.cos(theta)+data['n']*np.cos(data['thetaT']))))
data['rp'] = np.abs(np.real((data['n']*np.cos(theta)-np.cos(data['thetaT']))/(data['n']*np.cos(theta)+np.cos(data['thetaT']))))


plt.plot(data['rs'])


plt.ylim(0,1)
plt.show()