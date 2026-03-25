import pandas as pd
import matplotlib.pyplot as plt
# from package_DBR import myRound, FOPDT, FOPDT_cost, SOPDT, SOPDT_cost
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))

parent_dir = os.path.dirname(current_dir)

if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)
# etant donné que package dbr est dans le dossier parent on peut pas directement import
from package_DBR import *


nameFile = 'Cleaned_data_Open_loop_experiment_on_DV_2026-03-04-18h52.txt'

  
titleName = nameFile.split('.')[0]    
data = pd.read_csv('Data/' + nameFile)

if 'MV' in nameFile:
    ExpVariable = 'MV'
    tm = data['tm'].values
    MVm = data['MVm'].values
    PVm = data['PVm'].values    
else:    
    ExpVariable = 'DV'
    tm = data['tm'].values
    DVm = data['DVm'].values 
    PVm = data['PVm'].values
     


#SOPDT DV
K_DV_SOPDT = 0.393368177568377
T1_DV_SOPDT = 161.75194144935833
T2_DV_SOPDT = 36.28794803765712
theta_DV_SOPDT = 9.570305019620053

#FOPDT DV

K_DV_FOPDT = 0.3956172303278946
T_DV_FOPDT = 174.43578111991812
theta_DV_FOPDT = 38.3308784234789

Ts = 1



PVSim_FOPDT = FOPDT(DVm,K_DV_FOPDT,T_DV_FOPDT,theta_DV_FOPDT,Ts)    
PVSim_SOPDT = SOPDT(DVm,K_DV_SOPDT,T1_DV_SOPDT, T2_DV_SOPDT,theta_DV_SOPDT,Ts)    


plt.figure(figsize = (22,12))

plt.subplot(2,1,1)
if ExpVariable == 'MV':
    plt.step(tm,MVm,'b-',linewidth=2,label='Cleaned MV',where='post')
    plt.ylabel('Value of MV')
else:
    plt.step(tm,DVm,'b-',linewidth=2,label='Cleaned DV',where='post')
    plt.ylabel('Value of DV') 
    
plt.legend(loc='best')
plt.xlim([0, tm[-1]])

plt.subplot(2,1,2)
plt.step(tm,PVm,'g-',linewidth=2,label='Cleaned experimental data',where='post')
plt.step(tm,PVSim_FOPDT,'r-',linewidth=2,label='Identified first order plus delay response',where='post')
plt.step(tm,PVSim_SOPDT,'blue',linewidth=2,label='Identified first order plus delay response',where='post')
plt.ylabel('Value of PV')
plt.xlabel('Time [s]')
plt.legend(loc='best')
plt.xlim([0, tm[-1]])



plt.step(tm,PVSim_FOPDT,'r-',linewidth=2,label='Identified first order plus delay response',where='post')
plt.title('Experiment and identified first order plus delay response')

plt.step(tm,PVSim_SOPDT,'blue',linewidth=2,label='Identified second order plus delay response',where='post')
plt.title('Experiment and identified second order plus delay response')

plt.show()