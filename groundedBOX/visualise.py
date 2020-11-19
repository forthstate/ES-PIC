#!/usr/bin/python3
"""
This file contains the code to visualise
the diagsnostics output 
of the program :
====================================================
A 3D ElectroStatic PIC simulation
====================================================
The code simulates plasma confined in a reflecting box.
====================================================
NOTE : choose between normalised and regular energy plots
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

#---open and read file----
filename = "diagnostics.csv"
data = pd.read_csv(filename)

#---inistialse figure----
fig,axs = plt.subplots(1,2,figsize=(20,10))

"""
#--- normalised energy plot----
#normalised data to plot
ts = data['Ts'].values
#kinetic energy
ke_o = data['KE.O+'].values
ke_o /= ke_o.max()
ke_e = data['KE.e-'].values
ke_e /= ke_e.max()
#potential energy
pe = data['PE'].values
pe /= pe.max()
#total energy
te = data['TE'].values
te /= te.max()

axs[0].plot(ts,ke_o,label='KE.O+')
axs[0].plot(ts,ke_e,label='KE.e-')
axs[0].plot(ts,pe,label='PE')
axs[0].plot(ts,te,label='total_Energy')

axs[0].set_xlabel('Ts * 2e-10s')
axs[0].set_ylabel('normalised Energy')
axs[0].set_title('Energy vs Time')

axs[0].grid()
axs[0].legend()
"""
#---regular energy plot----
axs[0].plot(data['Ts'].values,data['KE.O+'].values,label='KE.O+')
axs[0].plot(data['Ts'].values,data['KE.e-'].values,label='KE.e-')
axs[0].plot(data['Ts'].values,data['PE'].values,'--',label='PE')
axs[0].plot(data['Ts'].values,data['TE'].values,label='total_Energy')

axs[0].set_xlabel('Ts * 2e-10s')
axs[0].set_ylabel('Energy / J')
axs[0].set_title('Energy vs Time')

axs[0].grid()
axs[0].legend()

#---momentum plot----
axs[1].plot(data['Ts'].values,data['px.e-'].values,label='px.e-')
axs[1].plot(data['Ts'].values,data['py.e-'].values,label='py.e-')
axs[1].plot(data['Ts'].values,data['pz.e-'].values,label='pz.e-')

axs[1].set_xlabel('Ts * 2e-10s')
axs[1].set_ylabel('Momentum / kg-m/s')
axs[1].set_title('Momentum vs Time')

axs[1].grid()
axs[1].legend()

plt.show()


