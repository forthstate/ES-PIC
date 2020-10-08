"""
=================================
A 1D ElectroStatic PIC Simulation
=================================

This file contains visualisation code for the data
contained in the 'espic1d_data.csv', which is the output
from the 'espic1d.py'.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#---Load csv Data file

data = pd.read_csv('espic1d_data.csv')

#---Figure initialisation

fig, axs = plt.subplots(1,2,figsize=(20, 4))
  
#---position Animation in Subplot_0

dt = 1e-9 #time-step
niter = 5000 #number of iterations

  #set x-axis limits
axs[0].set_xlim(0,0.1)  
  #axes labels
axs[0].set_xlabel('Position(X), m')
  #drop y-ticks 
axs[0].set_yticks([])
  #title
axs[0].set_title('Position')
axs[0].grid()

line_x, = axs[0].plot([], [], 'o', lw=10) #particle position plot
time_template = 'time = %0.1fs ns' #time
time_text = axs[0].text(0.05, 0.9, '', transform=axs[0].transAxes)

  #init
def init_x():
    line_x.set_data([], [])
    time_text.set_text('')
    return line_x, time_text

  #animation function
def animate_x(i):
    x = data['X'].values[i]

    line_x.set_data(x,0)
    time_text.set_text(time_template % (1e9*data['Time'].values[i]))
    return line_x, time_text #, line_ke, line_pe   
    
anim_x = animation.FuncAnimation(fig, animate_x, interval=100, blit=True, init_func=init_x) 

#---------------------------------------------------------------------------------
#---energy statc-plot in Subplot_1

t_max = 200 #plot just for first t_max ns
  #axes labels
axs[1].set_xlabel('Time, ns')
axs[1].set_ylabel('Energy, eV')

  #axes limits for Energy plot
E_min = min(data['KE'].values.min(),data['PE'].values.min())
E_max = max(data['KE'].values.max(),data['PE'].values.max())
axs[1].set_ylim(E_min-1,E_max+1) #y-axis
axs[1].set_xlim(0,t_max) #x-axis

  #Plot
axs[1].plot(1e9*data['Time'].values[:t_max],data['KE'].values[:t_max],label='Kinetic Energy')
axs[1].plot(1e9*data['Time'].values[:t_max],data['PE'].values[:t_max],label='Potential Energy')
E_tot = data['KE'].values + data['PE'].values
axs[1].plot(1e9*data['Time'].values[:t_max],E_tot[:t_max],label='Total Energy')

axs[1].legend()
axs[1].grid()
axs[1].set_title('Energy vs Time')
#----------------------------------------------------------------------------------
"""#----------------------------------------------------------------------------------
###---energy Animation in Subplot_1

  #axes labels
axs[1].set_xlabel('Time, ns')
axs[1].set_ylabel('Energy, eV')

  #y-axis limits for Energy plot
E_min = min(data['KE'].values.min(),data['PE'].values.min())
E_max = max(data['KE'].values.max(),data['PE'].values.max())
axs[1].set_ylim(E_min-1,E_max+1)

line_ke, = axs[1].plot([], [],label='Kinetic Energy') #KE plot
line_pe, = axs[1].plot([], [],label='Potential Energy') #PE plot

axs[1].set_title('Energy vs Time')
axs[1].legend()
axs[1].grid()

  #init
def init_en():
    line_ke.set_data([], [])
    line_pe.set_data([], [])    
    return line_ke, line_pe

t=[] ; ke=[]; pe=[];
  #animation function
def animate_en(i):
    t.append(1e9*data['Time'].values[i])
    ke.append(data['KE'].values[i])
    pe.append(data['PE'].values[i])

    axs[1].set_xlim(i-10,i+4)
    line_ke.set_data(t,ke)
    line_pe.set_data(t,pe)

    return line_ke, line_pe
    
anim_en = animation.FuncAnimation(fig, animate_en, interval=100, blit=True,init_func=init_en) 
"""#-------------------------------------------------------------------------------------------                          
"""#-------------------------------------------------------------------------------------------
#---Velocity static plot in Subplot_1

t_max = 100 #plot for first t_max ns

  #axes labels
axs[1].set_xlabel('Time, ns')
axs[1].set_ylabel('Velocity, m/s')

  #axes-limits
v_min = data['V'].values.min()
v_max = data['V'].values.max()
axs[1].set_ylim(v_min-(v_max/10),v_max+(v_max/10)) #y-limits
axs[1].set_xlim(0,t_max) #x-limits
  #Plot
axs[1].plot(1e9*data['Time'].values[:t_max],data['V'].values[:t_max],label='Velocity')

  #title
axs[1].set_title('Velocity vs Time')
axs[1].grid()
"""#-----------------------------------------------------------------------------------------------


plt.show()
 
 





