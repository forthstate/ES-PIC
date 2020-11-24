#!/usr/bin/python3
"""
still in Trial phase

====================================================
A 3D ElectroStatic PIC simulation
====================================================
The code simulates plasma confined in a reflecting box.
====================================================
"""

import numpy as np
import time
import csv
from world import *
from fieldSolver import *
from species import *
from output import *

start = time.time()

world = World(21,21,21)
world.setExtents(0,0,0,0.1,0.1,0.1)
world.setTime(2e-10,2000)

species = []
species.append(Species('O+',16*AMU,QE,world))
species.append(Species('e-',ME,-1*QE,world))

#specify particles
np_ions = 60000
np_electrons = 8000
num_den = 1e11 #number density

#load Particles
species[0].loadParticlesBox(world,world.r_min,world.r_max,num_den,[41,41,41],EqualSpace=True)
species[1].loadParticlesBox(world,world.r_min,world.r_cen,num_den,[21,21,21],EqualSpace=True)

for sp in species:
    print(sp.name,' has ',len(sp.particles),' particles.')

t_s = time.time()-start
print('time taken for species allotment: ',t_s)

#initialise diagnostics output file
output_diag_init(species)

#solve initial Potential
computePotential(world,10000,1e-4)
t_pot = time.time()-start
print('time taken for potential: ',t_pot-t_s,'\nphi[6,7,19]: ',world.phi[6,7,19])
#get initial Electric Field
computeEF(world)
t_ef = time.time()-start
print('time taken for electricField: ',t_ef-t_pot)

#main loop for advancing particles
print('loop started!')
while(world.advanceTime()):
        
    t1=time.time()
        #move particles
    for sp in species:
        sp.advance(world)
        sp.computeNumDen(world)
    
    t2=time.time()
    
        #compute charge density
    world.computeChargeDensity(species)
    
    t3=time.time()

        #update Potential
    computePotential(world,10000,1e-4)
    
    t4=time.time()
    
        #update Electric Field
    computeEF(world)  
    
    t5=time.time() 
    
    #diagnostics to file
    if (world.ts==1)or(world.ts%250==0):
        output_diag(world,species)
        t6 = time.time()
        print('Time: ',world.ts,'time taken for Energy+Momentum: ',t6-t5,'\n')     
        
    #runtime diagnostics      
    if (world.ts==1)or(world.ts%100==0):
        print('Time: ',world.ts,' PHI[10,11,15]: ',world.phi[10,11,15])
        print('O+,',' pos: ',species[0].particles[20000].pos)
        print('e-,',' pos: ',species[1].particles[2000].pos)
        print('time taken for ADVANCE+NUM_DEN: ',t2-t1)
        print('time taken for CHARGE_DEN: ',t3-t2)
        print('time taken for POTENTIAL: ',t4-t3)
        print('time taken for ELECTRIC_FIELD: ',t5-t4)
        print('time taken per loop: ',time.time()-t1)
        print('TOTAL time taken since START: ',time.time()-start,'\n')
        #output data to a VTK file to visualisation
        outputVTK(world.ts,world,species)
       
print('\nTotal time taken: ',time.time()-start,'\nAll is Well!!')



