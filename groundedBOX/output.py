#!/usr/bin/python3
"""
This file contains the Output fucntions 
to create data files to run diagnostics
for the program :
====================================================
A 3D ElectroStatic PIC simulation
====================================================
The code simulates plasma confined in a reflecting box.
====================================================
"""
from pyevtk.hl import imageToVTK
import numpy as np
import csv
from world import *

#---fucntion to write to vti/vtk file----------------- 
def outputVTK(ts,world,species):
    filename = "results/fields_"+str(ts)
    
    origin = (world.r_min[0],world.r_min[1],world.r_min[2])
    spacing = (world.dr[0],world.dr[1],world.dr[2])
    
    #variables
    ndO = species[0].den
    nde = species[1].den

    imageToVTK(filename, origin=origin, spacing=spacing, cellData = None, pointData={"NodeVol":world.node_vol, "Rho":world.rho, \
              "Phi":world.phi, "EFx":world.efx, "EFy":world.efy, "EFz":world.efz, "nd.O+":ndO, "nd.e-":nde })
              
    #--- end outputVTK-----------------
    
#---function to initialise diagnostics file------
def output_diag_init(species):
    #open file to write
    filename = "diagnostics.csv"
    csvfile = open(filename,'w')
    file = csv.writer(csvfile)
    #set labels
    labels=['Ts','PE'] 
    for sp in species:
        labels.append('KE.'+str(sp.name))
        labels.append('px.'+str(sp.name))
        labels.append('py.'+str(sp.name))
        labels.append('pz.'+str(sp.name))
    labels.append('TE')
    
    #write labels to file
    file.writerow(labels)
    
    #close file
    csvfile.close()
    
    #--- end output_diag_init--------
    
#---function to write diagnostics to file---------
def output_diag(world,species):
    #open file to write
    filename = "diagnostics.csv"    
    csvfile = open(filename,'a')
    file = csv.writer(csvfile)
    
    #diagnostics data
    diags=[world.ts] #time
    e = world.getPE() #PotentialEnergy
    diags.append(e)
    for sp in species:
        ke = sp.getKE() #kineticEnergy
        diags.append(ke)
        e += ke
        mom = sp.getMOM() #Momentum
        diags.append(mom[0])
        diags.append(mom[1])
        diags.append(mom[2])
    diags.append(e) #totalEnergy
    
    #write data to file
    file.writerow(diags)
    
    #close file
    csvfile.close()
    
    #--- end output_diag -----------
        
        
    
        
             
              
