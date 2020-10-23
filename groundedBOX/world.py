#!/usr/bin/python3
"""
This file contains the necessary constants and the
Class 'World' which ecapsulates the information 
about the computation domain for the program :
====================================================
A 3D ElectroStatic PIC simulation
====================================================
The code simulates plasma confined in a reflecting box.
====================================================

"""

import numpy as np

#-----------------------------------------------------
#constants
EP0 = 8.85418782e-12 #C/V/m , permittivity of free space
QE = 1.602176565e-19 #C, electron charge
AMU = 1.660538921e-27 #kg , atomic mass unit
ME = 9.10938215e-31 #kg , electron mass
KB = 1.38e-23 #J/K , Boltzmann's constant

#-----------------------------------------------------
# the class encapsulating information about the computational domain
class World:
    
    #constructor
    def __init__(self,ni,nj,nk):
 
        self.nn= [ni, nj, nk] #number of nodes in x,y,z-direction
        
        #node volumes declaration
        self.node_vol = np.zeros((ni,nj,nk))
        
        #electric potential at nodes
        #also set Dirichlet-Boundary condition for grounded box
        self.phi = np.zeros((ni,nj,nk)) 
        
        #charge density at nodes
        self.rho = np.zeros((ni,nj,nk))
        
        #electric field at nodes
        self.efx = np.zeros((ni,nj,nk)) #x-direction
        self.efy = np.zeros((ni,nj,nk)) #y-direction        
        self.efz = np.zeros((ni,nj,nk)) #z-direction
        #self.ef = np.asarray([self.efx,self.efy,self.efz])
    #--- end __init__ ---------------------                   
    
    #sets the mesh span,computes cell-size, compute node-volumes
    def setExtents(self,x1,y1,z1,x2,y2,z2):
        #mesh origin, r_min
        self.r_min = np.asarray([x1,y1,z1])
        
        #diagonally opposite end, r_max
        self.r_max = np.asarray([x2,y2,z2])
        
        #self.dr=[] ; self.r_cen = []
        self.dr = np.zeros(3); self.r_cen = np.zeros(3)
        self.dr = (self.r_max - self.r_min)/(self.nn-np.ones(3))
        self.r_cen = 0.5*(self.r_max + self.r_min) 
        r = np.asarray(self.dr)
           
        #compute node-volumes
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]): #loop over nodes
                    dv = self.dr[0]*self.dr[1]*self.dr[2] #standard volume
                    if (i==0 or i==self.nn[0]-1): dv*=0.5
                    if (j==0 or j==self.nn[1]-1): dv*=0.5
                    if (k==0 or k==self.nn[2]-1): dv*=0.5
                    self.node_vol[i,j,k] = dv
    #---end setExtents-------------------
                    
    #set time-step and number of iterations
    def setTime(self,dt,num_dt):
        
        self.dt = dt #time-step 
        self.num_dt = num_dt #number of time-steps/iterations
        self.ts = 0 #iteration number
    #--- end setTime --------------------
    
    #advance time for looping 
    def advanceTime(self):
        self.ts += 1
        return self.ts <= self.num_dt    
    #--- end advanceTime -----------------
    
    #check if postion is within the computational domain
    def inBounds(self,pos):
        for i in range(3):
            if (pos[i]<self.r_min[i] or pos[i]>=self.r_max[i]):
                return False
        return True        
    #--- end inBounds --------------------
        
    #compute float-point-node,lc, from position x
    def XtoL(self,x):
        
        lc = (x - self.r_min)/self.dr
        
        return lc
    #--- end XtoL ------------------------
        
    #compute charge density at nodes, rho
    def computeChargeDensity(self,species):
        self.rho *= 0 #set charge density to zero
        
        for sp in species: #loop over species
            if (sp.charge==0):
                continue  #ignore neutral particles
        
            self.rho += sp.charge*np.ascontiguousarray(sp.den) #density scales by charge
    #--- end computeChargeDensity -----------    
           
                                                        
            
            
            
            
#------------------------------------------------------           
        
        
        
