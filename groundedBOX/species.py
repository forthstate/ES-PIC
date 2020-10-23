#!/usr/bin/python3
"""
This file contains the necessary classes and methods that
encapsulate all the information about the charged/neutral particles 
that populate the computational domain
Class 'Species' - classifies type of particle
Class 'Particle' - properties of every particle
for the program :
====================================================
A 3D ElectroStatic PIC simulation
====================================================
The code simulates plasma confined in a reflecting box.
====================================================
"""
import numpy as np
import random
from world import *
import data_utils

#-----------------------------------------------------------------
#class of every particle's properties
class Particle:

    #constructor
    def __init__(self,pos,vel,mpw):
        self.pos = pos #position vector
        self.vel = vel #velocity vector
        self.mpw = mpw #macro-particle-weight
        
#-----------------------------------------------------------------
#class to group particles by gas species
class Species:
    
    #constructor
    def __init__(self,name,mass,charge,world):
        self.name = name #name of species
        self.mass = mass #mass of species-particle 
        self.charge = charge #charge of species-particle
        
        #number density
        self.den = np.zeros((world.nn[0],world.nn[1],world.nn[2]))
        self.den = np.asfortranarray(self.den)
        
        #list of particles of class Particles of the species
        self.particles = []
    
    #--end __init__----------
        
                    
    #add a new particle to the particles list
    def addParticle(self,pos,vel,mpw,world,efx,efy,efz):
        
        #boundary check
        if (world.inBounds(pos)==False): return
        
        #get logical coordinate
        lc = world.XtoL(pos)
        
        #compute Electric field at particle position   
        #interpolate EF data from 8 surrounding nodes to particle
        efx_par,efy_par,efz_par = data_utils.gatheref(lc,efx,efy,efz) 
                
        #rewind velocity by 0.5dt for the LEAPFROG METHOD
        dt = world.dt
        vel[0] -= efx_par*(0.5*dt*self.charge/self.mass)
        vel[1] -= efy_par*(0.5*dt*self.charge/self.mass)            
        vel[2] -= efz_par*(0.5*dt*self.charge/self.mass)
        
        self.particles.append(Particle(pos,vel,mpw))
    
    #--end addParticle----------
    
    #loads num_mp particles in the BOX
    def loadParticlesBox(self,world,x1,x2,num_den,num_mp,EqualSpace=False):
        """
        INPUT : r1, r2 - extent of box
                num_den - given number density
                num_mp - number of macro-particles
        """
        #box volume
        box_vol = (x2[0] - x1[0])*(x2[1] - x1[1])*(x2[2] - x1[2])
        #real number of particles
        num_real = num_den*box_vol
        
        #EF values as fortran arrays to add particle
        efx = np.asfortranarray(world.efx)
        efy = np.asfortranarray(world.efy)
        efz = np.asfortranarray(world.efz)
        
        if EqualSpace==True :
            #load particles in Equally Spaced Intervals
            
            num_mp_tot = (num_mp[0]-1)*(num_mp[1]-1)*(num_mp[2]-1)
            #macro-particle-weight
            mpw = num_real/num_mp_tot
            
            #compute particle grid spacing
            dp = []
            for i in range(3):
                dp.append((x2[i]-x1[i])/(num_mp[i]-1))
                
            #load particles 
            for i in range(num_mp[0]):
                for j in range(num_mp[1]):
                    for k in range(num_mp[2]):
                       
                        pos = np.zeros(3)
                        pos[0] = x1[0] + i*dp[0]
                        pos[1] = x1[1] + j*dp[1]
                        pos[2] = x1[2] + k*dp[2]
                        
                        #shift particles on boundary back into domain
                        if (pos[0]==x2[0]): pos[0] -= 1e-4*dp[0]
                        if (pos[1]==x2[1]): pos[1] -= 1e-4*dp[1]                        
                        if (pos[2]==x2[2]): pos[2] -= 1e-4*dp[2]
                        
                        #relative-weight for boundary particles
                        #faces:0.5w, edges:0.25w, corners:0.125w
                        w = 1
                        if (i==0 or i==num_mp[0]-1): w*0.5
                        if (j==0 or j==num_mp[1]-1): w*0.5                
                        if (k==0 or k==num_mp[2]-1): w*0.5
                        
                        #initial velocity
                        vel = np.zeros(3) #stationary
                        
                        #add the new particle to the array
                        self.addParticle(pos,vel,mpw*w,world,efx,efy,efz)                        
            
        else :
            #load particles in UNIFORM RANDOM Distribution
            
            #macro-particle-weight
            mpw = num_real/num_mp
            
             
            for i in range(num_mp):
                #add by sampling position
                pos = np.zeros(3)
                pos[0] = x1[0] + random.uniform(0,1)*(x2[0]-x1[0])
                pos[1] = x1[1] + random.uniform(0,1)*(x2[1]-x1[1])            
                pos[2] = x1[2] + random.uniform(0,1)*(x2[2]-x1[2])
            
                #initial velocity
                vel = np.zeros(3) #stationary
            
                #add the new particle to the array
                self.addParticle(pos,vel,mpw,world,efx,efy,efz)
    
    #--end loadParticlesBox--------        
        
    #compute number density
    def computeNumDen(self,world):
        self.den *= 0 #set density to zero
        
        #loop over particles
        for par in self.particles:
            lc = world.XtoL(par.pos)
            
            #SCATTER function: particle-data to mesh interpolation
            """
            #make sure to be in domain
            if (lc[0]<0 or lc[0]>world.nn[0]-1 or lc[1]<0 or \
                lc[1]>world.nn[1]-1 or lc[2]<0 or lc[2]>world.nn[2]-1):
                return
            """
            self.den = data_utils.scattermpw(lc,self.den,par.mpw)
            
        self.den /= np.asfortranarray(world.node_vol)
            
    #--end computeNumDen--------
    
    #move particles using electric field, ef
    def advance(self,world):
        
        #time-step
        dt = world.dt
        
        #EF values as fortran arrays to add particle
        efx = np.asfortranarray(world.efx)
        efy = np.asfortranarray(world.efy)
        efz = np.asfortranarray(world.efz)
        
        #pre-Compute
        dcbym=dt*self.charge/self.mass
        #loop over particles
        for par in self.particles:
            
            #get logical coordinate
            lc = world.XtoL(par.pos)
            
            #computer electric field at particle position
            #interpolate EF data from 8 surrounding nodes to particle
            efx_par,efy_par,efz_par = data_utils.gatheref(lc,efx,efy,efz)
            
            #update velocity using, a=qE/m
            par.vel[0] += efx_par*dcbym#(dt*self.charge/self.mass)
            par.vel[1] += efy_par*dcbym#(dt*self.charge/self.mass)            
            par.vel[2] += efz_par*dcbym#(dt*self.charge/self.mass)
            
            #update position, x=v*dt
            par.pos += par.vel*dt    
            
            #reflect particles at walls back into domain
            for i in range(3):
                if (par.pos[i] < world.r_min[i]) :
                    par.pos[i] = 2*world.r_min[i] - par.pos[i]
                    par.vel[i] *= -1.0 #reverse velocity
                elif (par.pos[i] >= world.r_max[i]):
                    par.pos[i] = 2*world.r_max[i] - par.pos[i]
                    par.vel[i] *= -1.0 #reverse velocity 
        
        #--end advance----------
        
        
        
        
        
        
