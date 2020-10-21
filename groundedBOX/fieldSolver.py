#!/usr/bin/python3
"""
This file ecapsulates the methods necessary to solve for 
various field values - Potential, Electric Field for 
the program:
====================================================
A 3D ElectroStatic PIC simulation
====================================================
The code simulates plasma confined in a reflecting box.
====================================================

"""

import numpy as np
from world import *
import utils
#--------------------------------------------------------------------
#function to compute Potential
def computePotential(world,max_iter,tol):
    """
    comptutes potential at the nodes of the computational
    domain by solving the PoissonsEquation using GaussSiedel 
    method with -
    Input : world - a World object
            max_iter - maximum number of iterations
            tol - tolerance value to check convergence
    Return : boolean convergence check
    
    """
    w = 1.4 #SOR weight 
    l2 = 0; #least-square-error
    converged = False #boolean value to return
    
    #ep0 = world.EP0
    #precompute 1/dx^2 , 1/dy^2, 1/dz^2
    idx2 = 1/(world.dr[0]**2)
    idy2 = 1/(world.dr[1]**2)
    idz2 = 1/(world.dr[2]**2)        
    
    rho = np.asfortranarray(world.rho)
    phi = np.asfortranarray(world.phi)
        
    for iter in range(max_iter):
        
        phi = utils.potential_loop(world.rho,world.phi,idx2,idy2,idz2,EP0,w)
                    
        #convergence check
        if (iter%25 == 0) :
        
            r2sum = utils.convergence_loop(world.rho,world.phi,idx2,idy2,idz2,EP0)
            
            l2 = np.sqrt(r2sum/(world.nn[0]*world.nn[1]*world.nn[2]))
            if (l2 < tol) :
                converged = True
                break
                
    if (converged==False):
        print('GaussSiedel failed to converge, L2 = ',l2)
    
    world.phi = np.ascontiguousarray(phi)
    return converged        

#-------------------------------------------------------------------       
#function to compute Electric Field = -gradient(phi)
def computeEF(world):
    """
    computer electric field to by finding the derivative of the
    Potential using second order finite-difference-methods with:
    INPUT : world object
    """
    phi = np.asfortranarray(world.phi)
    efx = np.asfortranarray(world.efx)
    efy = np.asfortranarray(world.efy)
    efz = np.asfortranarray(world.efz)
    
    efx,efy,efz = utils.electric_field(world.phi,efx,efy,efz,world.dr[0],world.dr[1],world.dr[2])
    
    world.efx = np.ascontiguousarray(efx)
    world.efy = np.ascontiguousarray(efy)
    world.efz = np.ascontiguousarray(efz)
    
#-------------------------------------------------------------



