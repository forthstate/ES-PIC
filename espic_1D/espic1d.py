"""
=================================
A 1D ElectroStatic PIC Simulation
=================================

Simulation of an test-particle in a 1D-potential well containing 
a plasma consisting of singly-positively-charged ions of 
uniform and constant density. The walls of infinite depth and 
height are kept at ZERO-potential. This model ignores any 
collision losses as well as any perturbance in the potential
created by the single test-particle. Also B=0.
NOTE: Check the "visualise.py" file for the corresponding plots 
      and animations.
"""
import numpy as np
from matplotlib import pyplot as plt
import csv
import time 

start = time.time() #start point

#---constants 
QE = 1.602e-19  # C , electron charge
EP0 = 8.85e-12 # C/V/m , permittivity of free space
ME = 9.1e-31 # kg , mass of en electron

NI = 1e12 # m^-3 , given ion number density

#---spatial grid
n = 21 #number of nodes
x0 = 0 #start-node(origin)
xl = 0.1 #end-node 
dx = (xl-x0)/(n-1) #step-size
dx2 = dx*dx

#---fields initialization
rho = QE*NI*np.ones(n) #charge density, constant and uniform 
phi = np.zeros(n) #potential with intial values: zero
ef = np.zeros(n) #electric field with intial values: zero


#---Potential field
"""
Solves the 1D-Poissons equation to compute the Potential profile 
along grid using finite-central-difference-method:
d2/dx2(phi) = (phi[i-1] - 2*phi[i] + phi[i+1])/dx2
TriDiagonal algorithm is used for solving the subsequent matrix.
"""
poisson_time = time.time() #time stamp to analyse poisson solver

  #coefficients and RHS
a = np.zeros(n) 
b = np.zeros(n)
c = np.zeros(n)
d = np.zeros(n) 
  #set coefficients
a[1:n-1] = 1/dx2 #coefficient of phi[i-1]
b[1:n-1] = -2/dx2 #coefficient of phi[i]
c[1:n-1] = 1/dx2 #coefficient of phi[i+1]
d[1:n-1] = -rho[1:n-1]/EP0 #RHS
b[0]=b[n-1]=1 ; d[0]=d[n-1]=0 #Dirichlet boundary condition
  
  #forward substitution
c[0] = c[0]/b[0] ; d[0] = d[0]/b[0]
for i in range(1,n):
    if i < n-1:
        c[i] = c[i]/(b[i] - a[i]*c[i-1])
    d[i] = (d[i] -a[i]*d[i-1])/(b[i] - a[i]*c[i-1])

  #backward substitution
phi[n-1] = d[n-1]
for i in range(n-2,-1,-1):
    phi[i] = d[i] - c[i]*phi[i+1]

print('Time taken by the Poisson solver:', time.time()-poisson_time)

  #maximum potential 
phi_max = phi.max() #maximum potential value
  

#---Electric field
""" 
computes the electric field by differentiating the Potential: 
ef = -d/dx(phi)
"""
  #with finite central difference for internal nodes
for i in range(1,n-1):
    ef[i] = -(phi[i+1] - phi[i-1])/(2*dx)

  #second order finite difference for boundary nodes
ef[0] = (3*phi[0] - 4*phi[1] + phi[2])/(2*dx) #start-node
ef[n-1] = (-3*phi[n-1] + 4*phi[n-2] - phi[n-3])/(2*dx) #end-node


#---test particle motion

niter = 5000 #number of iterations
dt = 1e-9 #time-step
  
  #load particle
q = -QE #charge of test_particle
m = ME #mass of the test_particle  
x = 1.5*dx  #initial position
v = 0  #initial velocity
  
  #potential/electric field interpolation
def gather(x,field):
    """ 
    The function utilises linear-gather-method to find the 
    electric-field/potential value acting on the particle 
    corresponding to its positon with -
    input: x - current position
           field - field array (potential/electric field)
    output: field(x) - field at x 
    """
    l = (x-x0)/dx #floating-point node position
    i = int(l) #grid node to the left
    di = l-i #fractional distance between nodes i and i+1
      
    return (1-di)*field[i] + di*field[i+1]
      
  #leap-frog method
v -= 0.5*(q/m)*gather(x,ef)*dt #rewind veloctiy by 0.5dt steps
  
  #csv file to output the data for further precessing/visualisation
filename = "espic1d_data.csv"
with open(filename, 'w') as csvfile:  
    csvwriter = csv.writer(csvfile)  
    #column names    
    csvwriter.writerow(['Time','X','V','KE','PE'])
    print('Time , X ,  V ,  KE , PE \t')
  
    #motion loop
    x_old = x
    for t in range(1,niter):
        ef_t = gather(x,ef) #electric field at position x
    
        x_old = x
        v += (q/m)*ef_t*dt #integrate velocity
        x += v*dt        #integrate position
        
        #potential value at x with x and V at same intant
        phi_t = gather(0.5*(x_old+x),phi) #potential at x(t-0.5dt)
    
        #energy of the particle in electron volts 
        pe = q*(phi_t-phi_max)/QE #eV, potentinal energy(PE) 
        ke = 0.5*m*v*v/QE #eV, kinetic energy(KE)

        #output to file
        csvwriter.writerow([t*dt,x,v,ke,pe])
        #screen-output every 500 iterations
        if t==1 or t%500==0:
            print(t*dt,",",x,",",v,",",ke,",",pe,"\t")


print('total time taken: ', time.time()-start)

 

 
