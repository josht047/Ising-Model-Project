import matplotlib.pyplot as plt
import numpy as np
from numpy.random import rand
import time


#defining function to randomly generate initial spin configuration
def initialise_array(n):
    init = np.random.randint(2,size=(n,n))
    return 2*init -1   #as spin can be -1 or 1, but the previous line gives a random matrix of 0 and 1

#function to use metropolis algorithm on lattice
def metropolis(state, beta):
    a = np.random.randint(0,n)
    b = np.random.randint(0,n)
    spin_ij =  state[a,b]
    spin_nn = state[(a+1)%n,b] + state[a,(b+1)%n] + state[(a-1)%n,b] + state[a,(b-1)%n] #sum of spin of nearest neighbours
    dE = 2*spin_ij*spin_nn
    if dE < 0:
        spin_ij *= -1
    elif np.exp(-dE*beta) > rand():
        spin_ij *= -1
    state[a,b] = spin_ij
    return state

#function that plots the current spin state of the lattice when called

#simple function to return total magnetisation of the lattice
def find_mag(state):
    mag = np.sum(state)
    return mag    
       
    
#defining necessary variables etc
T = np.random.normal(2.26,0.5,2**5) #creating normal distribution around T_c to easier see change in magnetisation
T = T[(T>1.1) & (T<3.8)] #taking only values of T in a certain range

n = 2**4
eqsteps = 500000 #number of steps to bring lattice to equilibrium
compsteps = 100000 #number of steps for calculation of thermodynamic properties
magnetisation = []

#main part of the code
start_time = time.time() #timing code out of interest

#loop to repeatedly apply the algorithm to the lattice
for j in range(len(T)):
    mag = 0
    state = initialise_array(n)
    b = 1.0/T[j]
    
    for i in range(eqsteps): #system is brought to equilbrium 
        metropolis(state,b)
            
    for i in range(compsteps):
        metropolis(state,b)
        m1 = find_mag(state)
        mag += m1
    magnetisation.append(mag/(compsteps*n*n))
    
print ("for loop took", (time.time() - start_time)/60, "minutes")

plt.plot(T,np.fabs(magnetisation),'ro')
plt.title('Magnetization against Temperature')
plt.xlabel('Temperature (T) in J/k')
plt.ylabel('Magnetisation')
plt.ylim(-0.05,1.05)
plt.show()
    
    
    
    
    
    