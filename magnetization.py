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
  

def find_energy(state):
    energy = 0
    for i in range(0,n):
        for j in range(0,n):
            spin_ij = state[i,j]
            neighbours = state[(i+1)%n,j] + state[(i-1)%n,j] + state[i,(j+1)%n] + state[i, (j-1)%n]
            energy += -spin_ij*neighbours
    return energy
    
#defining necessary variables etc
T = np.random.normal(2.26,0.5,2**5) #creating normal distribution around T_c to easier see change in magnetisation
T = T[(T>1.1) & (T<3.8)] #taking only values of T in a certain range

n = 2**4
eqsteps = 500000 #number of steps to bring lattice to equilibrium
compsteps = 100000 #number of steps for calculation of thermodynamic properties
magnetisation = energy = susceptibility = specific_heat = []


#main part of the code
start_time = time.time() #timing code out of interest

#loop to repeatedly apply the algorithm to the lattice
for j in range(len(T)):
    mag = ene = sus = spe = 0
    state = initialise_array(n)
    beta = 1.0/T[j]
    
    for i in range(eqsteps): #system is brought to equilbrium 
        metropolis(state,beta)
            
    for i in range(compsteps):
        metropolis(state,beta)
        m = find_mag(state)
        e = find_energy(state)
        mag += m
        ene += e
        sus += m**2
        spe += e**2
        
    magnetisation.append(mag/(compsteps*n*n))
    energy.append(ene/(compsteps*n*n))
    susceptibility.append(((mag/(compsteps*n*n))-(((mag)**2)/(compsteps*compsteps*n*n)))*beta)
    specific_heat.append(((ene/(compsteps*n*n))-(((ene)**2)/(compsteps*compsteps*n*n)))*(beta**2))

print ("for loop took", (time.time() - start_time)/60, "minutes")

fig = plt.figure()

fig.add_subplot(122)
plt.plot(T,np.fabs(magnetisation),'ro')
plt.title('Magnetization against Temperature')
plt.xlabel('Temperature (T) in J/k')
plt.ylabel('Magnetisation')
plt.ylim(-0.05,1.05)

fig.add_subplot(222)
plt.plot(T,np.fabs(energy),'ko')
plt.title('Energy against Temperature')
plt.xlabel('Temperature (T) in J/k')
plt.ylabel('Energy')

fig.add_subplot(322)
plt.plot(T,np.fabs(susceptibility),'ro')
plt.title('Susceptibility against Temperature')
plt.xlabel('Temperature (T) in J/k')
plt.ylabel('Susceptibility')

fig.add_subplot(422)
plt.plot(T,np.fabs(specific_heat),'ko')
plt.title('Specific Heat against Temperature')
plt.xlabel('Temperature (T) in J/k')
plt.ylabel('Specific Heat')

plt.show()
    
    
    
    
    
    