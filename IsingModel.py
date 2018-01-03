import matplotlib.pyplot as plt
import numpy as np
from numpy.random import rand
import scipy.constants


#defining function to randomly generate initial spin configuration
def initialise(n):
    state = np.random.randint(2,size=(n,n))
    return 2*state -1   #as spin can be -1 or 1, but the previous line gives a random matrix of 0 and 1

#function to use metropolis algorithm on lattice
def metropolis(state,beta,n):
    a = np.random.randint(0,n)
    b = np.random.randint(0,n)
    spin_ij = state[a,b] #considering the spin of a random point in the lattice
    spin_neighbours = state[(a+1)%n,b] + state[(a-1)%n,b] + state[a,(b+1)%n] + state[a,(b+1)%n] #considering spin of the neighbours #by using modulo n, periodic boundary conditions can be established
    hamil = -spin_ij*spin_neighbours #- spin_ij
    dE = -2*hamil #as dE = flipped Hamiltonian - Hamiltonian = -2*Hamiltonian
    if dE < 0:
        spin_ij *= -1.0
    #here we compare the probability to a random number between 0 and 1 to determine whether or not the spin will flip
    elif np.exp(-dE*beta) > rand():
        spin_ij *= -1.0
    state[a,b] = spin_ij
    return state

#function that plots the current spin state of the lattice when called
def plot(state,i):
    plt.figure(figsize=(5,5))
    plt.imshow(state)
    plt.title("n = %s" %i)

#simple function to return total magnetisation of the lattice
def find_mag(state):
    mag = np.sum(state)
    return mag


#defining necessary variables etc
k = scipy.constants.k
T_arr = np.linspace(1.1,4.0,58)
T = T_arr/k

n = 16
eqsteps = 1000000  #number of steps to bring lattice to equilibrium
calcsteps = 500000 #number of steps for calculation of thermodynamic properties
magnetisation = []

#main part of the code

#loop to repeatedly apply the algorithm to the lattice
for j in range(len(T)):
    mag = 0
    state = initialise(n)
    beta = 1/(k*T[j])
    
    for i in range(eqsteps): #system is brought to equilbrium 
        metropolis(state,beta,n)
            
    for i in range(calcsteps):
        metropolis(state,beta,n)
        mag += find_mag(state)
    magnetisation.append(mag/(calcsteps*n*n))
    
print(magnetisation)
fig = plt.plot(T,np.fabs(magnetisation),'ro')
plt.xlabel('Temperature (T)')
plt.ylabel('Magnetisation')
plt.show()
    
    
    
    
    
    