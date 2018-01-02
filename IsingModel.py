import matplotlib.pyplot as plt
import numpy as np
from numpy.random import rand
import scipy.constants

T = 273.14
k = scipy.constants.k
beta = 1.0/k*T

#defining function to randomly generate initial spin configuration
def initialise(n):
    state = np.random.randint(2,size=(n,n))
    return 2*state -1   #as spin can be -1 or 1, but the previous line gives a random matrix of 0 and 1

#function to use metropolis algorithm on lattice
def metropolis(state,beta,n):
    for i in range(n):
        for j in range(n):
            a = np.random.randint(0,n)
            b = np.random.randint(0,n)
            spin_ij = state[a,b] #considering the spin of a random point in the latticr
            spin_neighbours = state[(a+1)%n,b] + state[(a-1)%n,b] + state[a,(b+1)%n] + state[a,(b+1)%n] #considering spin of the neighbours #by using modulo n, periodic boundary conditions can be established
            dE = 2*spin_ij*spin_neighbours
            if dE < 0:
                spin_ij *= -1.0
            #here we compare the Boltzmann distribution to a random integer to determine whether or not the spin will flip
            elif np.exp(-dE*beta) > rand():
                spin_ij *= -1.0
            state[a,b] = spin_ij
    return state


#main function                      
n = 2**6
nsteps = 2*11

state = initialise(n)
print (state)

#loop to repeatedly apply the algorithm to the lattice
for i in range(nsteps):  
    result = metropolis(state,beta,n)
    
print (result)
f = plt.figure(figsize=(5,5))
plt.imshow(result)

