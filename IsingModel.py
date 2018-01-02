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
    a = np.random.randint(0,n)
    b = np.random.randint(0,n)
    spin_ij = state[a,b] #considering the spin of a random point in the lattice
    spin_neighbours = state[(a+1)%n,b] + state[(a-1)%n,b] + state[a,(b+1)%n] + state[a,(b+1)%n] #considering spin of the neighbours #by using modulo n, periodic boundary conditions can be established
    dE = 2*spin_ij*spin_neighbours
    if dE < 0:
        spin_ij *= -1.0
    #here we compare the probability to a random number between 0 and 1 to determine whether or not the spin will flip
    elif np.exp(-dE*beta) > rand():
        spin_ij *= -1.0
    state[a,b] = spin_ij
    return state

def plot(result):
    plt.figure(figsize=(5,5))
    plt.imshow(result)
    
#main function                      
n = 100
nsteps = 1000000


state = initialise(n)
print (state)

#loop to repeatedly apply the algorithm to the lattice
for i in range(nsteps):  
    result = metropolis(state,beta,n)
    if i == 1:
        plot(result)
    if i == 100:
        plot(result)
    if i == 1000:
        plot(result)
    if i == 10000:
        plot(result)
    if i == 100000:
        plot(result)
    if i == 1000000:
        plot(result)

    
print (result)
plot(result)