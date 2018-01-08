import matplotlib.pyplot as plt
import numpy as np
from numpy.random import rand
import time


#defining function to randomly generate initial spin configuration
def initialise_array(n):
    init = np.random.randint(2,size=n) #initialising 1-dimensional array of length n with 0's and 1's
    return 2*init -1   #as spin can be -1 or 1, but the previous line gives a random matrix of 0 and 1


#function to use metropolis algorithm on lattice
def metropolis(state, beta):
    a = np.random.randint(0,n) #choosing random coordinate in lattice
    spin = state[a] #spin of a given point
    neighbours = state[(a+1)%n] + state[(a-1)%n] #sum of spin of nearest neighbours
    dE = 2*spin*neighbours #change in energy if the spin is flipped
    if dE < 0:
        spin *= -1
    elif np.exp(-dE*beta) > rand(): #allowing spin to flip based on a probability by comparison to a random number
        spin *= -1
    state[a] = spin
    return state


#simple function to return total magnetisation of the lattice
def find_mag(state):
    mag = np.sum(state) #summing all the spins in the lattice in a given configuration
    return mag    
  
#defining function to calculate energy of a given configuration
def find_energy(state):
    energy = 0
    for i in range(0,n):
        spin = state[i] #spin of a given point
        neighbours = state[(i+1)%n] + state[(i-1)%n] #spin of the points two neighbours
        energy += -spin*neighbours
    return energy/4
  

#defining necessary variables etc
T = np.random.normal(0.5,0.5,2**7) #creating normal distribution around T_c to show change in properties
T = T[(T>0) & (T<1.0)] #taking only values of T in a certain range
T = np.append(T,0.000001) #appending value to T array to show critical temperature at T=0
print(T)

n = 2**7 #dimension of lattice
eqsteps = 1000000 #number of steps to bring lattice to equilibrium
compsteps = 500000 #number of steps for calculation of thermodynamic properties

#creating arrays to store values of physical properties
magnetisation = np.zeros(len(T))
energy = np.zeros(len(T))
susceptibility = np.zeros(len(T))
specific_heat = np.zeros(len(T))

#defining constants for use in for loop below
c1 = 1.0/(compsteps*n)
c2 = 1.0/(compsteps*compsteps*n)

#main part of the code
start_time = time.time() #timing code out of interest

#loop to repeatedly apply the algorithm to the lattice
for j in range(len(T)):
    #creating variables to store values of physical properties
    mag = 0
    ene = 0
    sus = 0
    spec = 0
    state = initialise_array(n) #initialising array
    beta = 1.0/T[j] 
    
    for i in range(eqsteps): #system is brought to equilbrium 
        metropolis(state,beta)
            
    for i in range(compsteps): #properties can be examined now as system is in equilibrium
        metropolis(state,beta)
        #calculating magnetism and energy for each configuration
        m = find_mag(state)
        e = find_energy(state)
        #adding to the variables to find property
        mag += m
        ene += e
        sus += m*m
        spec += e*e
    
    #appending final value to the arrays for each property
    magnetisation[j] = c1*mag
    energy[j] = c1*ene
    #implementing formulas to calculate susceptibility & specific heat from values of magnetism & energy respectively 
    susceptibility[j] = (c1*sus - c2*mag*mag)/T[j] 
    specific_heat[j] = (c1*spec - c2*ene*ene)/(T[j]**2)

print ("for loop took", (time.time() - start_time)/60, "minutes") #printing time taken for loop

#creating one figure with all 4 plots using subplot from pyplot package 
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(18, 10)); # plot the calculated values    

#defining subplots to display behaviour of each of the properties wrt T
ax1.plot(T, energy, 'bo');
ax1.set_xlabel("Temperature (T)", fontsize=15);
ax1.set_ylabel("Energy ", fontsize=15);
ax1.set_xticks(np.arange(0,1.0,0.2))

ax2.plot(T, np.fabs(magnetisation), 'ro');
ax2.set_xlabel("Temperature (T)", fontsize=15);
ax2.set_ylabel("Magnetization", fontsize=15);
ax2.set_xticks(np.arange(0,1.0,0.2))

ax3.plot(T, specific_heat, 'bo');
ax3.set_xlabel("Temperature (T)", fontsize=15);
ax3.set_ylabel("Specific Heat ", fontsize=15);
ax3.set_xticks(np.arange(0,1.0,0.2))

ax4.plot(T, susceptibility, 'ro');
ax4.set_xlabel("Temperature (T)", fontsize=15);
ax4.set_ylabel("Susceptibility", fontsize=15);
ax4.set_xticks(np.arange(0,1.0,0.2))

plt.show() 
    
    
    
    