#hexagonal/triangular lattice
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import rand
import random
import time

#in order to construct the triangular lattice it was noted that it consists of two types of "layer" between each other
#thus if we make a 2n*n lattice of each type we can use the two of them to create a single lattice 

#defining function to randomly generate initial spin configuration
#defining two arrays: one for each "type" of row
def initialise_array(n):
    init_row1 = np.random.randint(2,size=(n,2*n))
    init_row2 = np.random.randint(2,size=(n,2*n))

    init_row1 = 2*init_row1 -1#as spin can be -1 or 1, but the previous lines give a random matrix of 0s and 1s
    init_row2 = 2*init_row2 -1
    return init_row1, init_row2   


#function to use metropolis algorithm on lattice
#each iteration performs the algorithm on one spin from each row type
def metropolis(state1,state2,beta):
    #defining random variables to be used as lattice coordinates
    a = random.randrange(0,n,2)
    b = random.randrange(0,2*n,2)
    c = random.randrange(1,n,2)
    d = random.randrange(1,2*n,2)
    
    spin_row1 =  state1[a,b]
    spin_nn1 = state2[(a+1)%n,(b+1)%n] + state2[(a-1)%n,(b+1)%n] + state2[(a+1)%n,(b-1)%n] + state2[(a-1)%n,(b-1)%n] + state1[(a+2)%n,b] + state1[(a-2)%n,b]#sum of spin of nearest neighbours
    dE1 = 2*spin_row1*spin_nn1 #change in energy if spin is flipped
    if dE1 < 0: #spin flips if energy change is negative
        spin_row1 *= -1
    elif np.exp(-dE1*beta) > rand(): #allowing positive energy change based on a probability by comparing to random number
        spin_row1 *= -1
    state1[a,b] = spin_row1 #updating lattice with new spin
    
    spin_row2 = state2[c,d]
    spin_nn2 = state1[(c+1)%n,(d+1)%n] + state1[(c-1)%n,(d+1)%n] + state1[(c+1)%n,(d-1)%n] + state1[(c-1)%n,(d-1)%n] + state2[(c+2)%n,d] + state2[(c-2)%n,d]
    dE2 = 2*spin_row2*spin_nn2 #change in energy if spin is flipped
    if dE2 < 0: #spin flips if energy change is negative
        spin_row2 *= -1
    elif np.exp(-dE2*beta) > rand(): #allowing positive energy change based on a probability by comparing to random number
        spin_row2 *= -1
    state2[c,d] = spin_row2 #updating lattice with new spin
    
    return state1,state2


#simple function to return total magnetisation of the lattice
def find_mag(state1,state2):
    #defining variables to store magnetism for each row "type"
    mag_row1 = 0
    mag_row2 = 0
    
    #summing all spins in first row "type"
    for i in range(0,n,2):
        for j in range(0,2*n,2):
            mag_row1 += state1[i,j]
            
    #summing all spins in second row "type"
    for i in range(1,n,2):
        for j in range(1,2*n,2):
            mag_row2 += state2[i,j]
        
    #summing magnetisation for both row types to find total magnetisation
    mag = mag_row1 + mag_row2
    return mag    
  

def find_energy(state1,state2):
    #defining energy variable for each row type
    energy_row1 = 0
    energy_row2 = 0
    #nested for loops to consider each point of first row type once
    for i in range(0,n,2):
        for j in range(0,2*n,2):
            spin_row1 = state1[i,j]
            #sum of spin of nearest neighbours
            neighbours1 = state2[(i+1)%n,(j+1)%n] + state2[(i-1)%n,(j+1)%n] + state2[(i-1)%n,(j+1)%n] + state2[(i-1)%n,(j-1)%n] + state1[(i+2)%n,j] + state1[(i-2)%n,j] 
            energy_row1 += -spin_row1*neighbours1
            
    #nested for loops to consider each point of second row type once       
    for i in range(1,n,2):
        for j in range(1,2*n,2):
            spin_row2 = state2[i,j]
            #sum of spin of nearest neighbours
            neighbours2 = state1[(i+1)%n,(j+1)%n] + state1[(i-1)%n,(j+1)%n] + state1[(i-1)%n,(j+1)%n] + state1[(i-1)%n,(j-1)%n] + state2[(i+2)%n,j] + state2[(i-2)%n,j]
            energy_row2 += -spin_row2*neighbours2
            
    energy = energy_row1 + energy_row2
    return energy/2
  

#defining necessary variables etc
T = np.random.normal(3.6,0.6,2**6) #creating normal distribution around T_c to show change in properties
T = T[(T>2.4) & (T<4.8)] #taking only values of T in a certain range

n = 2**3 #dimension of lattice
eqsteps = 200000 #number of steps to bring lattice to equilibrium
compsteps = 150000 #number of steps for calculation of thermodynamic properties

#defining arrays of the same size as T to store properties
magnetisation = np.zeros(len(T))
energy = np.zeros(len(T))
susceptibility = np.zeros(len(T))
specific_heat = np.zeros(len(T))

#defining constants for use in for loop below
c1 = 1.0/(compsteps*n*n)
c2 = 1.0/(compsteps*compsteps*n*n)

#main part of the code
start_time = time.time() #timing code out of interest

#loop to repeatedly apply the algorithm to the lattice
for j in range(len(T)):
    #initialising variables for various properties
    mag = 0
    ene = 0
    sus = 0
    spec = 0
    
    #initialising arrays for each row "type"
    state1,state2 = initialise_array(n)
    beta = 1.0/T[j] 
    
    #system is brought to equilbrium 
    for i in range(eqsteps): 
        metropolis(state1,state2,beta)
     
    #properties are examined as system is now in equilibirum
    for i in range(compsteps):
        metropolis(state1,state2,beta)
        m = find_mag(state1,state2)
        e = find_energy(state1,state2)
        
        #adding to the variables to find property
        mag += m
        ene += e
        sus += m*m
        spec += e*e
        
    #using formulas to determine each property for the given temperature    
    magnetisation[j] = c1*mag
    energy[j] = c1*ene
    susceptibility[j] = (c1*sus - c2*mag*mag)/T[j]
    specific_heat[j] = (c1*spec - c2*ene*ene)/(T[j]**2)

print ("for loop took", (time.time() - start_time)/60, "minutes")

#creating one figure with all 4 plots using subplot from pyplot package 
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(18, 10)); # plot the calculated values    

ax1.plot(T, energy, 'bo');
ax1.set_xlabel("Temperature (T)", fontsize=15);
ax1.set_ylabel("Energy ", fontsize=15);
ax1.set_xticks(np.arange(2.4,4.8,0.25))
ax1.axvline(x=3.6,color='k',alpha=0.7,linestyle='dashed')

ax2.plot(T, np.fabs(magnetisation), 'ro');
ax2.set_xlabel("Temperature (T)", fontsize=15);
ax2.set_ylabel("Magnetization", fontsize=15);
ax2.set_xticks(np.arange(2.4,4.8,0.25))
ax2.axvline(x=3.6,color='k',alpha=0.7,linestyle='dashed')

ax3.plot(T, specific_heat, 'bo');
ax3.set_xlabel("Temperature (T)", fontsize=15);
ax3.set_ylabel("Specific Heat ", fontsize=15);
ax3.set_xticks(np.arange(2.4,4.8,0.25))
ax3.axvline(x=3.6,color='k',alpha=0.7,linestyle='dashed')

ax4.plot(T, susceptibility, 'ro');
ax4.set_xlabel("Temperature (T)", fontsize=15);
ax4.set_ylabel("Susceptibility", fontsize=15);
ax4.set_xticks(np.arange(2.4,4.8,0.25))
ax4.axvline(x=3.6,color='k',alpha=0.7,linestyle='dashed')

plt.show() 
    
    
    
    