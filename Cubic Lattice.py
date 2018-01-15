import matplotlib.pyplot as plt
import numpy as np
from numpy.random import rand
import time


#defining function to randomly generate initial spin configuration
def initialise_array(n):
    init = np.random.randint(2,size=(n,n,n)) #initialising 3d array (rank-3 tensor)
    return 2*init -1   #as spin can be -1 or 1, but the previous line gives a random tensor of 0's and 1's


#function to use metropolis algorithm on lattice
def metropolis(state, beta):
    a = np.random.randint(0,n)
    b = np.random.randint(0,n)
    c = np.random.randint(0,n)
    spin_ij =  state[a,b,c] #randomly choosing point in the lattice
    spin_nn = state[(a+1)%n,b,c] + state[a,(b+1)%n,c] + state[a,b,(c+1)%n] + state[(a-1)%n,b,c] + state[a,(b-1)%n,c] + state[a,b,(c-1)%n] #sum of spin of nearest neighbours
    dE = 2*spin_ij*spin_nn #change in energy if spin is flipped
    if dE < 0:
        spin_ij *= -1
    elif np.exp(-dE*beta) > rand(): #allowing positive energy change based on a probability by comparing to random number
        spin_ij *= -1
    state[a,b] = spin_ij
    return state


#simple function to return total magnetisation of the lattice
def find_mag(state):
    mag = np.sum(state) #summing all spins in the lattice
    return mag    
  
#function to find total energy of the lattice
def find_energy(state):
    energy = 0 #defining energy variable
    #nested for loops to consider each point of lattice once
    for i in range(0,n):
        for j in range(0,n):
            for k in range(0,n):
                spin_ij = state[i,j,k]
                #total spin of nearest neighbours
                neighbours = state[(i+1)%n,j,k] + state[(i-1)%n,j,k] + state[i,(j+1)%n,k] + state[i,(j-1)%n,k] + state[i,j,(k+1)%n] + state[i,j,(k-1)%n] 
                energy += -spin_ij*neighbours #adding energy of each point to the energy
    return energy/2
  

#defining necessary variables etc
T = np.random.normal(4.5,0.65,2**6) #creating normal distribution around T_c to show change in properties
T = T[(T>3.5) & (T<6.0)] #taking only values of T in a certain range

n = 2**3 #dimension of lattice
eqsteps = 200000 #number of steps to bring lattice to equilibrium
compsteps = 100000 #number of steps for calculation of thermodynamic properties

#defining arrays to store values of physical properties
magnetisation = np.zeros(len(T))
energy = np.zeros(len(T))
susceptibility = np.zeros(len(T))
specific_heat = np.zeros(len(T))

#defining constants for use in for loop below
c1 = 1.0/(compsteps*n*n*n)
c2 = 1.0/(compsteps*compsteps*n*n*n)

#main part of the code
start_time = time.time() #timing code out of interest

#loop to repeatedly apply the algorithm to the lattice
for j in range(len(T)):
    #initialising variables for various properties
    mag = 0
    ene = 0
    sus = 0
    spec = 0
    state = initialise_array(n) #iniitialising array
    beta = 1.0/T[j]
    
    #system is brought to equilbrium
    for i in range(eqsteps):  
        metropolis(state,beta)
     
    #now system is in equilibirium so properties can be investigated
    for i in range(compsteps): 
        metropolis(state,beta)
        #calculating magnetism and energy for each configuration
        m = find_mag(state)
        e = find_energy(state)
        
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

print ("for loop took", (time.time() - start_time)/60, "minutes") #printing time taken for loop

#creating one figure with all 4 plots using subplot from pyplot package 
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(18, 10)); # plot the calculated values    

#defining subplots to display behaviour of each of the properties wrt T
ax1.plot(T, energy, 'bo');
ax1.set_xlabel("Temperature (T)", fontsize=15);
ax1.set_ylabel("Energy ", fontsize=15);
ax1.set_xticks(np.arange(3.5,6.0,0.25))
ax1.axvline(x=4.5,color='k',alpha=0.7,linestyle='dashed')

ax2.plot(T, np.fabs(magnetisation), 'ro');
ax2.set_xlabel("Temperature (T)", fontsize=15);
ax2.set_ylabel("Magnetization", fontsize=15);
ax2.set_xticks(np.arange(3.5,6.0,0.25))
ax2.axvline(x=4.5,color='k',alpha=0.7,linestyle='dashed')

ax3.plot(T, specific_heat, 'bo');
ax3.set_xlabel("Temperature (T)", fontsize=15);
ax3.set_ylabel("Specific Heat ", fontsize=15);
ax3.set_xticks(np.arange(3.5,6.0,0.25))
ax3.axvline(x=4.5,color='k',alpha=0.7,linestyle='dashed')

ax4.plot(T, susceptibility, 'ro');
ax4.set_xlabel("Temperature (T)", fontsize=15);
ax4.set_ylabel("Susceptibility", fontsize=15);
ax4.set_xticks(np.arange(3.5,6.0,0.25))
ax4.axvline(x=4.5,color='k',alpha=0.7,linestyle='dashed')

plt.show() 
    
    
    
    