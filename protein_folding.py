# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 10:34:09 2018

@author: John
"""

import numpy as np
from matplotlib import pyplot as plt
from numba import jit
import math

np.set_printoptions(threshold=np.nan)





T = 10
N = 15
MCS = 10000000


# Interaction energies between different types of acids
# Interaction energy between type i and type j is J[i, j]
J = np.random.uniform(-4, -2, size=(20, 20))
J = np.hstack((np.zeros((J.shape[0],1)),J)) # Type 0 doesn't exist, has an interaction energy
J = np.vstack((np.zeros((1,J.shape[1])),J)) # of 0 with every type


# Initially, no acids are non-bonded nearest-neighbours
E_chain = 0


# Make 2N by 2N lattice to ensure that the protein has rooom to fold
lattice = np.zeros(shape=(2*N,2*N), dtype=int)


# The type of each acid in the chain
types = np.random.randint(1, 21, size=N, dtype=int)


# Create x and y coordinates such that the chain will be horizontal, and centered in the 2N x 2N lattice
x_start = int(N/2)
x_end = int(3*N/2)
pos_xs = np.arange(x_start, x_end, dtype=int)
pos_ys = np.empty(N, dtype=int); pos_ys.fill(N)


# Zip type, x_pos, and y_pos to create array defining each acid (this is the protein)
A = np.empty(shape=(N,3), dtype=int)
for i in range(A.shape[0]):
    A[i] = np.array([types[i], pos_xs[i], pos_ys[i]])
    
    
    
# Put protein on lattice
lattice[N, x_start:x_end] = types


# Pre-select random acid and neighbour picks for better performance
acid_picks = np.random.randint(0, A.shape[0], size=MCS) # 1 acid picked from A each MCS
curr_acid_pick = 0

neigh_picks = np.random.randint(0, 4, size=MCS) # 1 acid picked from A each MCS
curr_neigh_pick = 0

prob_picks = np.random.uniform(0, 1, size=MCS)
curr_prob_pick = 0

#left-down, left-up, right-up, right-down
move_neighs = np.array([(-1,1), (-1,-1), (1,-1), (1,1)], dtype=int)




@jit
def get_neighbours(x, y):
    left = (x-1, y) if not x == 0 else None
    bottom = (x, y+1) if not y == lattice.shape[1]-1 else None
    right = (x+1, y) if not x == lattice.shape[0]-1 else None
    top = (x, y-1) if not y == 0 else None
        
    return [left, bottom, right, top]


@jit
def energy(a_type, x, y, partners):
    
    neighbours = get_neighbours(x, y)
    
    #print("\nEnergy:\nFor " + str(x) + ", " + str(y))
        
    E = 0
    for neigh in neighbours:
        
        if neigh not in partners:
            
            #print(str(neigh) + " " + str(lattice[neigh[1], neigh[0]]))
            
            neigh_type = lattice[neigh[1], neigh[0]] # SLOW
            E += J[a_type, neigh_type]
      
    #print("E: " + str(E))
    return E




# Selects one acid, determines the change in energy from moving that acid, then either moves it or doesn't
#@jit
def sweep(acid_index, neigh_index, prob_index): 
        
    #global A 
    #global E_chain
    
    # Randomly select one acid in the chain
    acid_pos = acid_picks[acid_index]
    
    print(acid_pos)
    
    acid = A[acid_pos]
    acid_x = acid[1]
    acid_y = acid[2]
    
    print(acid)
    
    #print("Selected acid: " + str(acid))
    
    
    # Determine the acids the currently examined acid is bonded to 
    first_partner = tuple(A[acid_pos-1, 1:]) if not acid_pos == 0 else None
    second_partner = tuple(A[acid_pos+1, 1:]) if not acid_pos == len(A)-1 else None
    partners = (first_partner, second_partner)
    
    #print("Partners: " + str(partners))
        
    
    # Randomly select a nearest neighbour to move to
    neigh_type = neigh_picks[neigh_index]
    
    neigh_offset = move_neighs[neigh_type]
    
    neigh_pos = (acid_x+neigh_offset[0], acid_y+neigh_offset[1]) # creating a tuple is faster than assigning a list
    
    #print("Proposed new position: " + str(neigh_pos))
    
    
    # Determine if acid can move to selected position without breaking bonds  
    can_move = True
    
    if lattice[neigh_pos[1], neigh_pos[0]] != 0:
        can_move = False
        
    else:
        for partner in partners:
            if partner is None:
                continue
            elif not (abs(partner[0] - neigh_pos[0]) + abs(partner[1] - neigh_pos[1]) == 1):
                can_move = False
    
    #print(can_move)
    
    energy(acid[0], acid[1], acid[2], partners)
    
    if can_move:
        1
        
        E_i = energy(acid[0], acid_x, acid_y, partners)
        E_f = energy(acid[0], neigh_pos[0], neigh_pos[1], partners)
        delta_E = E_f-E_i
        #print("\nDELTA_E: " + str(delta_E))
        
        if delta_E < 0 or prob_picks[prob_index] < math.exp(delta_E/T):
            #print("Moving from (" + str(acid_x) + ", " + str(acid_y) + ") to (" + str(neigh_pos[0]) + ", " + 
                  #str(neigh_pos[1]) + ")")
            lattice[acid_y, acid_x] = 0
            lattice[neigh_pos[1], neigh_pos[0]] = acid[0]
            A[acid_pos] = [acid[0], neigh_pos[0], neigh_pos[1]]
            
            #blank = E_chain
            
            #E_chain += delta_E
            
            #print("changed A at " + str(acid_pos))
            #print(A)

            # Show the protein on the lattice
            #plt.imshow(lattice)
            #plt.show()
        
    #return A
    
    
    
    

for i in range(10):
    print(A)
    sweep(curr_acid_pick, curr_neigh_pick, curr_prob_pick)
    curr_acid_pick += 1
    curr_neigh_pick += 1
    curr_prob_pick += 1
    
# Show the protein on the lattice
#plt.imshow(lattice)
#plt.show()