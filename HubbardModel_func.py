# -*- coding: utf-8 -*-

from itertools import permutations
import numpy as np

# returns the basis
#
# __parameters__
# model : list of values, [ number of up particles,
#                           number of down particles,
#                           number of sites,
#                           category                  ]
#
# __returns__
# basis : list of all possible states
#
def get_basis(model):
    
    nUp, nDown, nSites, cat = model
    
    # list of nUp 1s and nSites-nUp 0s (nDown 1js)
    up = [complex(1,0) for n in range(nUp)]
    up += [complex() for n in range(nSites-nUp)]
    down  = [complex(0,1) for n in range(nDown)]
    down += [complex() for n in range(nSites-nDown)]

    # get permutations
    upPerms = [u for u in set(permutations(up,nSites))]
    downPerms = [d for d in set(permutations(down,nSites))]

    # combine permutations
    basis = [np.add(u,d).tolist() for u in upPerms for d in downPerms]
        
    return basis

#
# helper function for get_hamiltonian
# 
# given an initial state,
# returns a list of the possible states after one perturbation
#
# __parameters__
# initState : intial state
#
# __returns__
# possibleStates : list of possible states after one perturbation
#
def possible_states(initState):
    
    possibleStates = []
    
    for s,site in enumerate(initState):
        
        if site.real == 1:
            
            # move to the left
            if s != 0 and initState[s-1].real == 0:
                newState = initState.copy()
                newState[s] -= 1
                newState[s-1] += 1
                possibleStates.append(newState)

            # move to the right
            if s != len(initState)-1 and initState[s+1].real == 0:
                newState = initState.copy()
                newState[s] -= 1
                newState[s+1] += 1
                possibleStates.append(newState)
                
        if site.imag == 1:
                
            # move to the left
            if s != 0 and initState[s-1].imag == 0:
                newState = initState.copy()
                newState[s] -= 1j
                newState[s-1] += 1j
                possibleStates.append(newState)

            # move to the right
            if s != len(initState)-1 and initState[s+1].imag == 0:
                newState = initState.copy()
                newState[s] -= 1j
                newState[s+1] += 1j
                possibleStates.append(newState)
                
    return possibleStates

#
# helper function for get_hamiltonian_cyclic
#
# treats the system as cyclic
# 
# given an initial state,
# returns a list of the possible states after one perturbation
#
# __parameters__
# initState : intial state
#
# __returns__
# possibleStates : list of possible states after one perturbation
#
def possible_states_cyclic(initState):
    possibleStates = []
    for s,site in enumerate(initState):
        
        if site.real == 1:
            
            # move to the left
            if s != 0 and initState[s-1].real == 0:
                newState = initState.copy()
                newState[s] -= 1
                newState[s-1] += 1
                possibleStates.append(newState)
                
            if s == 0 and initState[len(initState)-1].real == 0:
                newState = initState.copy()
                newState[s] -= 1
                newState[len(initState)-1] += 1
                possibleStates.append(newState)

            # move to the right
            if s != len(initState)-1 and initState[s+1].real == 0:
                newState = initState.copy()
                newState[s] -= 1
                newState[s+1] += 1
                possibleStates.append(newState)
                
            if s == len(initState)-1 and initState[0].real == 0:
                newState = initState.copy()
                newState[s] -= 1
                newState[0] += 1
                possibleStates.append(newState)
                
        if site.imag == 1:
                
            # move to the left
            if s != 0 and initState[s-1].imag == 0:
                newState = initState.copy()
                newState[s] -= 1j
                newState[s-1] += 1j
                possibleStates.append(newState)
                
            if s == 0 and initState[len(initState)-1].imag == 0:
                newState = initState.copy()
                newState[s] -= 1j
                newState[len(initState)-1] += 1j
                possibleStates.append(newState)

            # move to the right
            if s != len(initState)-1 and initState[s+1].imag == 0:
                newState = initState.copy()
                newState[s] -= 1j
                newState[s+1] += 1j
                possibleStates.append(newState)
                
            if s == len(initState)-1 and initState[0].imag == 0:
                newState = initState.copy()
                newState[s] -= 1
                newState[0] += 1
                possibleStates.append(newState)
                
    return possibleStates

# 
# returns the Hsubbard hamiltonian
#
# __parameters__
# U       : value of U
# t       : value of t
# basis   : list of all possible states
# nStates : number of states in the basis
#
# __returns__
# hamiltonian : Hubbard hamiltonian matrix
#
def get_hamiltonian(U,t,basis,nStates):

    # allocate hamiltonian
    hamiltonian = np.zeros((nStates,nStates),dtype=float)
    
    mu = U / 2
    particles = sum([ site.real + site.imag for site in basis[0] ])
    
    # iterate through basis
    for m in range(nStates):
        
        # U * double occupancies
        hamiltonian[m][m] = U * basis[m].count(1+1j) - mu * particles
        
        # add -t
        possibleStates = possible_states(basis[m])
        for n in range(m+1,nStates):
            if list(basis[n]) in possibleStates:
                hamiltonian[m][n] = -t
                hamiltonian[n][m] = -t

    return hamiltonian

# 
# treats the system as cyclic
# 
# returns the Hsubbard hamiltonian
#
# __parameters__
# U       : value of U
# t       : value of t
# basis   : 
# nStates : number of states in the basis
#
# __returns__
# hamiltonian : Hubbard hamiltonian matrix
#
def get_hamiltonian_cyclic(U,t,basis,nStates):
    
    mu = U/2
    
    # number of particles
    particles = sum([site.real + site.imag for site in basis[0]])
    
    # allocate hamiltonian matrix
    hamiltonian = np.zeros((nStates,nStates),dtype=float)
    
    # iterate through basis
    for m in range(nStates):
        
        # U * double occupancies - mu * number of particles
        hamiltonian[m][m] = U * basis[m].count(1+1j) - mu * particles
        
        # add -t
        possibleStates = possible_states_cyclic(basis[m])
        for n in range(m+1,nStates):
            if list(basis[n]) in possibleStates:
                hamiltonian[m][n] = -t
                hamiltonian[n][m] = -t

    return hamiltonian

# 
# returns a list of basis vectors in category
# (starts with repeated antiferromagnetic state)
#
# __parameters__
# model : list of values, [ number of up particles,
#                           number of down particles,
#                           number of sites,
#                           category                  ]
# basis   : list of all possible states
# nStates : number of states in the basis
#
# __returns__
# psi : list of basis vectors
#       (starting with extra antiferromagnetic state)
#
def get_psi(model,basis,nStates):
    
    nUp, nDown, nSites, category = model
    
    # convert category to int
    cat = {"anti": 0, "single": 1, "double": 2}
    category = cat[category]
    
    # return state indexes
    state_indexes = []
    
    # antiferromagnetic states
    half = [complex()] * nSites
    mirror = [complex()] * nSites
    
    # fill half state
    if nUp > nDown:
        half = [complex(1,0) if n%2 == 0 else complex(0,1) for n in range(nSites)]
    else:
        half = [complex(0,1) if n%2 == 0 else complex(1,0) for n in range(nSites)]
    
    # add half filled at state_indexes[0]      
    state_indexes.append(basis.index(half))
    
    if category == 0:
        
        state_indexes.append(basis.index(half))
        if nSites % 2 == 0:
            mirror = [complex(h.imag,h.real) for h in half]
            state_indexes.append(basis.index(mirror))
            
    elif category == 1 or category == 2:
        
        # add states of interest
        for i in range(nStates):
            if 1+1j not in basis[i]:
                state_indexes.append(i)
        
    else:
        print('category must equal 0, 1, or 2')

    # basis vectors
    psi = []
    for i in state_indexes:
        psi.append([[1] if j == i else [0] for j in range(nStates)])
                
    return psi
