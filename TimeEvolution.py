# -*- coding: utf-8 -*-


from anytree import Node
import HubbardModel_func as hm
import matplotlib.pyplot as plt
from multiprocessing import cpu_count, Manager, Pool
import numpy as np
from scipy.linalg import expm
from time import perf_counter as timer

#
# returns the basis of a half filled system and 
#   basis vectors of the category
#
# __parameters__
# nSites   : number of states in the system
# category : states of interest,
#            { "anti": antiferromagnetic states, 
#              "single": singly occupied states,
#              "double": doubly occupied states  }
#
# __return__
# basis : number of sites in the system
# psi   : list of basis vectors in category (always starts with one antiferromagnetic state)
# 
#
def setup_model(nSites,category):
    
    # number of particles for half filled (nUp > nDown if odd)
    nDown = int(nSites / 2)
    if nSites % 2 == 0:
        nUp = int(nSites / 2)
    else:
        nUp = int(nSites / 2) + 1

    model = [nUp,nDown,nSites,category]

    basis = hm.get_basis(model)
    
    nStates = len(basis)

    # states of interest
    psi = hm.get_psi(model,basis,nStates)
    
    return basis, psi


#
# off diagonals of Heisenberg
# given one state and the basis,
# returns the indexes where particles can swap
#
# __parameters__
# state : list of complex numbers that represents a state
# basis : list of all possible states
#
# __returns__
# idxs : indexes of off diagonals in Heisenberg
#
def exchanged_states(state: list, basis: list):
    
    idxs = []
    
    for i,site in enumerate(state[:-1]):
        if site == 1 and state[i+1] == 1j:
            add = state[:]
            add[i] = 1j
            add[i+1] = 1
            idxs.append(basis.index(add))
            
        if site == 1j and state[i+1] == 1:
            add = state[:]
            add[i] = 1
            add[i+1] = 1j
            idxs.append(basis.index(add))
            
    return idxs


#
# Returns the Heisenberg matrix
#
# __parameters__
# basis : list of all possible states
#
# __returns__
# heisenberg : Heisenberg matrix
#
def get_heisenberg(basis):
    
    nStates = len(basis)
    
    # allocate hamiltonian
    heisenberg = np.zeros((nStates,nStates),dtype=float)
    
    # iterate through basis
    for m in range(nStates):
        
        state = basis[m]
        
        # diagonal                                  
        n_sum = 0
        # sum([(site.real - site.imag) * (state[i+1].real - state[i+1].imag) if site != 1+1 or state[i+1] != 1+1j for i,site in enumerate(state[:-1])])
        for i,site in enumerate(state[:-1]):
            if site == 1+1 or state[i+1] == 1+1j:
                n_sum += 0
            else:
                n_sum += (site.real - site.imag) * (state[i+1].real - state[i+1].imag)
        
        heisenberg[m][m] = n_sum
        
        # off diagonals
        exchangedStates = exchanged_states(state, basis)
        for i in exchangedStates:
            if i > m:
                heisenberg[m][i] = 2
                heisenberg[i][m] = 2

    return heisenberg

#
# returns a list of time operators for each U
#
# __parameters__
# U     : list of U values
# t     : t value
# dtau  : delta tau
# basis : list of all possible states
#
# __returns__
# timeOperators : list of time operators
#
def get_timeOperators(U,t,dtau,basis):
    
    timeOperators = []
    for u in U:
        H = hm.get_hamiltonian(u,t,basis,len(basis))
        timeOperators.append( expm(-1j* H* dtau) )
        
    return timeOperators


#
# takes the next time step and returns all possible values
# Returns a list of values for each time operator
#
# __parameters__
# U       : list of U values
# psi_0   : initial state
# timeOperators : list of time operators
# Heis    : Heisenberg eigenvectors
#
# __returns__
# vals : tuple for each time operator,
#           ( U, probability to be singly occupied,
#             overlap with heisenberg vectors, index )
#
def timestep(U, psi_t0, timeOperators, Heis):
    #vals = dict(zip(range(len(timeOperators)),[np.zeros(timeOperators[0].shape)]*len(timeOperators)))
    vals = []
    for i,evol in enumerate(timeOperators):
        psi_t1  = evol @ psi_t0

        prob_t1 = sum( [vec@psi_t1 * np.conj(vec@psi_t1) for vec in Heis.T] )[0].real
        heis_t1 = (Heis.T[0]@psi_t1 * np.conj(Heis.T[0]@psi_t1))[0].real

        vals.append( (U[i], prob_t1, heis_t1, i) )
        
    return vals

#
# trims the tree by averaging uph value
#
# __parameters__
# root : root of the anytree
# uph  : U, probability, or heisenberg value,
#         {0: U, 1: probability, 2: heisenberg}
#     
def trim_tree(root, uph):
    
    # maximize average prob or heis value
    if uph:
        sorted_leaves = sorted( root.leaves,
                                key= lambda leaf: ( sum([node.vals[uph] for node in leaf.path[1:]]) / len(leaf.path) ),
                                reverse= True )
    # minimize average u value
    else:
        sorted_leaves = sorted( root.leaves,
                                key= lambda leaf: ( sum([node.vals[uph] for node in leaf.path[1:]]) / len(leaf.path) ) )

    # set root to beginning of lowest path
    root = sorted_leaves[0].path[0]
    
    # loop through path and set nodes in tree
    last = root
    for node in sorted_leaves[0].path[1:]:
        
        last.children = [node]
        last = last.children[0]
    
    return

#
# performs a combinatoial search of time evolution operators
# uses anytree to create a tree of search results
# uses a passed function to trim the tree after 5 time steps
# returns the root of the tree
# 
#
# __parameters__
# U       : range of U values to search
# dtau    : delta tau
# psi_0   : starting state (basis vector)
# timeOperators : list of time operators at each U
# Heis    : Heisenberg eigenvectors
#
# __optional parameters__
# uph_first : maximize the first step by u, prob, or heis 
# uph_grow  : grow the tree by u, prob, or heis           { 0: U, 1: Probability, 2: Heisenberg }
# uph_trim  : trim the tree by u, prob, or heis           
#
# min_prob  : minimum probability allowed when searching
# min_heis  : minimum heisenberg overlap allowed when searching
#
# __returns__
# root : the root node of the anytree solution
#
def find_solution(U, dtau, psi_0, timeOperators, Heis, first_step= 'single', grow= 'U', trim= 'U', min_prob= 0.99, min_heis= 0.5):

    params = {
        'U': 0,
        'single': 1,
        'heis': 2
    }

    uph_first = params[first_step]
    uph_grow = params[grow]
    uph_trim = params[trim]

    # tau = 0
    heis = ((Heis[:,0] @ psi_0) * np.conj(Heis[:,0] @ psi_0))[0] # Heisenberg value of starting state
    root = Node( '~', psi= psi_0, tau= dtau, vals= (0,1,heis,-1) )  # root of the tree 
                                                                 # vals= (U, prob, heis, idx)
    
    # find FIRST STEP
    next_vals = timestep(U,psi_0,timeOperators,Heis)
    if uph_first:
        next_vals.sort(key= lambda x: x[uph_first], reverse= True)  # maximize first step
    else:
        next_vals.sort(key= lambda x: x[uph_first])                 # minimize first step (uph_first == 0)
    
    
    # update psi
    psi_0 = timeOperators[ next_vals[0][3] ] @ psi_0
    
    # add to tree
    Node( str(next_vals[0][0]), psi= psi_0, tau= dtau, vals= next_vals[0], parent= root )

    # goes for 100 time steps
    s = timer()
    for i in range(20):

        for j in range(5):

            # loop through leaves
            for sprout in root.leaves:
                
                # test each time operator
                next_vals = timestep(U, sprout.psi, timeOperators, Heis)

                # reduce by constrainsts
                next_poss = [vals for vals in next_vals if vals[1] >= min_prob and vals[2] >= min_heis]

                # sort by U, prob, or heis
                if uph_grow:
                    next_poss.sort(key= lambda x: x[ uph_grow ], reverse= True)
                else:
                    next_poss.sort(key= lambda x: x[ uph_grow ])
                    

                # add next 5 paths
                for five in next_poss[:5]:
                    psi_k = timeOperators[ five[3] ] @ sprout.psi
                    tau_k = sprout.tau + dtau
                    Node( str(five[0]), psi= psi_k, tau= tau_k, vals= five, parent= sprout )
        
        # trim the tree and repeat
        trim_tree(root, uph_trim)

        if i == 0:
            print(f"Estimated time: {round((timer()-s)/60 * 20, 3)} minutes")
        elif i == 19:
            print(f"Actual time: {round((timer()-s)/60, 3)} minutes")
        elif (i+1) % 5 == 0:
            # print timer
            print(f"{(i+1)*5}% done in {round((timer()-s)/60, 3)} minutes")

    return root



#
# plots the tree returned by find_path
# y-axis is determined by uph parameter
# red point indicates U regime, 
#  black indicates t regime
#
# __parameters__
# root : root of tree returned by find_path
# uph  : y-axis (0: U, 1: Probability, 2: Heisenberg)
# 
# __optional parameters__
# x    : range of x (does not start at zero)
#
def plot_solution(root):
    
    path = root.leaves[0].path[1:]
    
    labels = {
        0: ('U', 'U (energy)'),
        1: ('Singly Occupied', 'Probability'),
        2: ('Heisenberg', 'Probability')
    }
    
    fig = plt.figure(figsize=(12,10), layout="tight")
    for i, plot_index in enumerate([221,212,222]):
        
        ax = fig.add_subplot(plot_index)
        
        y = [node.vals[i] for node in path]
        x = np.linspace(root.tau, root.tau*len(y), len(y))
        
        c = ["k" if node.vals[0] < 1 else "r" for node in path]
    
        ax.plot(x,y)
        ax.scatter(x,y, c= c, marker= ".")
        
        ax.set_xlabel("tau")
        ax.set_ylabel(labels[i][1])
        ax.set_title(labels[i][0])
    
    return


def find_solution_parallel(U, dtau, psi_0, timeOperators, Heis, first_step='single', grow='U', trim='U', min_prob= 0.99, min_heis= 0.5):

    params = {
        'U': 0,
        'single': 1,
        'heis': 2
    }

    uph_first = params[first_step]
    uph_grow = params[grow]
    uph_trim = params[trim]

    # tau = 0
    heis = ((Heis[:,0] @ psi_0) * np.conj(Heis[:,0] @ psi_0))[0] # Heisenberg value of starting state
    root = Node( '~', psi= psi_0, tau= dtau, vals= (0,1,heis,-1) )  # root of the tree 
                                                                 # vals= (U, prob, heis, idx)
    
    # find FIRST STEP
    next_vals = te.timestep(U,psi_0,timeOperators,Heis)
    if uph_first:
        next_vals.sort(key= lambda x: x[uph_first], reverse= True)  # maximize first step
    else:
        next_vals.sort(key= lambda x: x[uph_first])                 # minimize first step (uph_first == 0)
    
    
    # update psi
    psi_0 = timeOperators[ next_vals[0][3] ] @ psi_0
    
    # add to tree
    Node( str(next_vals[0][0]), psi= psi_0, tau= dtau, vals= next_vals[0], parent= root )

    # multiprocessing pool
    manager = Manager()
    n_cpus = 4
    pool = Pool(processes=n_cpus)
    print(f'Using {n_cpus} cpus...')

    # goes for 100 time steps
    s = timer()
    for i in range(20):

        for j in range(5):

            # loop through leaves
            for sprout in root.leaves:

                next_vals = []
                k = 0
                for chunk in np.array_split(timeOperators, n_cpus):
                    # test each time operator
                    values = pool.apply_async(te.timestep, args=(U[k:k+len(chunk)], sprout.psi, chunk, Heis))
                    next_vals += values.get()
                    #next_vals = timestep(U[i:i+len(chunk)], sprout.psi, timeOperators, Heis)
                    k += len(chunk)

                # reduce by constrainsts
                next_poss = [vals for vals in next_vals if vals[1] >= min_prob and vals[2] >= min_heis]

                # sort by U, prob, or heis
                if uph_grow:
                    next_poss.sort(key= lambda x: x[ uph_grow ], reverse= True)
                else:
                    next_poss.sort(key= lambda x: x[ uph_grow ])
                    

                # add next 5 paths
                for five in next_poss[:5]:
                    psi_k = timeOperators[ five[3] ] @ sprout.psi
                    tau_k = sprout.tau + dtau
                    Node( str(five[0]), psi= psi_k, tau= tau_k, vals= five, parent= sprout )
        
        # trim the tree and repeat
        te.trim_tree(root, uph_trim)

        if i == 0:
            print(f"Estimated time: {round((timer()-s)/60 * 20, 3)} minutes")
        elif i == 19:
            print(f"Actual time: {round((timer()-s)/60, 3)} minutes")
        elif (i+1) % 5 == 0:
            # print timer
            print(f"{(i+1)*5}% done in {round((timer()-s)/60, 3)} minutes")

    pool.close()
    pool.join()

    return root
