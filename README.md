# Hubbard-Model
Simulating the Hubbard model and dynamically implementing the Gutzwiller projection.

## Example_Script.ipynb
Sets up a model then shows different ways to use the find_solution function and plot the results.

## Dynamically Implementing Gutzwiller Projection.pptx
A presentation summarizing my first semester of research. I was focused on presenting my results and showing how this research overlapped with the Computational Modeling and Data Analytics degree. The code in this repository is talked about at the end of the presentation. We found an error in our results from that semester, but the plots still look cool!

## Hubbard_func.py
HubbardModel func.py
HubbardModel func.py sets up a Hubbard Model, using lists of complex numbers as the possible states (1’s represent up spin and i’s represent down spin). The function _get\_basis_ returns the list of all possible states and _get\_hamiltonian_ uses this list to return the Hamiltonian. Also present is a helper function for _get\_hamiltonian_, called _possible\_states_, it finds the off diagonals of the Hamiltonian. Functions _get\_hamiltonian_ and _possible\_states_ suffixed with _\_cyclic_, treat the end sites as neighbors.

## TimeEvolution.py
_Probability is defined as the probability of the system to be in a singly occupied state. Heisenberg value is defined as the overlap with the Heisenberg eigenvectors._\
TimeEvolution.py has functions to set up a half filled Hubbard system with constant spin (_setup\_model_) and conduct a heuristic combinatorial search algorithm (_find\_solution_) to find a combination of time evolution operators that minimize U or maximize probability or Heisenberg value. The search algorithm creates a tree from the package _anytree_ then allows this tree to grow by applying all time operators for the next time step (_timestep_) and taking a step in the direction of the top five operators that minimize U or maximize probability or Heisenberg value; this choice is indicated by the parameter _grow_. This process is repeated five times and then the tree is trimmed to the solution that has the lowest average U value or highest average probability or Heisenberg value; this choice is indicated by the parameter _trim_. This solution can then be plotted with _plot\_solution_ which displays U, probability, and Heisenber value over tau. There are also functions to get a list of time operators (_get\_timeOps_) and to get the Heisenberg matrix (_get\_heisenberg_) with helper function, _exchanged\_states_.


