# Hubbard-Model
Simulating the Hubbard model and dynamically implementing the Gutzwiller projection.

## Example_Script.ipynb
Sets up a model then shows different ways to use the find_solution function and plot the results.

## Dynamically Implementing Gutzwiller Projection.pptx
A presentation summarizing my first semester of research. I was focused on presenting my results and showing how this research overlapped with the Computational Modeling and Data Analytics degree. The code in this repository is talked about at the end of the presentation. We found an error in our results from the first semester, but the plots still look cool!

## Hubbard_func.py
HubbardModel func.py
HubbardModel func.py sets up a Hubbard Model, using lists of complex num- bers as the possible states (1’s represent up spin and i’s represent down spin). The function get basis returns the list of all possible states and get hamiltonian uses this list to return the Hamiltonian. Also present is a helper function for get hamiltonian, called possible states, it finds the off diagonals of the Hamilto- nian. Functions get hamiltonian and possible states suffixed with cyclic, treat the end sites as neighbors.

## TimeEvolution.py
_Probability is defined as the probability of the system to be in a singly occupied state. Heisenberg value is defined as the overlap with the Heisenberg eigenvectors._\
TimeEvolution.py sets up a half filled Hubbard system with constant spin (setup model) and conducts a heuristic combinatorial search algorithm (find solution) to find a combination of time evolution operators that minimize/maximize the U, probability, or Heisenberg values. The search algorithm creates a tree from the package anytree then allows this tree to grow by applying all time operators for the next time step (test timeOps) and taking a step in the direction of the top five operators that minimize U or maximize probability or heisenberg value; this choice is indicated by the parameter uph grow. This process is repeated five times and then the tree is trimmed (reduced) to the solution that has the lowest average U value or highest average probability or Heisenberg value; this choice is indicated by the parameter uph trim. This solution can then be plot- ted with plot solution with the parameter, uph, indicating U, probability, or Heisenberg as the y-axis. There are also functions to get a list of time operators (get timeOps) and to get the Heisenberg matrix (get heisenberg) with helper function, exchanged states.


