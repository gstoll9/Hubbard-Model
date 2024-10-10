# Hubbard-Model
Simulating the Hubbard model and dynamically implementing the Gutzwiller projection.

## Hubbard_func.py

HubbardModel func.py
HubbardModel func.py sets up a Hubbard Model, using lists of complex num- bers as the possible states (1’s represent up spin and i’s represent down spin). The function get basis returns the list of all possible states and get hamiltonian uses this list to return the Hamiltonian. Also present is a helper function for get hamiltonian, called possible states, it finds the off diagonals of the Hamilto- nian. Functions get hamiltonian and possible states suffixed with cyclic, treat the end sites as neighbors.
$$
H = −t\ \sum_{<i,j>,\sigma} c^\dagger_{i\sigma} c_{j\sigma} + U \sum_i n_{i\downarrow} n_{i\uparrow} − \mu \sum_j ( n_{i\downarrow} − n_{i\uparrow} )
$$$$
Heis = \frac{J}{4}\ \sum_{<i,j>} \big[ 2 ( c^\dagger_{i,\uparrow} c_{i,\downarrow} c^\dagger_{j,\downarrow} c_{j,\uparrow} + c^\dagger_{i,\downarrow} c_{i,\uparrow} c^\dagger_{j,\uparrow} c_{j,\downarrow} + ( n_{i,\uparrow} − n_{i,\downarrow} ) ( n_{j,\uparrow} − n_{j,\downarrow} \big]
<i,j >
$$$$
| < \psi_f |e^{−i\tau_0H_0(\Delta\tau)}...e^{−i\tau_NH_n(\Delta\tau)}|\psi_i > |^2
$$

## TimeEvolution.py
_Probability is defined as the probability of the system to be in a singly occupied state. Heisenberg value is defined as the overlap with the Heisenberg eigenvectors._\
TimeEvolution.py sets up a half filled Hubbard system with constant spin (setup model) and conducts a heuristic combinatorial search algorithm (find solution) to find a combination of time evolution operators that minimize/maximize the U, probability, or Heisenberg values. The search algorithm creates a tree from the package anytree then allows this tree to grow by applying all time operators for the next time step (test timeOps) and taking a step in the direction of the top five operators that minimize U or maximize probability or heisenberg value; this choice is indicated by the parameter uph grow. This process is repeated five times and then the tree is trimmed (reduced) to the solution that has the lowest average U value or highest average probability or Heisenberg value; this choice is indicated by the parameter uph trim. This solution can then be plot- ted with plot solution with the parameter, uph, indicating U, probability, or Heisenberg as the y-axis. There are also functions to get a list of time operators (get timeOps) and to get the Heisenberg matrix (get heisenberg) with helper function, exchanged states.

## ExampleScript.ipynb
Sets up a model then shows different ways to use the find_solution function and plot the results.

## Dynamically Implementing the Gutzwiller Projection.ppt
A presentation summarizing my first semester of research. I was focused on presenting my results and showing how this research overlapped with the Computational Modeling and Data Analytics degree. The code in this repository is talked about at the end of the presentation. We found an error in our results from the first semester, but the plots still look cool!
