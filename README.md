# Hubbard-Model
Simulating the Hubbard model and dynamically implementing the Gutzwiller projection.

1 HubbardModel func.py
HubbardModel func.py sets up a Hubbard Model, using lists of complex numbers as the possible states (1’s represent up spin and i’s represent down spin).
The function get basis returns the list of all possible states and get hamiltonian
uses this list to return the Hamiltonian. Also present is a helper function for
get hamiltonian, called possible states, it finds the off diagonals of the Hamiltonian. Functions get hamiltonian and possible states suffixed with cyclic, treat
the end sites as neighbors.
H = −t
X
<i,j>,σ
c
†
iσcjσ + U
X
i
ni↓ni↑ − µ
X
i
(ni↓ − ni↑)
Heis =
J
4
X
<i,j>
[2(c
†
i,↑
ci,↓c
†
j,↓
cj,↑ + c
†
i,↓
ci,↑c
†
j,↑
cj,↓ + (ni,↑ − ni,↓)(nj,↑ − nj,↓]
| < ψf |e
−iτ0H0(∆τ)
...e−iτN Hn(∆τ)
|ψi > |
2
2 TimeEvolution.py
Probability is defined as the probability of the system to be in a singly occupied
state. Heisenberg value is defined as the overlap with the Heisenberg eigenvectors.
TimeEvolution.py sets up a half filled Hubbard system with constant spin
(setup model) and conducts a heuristic combinatorial search algorithm (find solution)
to find a combination of time evolution operators that minimize/maximize the
U, probability, or Heisenberg values. The search algorithm creates a tree from
the package anytree then allows this tree to grow by applying all time operators
for the next time step (test timeOps) and taking a step in the direction of the
top five operators that minimize U or maximize probability or heisenberg value;
this choice is indicated by the parameter uph grow. This process is repeated
