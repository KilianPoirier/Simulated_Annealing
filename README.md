# Simulated Annealing

This project is about solving the number partitioning problem using simulated annealing. Simulated annealing is a classical metaheuristic used to find approximate solutions to optimisation problems in a large search space. It uses stochastic sampling to mimic the way materials cool down to more stable and ordered configuration.
The proposed solution is written in python.
Two files make up the project: in `SA_toolbox.py` you will find the functions needed to build the algorithm, as well as the algorithm itself. In `SA_results.ipynb`, there are a set of plots to interpret the results and performance of the algorithm.

The libraries used in this project are the following:

`numpy` : numerical library used to make computations.

`matplotlib.pyplot` : graphic library for visualization.

`time` : library used to track time performance.

## Overall structure
Our approach of simulated annealing highly relies on the definition of the Ising Hamiltonian (definition is given in the results notebook). It also makes use of the metropolis criterion (definition in the notebook) that contributes to the introduction of randomness in the search. Using this criterion, the search does not get stuck in local minima. This criterion depends on a hyperparameter (temperature) that controls how much random changes are accepted.

### Functions
The functions must be imported in the notebook using `import SA_toolbox *`

`compute_energy`: Computes the value of the Ising Hamiltonian of a partition.

`test_metropolis_crit`: Tests for the metropolis criterion.

`get_temperature`: Computes the temperature depending on the time (selected between three schemes).

`init_random_candidate`: Initialize the partition with a random partition.

`find_candidate`: Find a new candidate partition by switching one number from one set to the other.

`record_energy`: Function used to record the energy for discussion.

`simulate_loop_annealing`: Implementation of the simulated annealing algorithm (details in the algorithm section)

### Description of the algorithm
The program begins with a random partition of the number set. A first candidate is found using `init_random_candidate`. The temperature is initialized at T0. A candidate partition is found using `find_candidate` and we compute the difference between the current and candidate partition's energies. If the candidate energy is lower, the partition is updated. If not, the metropolis criterion is computed, it is based on Boltzmann thermal distribution. The lower the temperature, the lower the criterion gets. A random probability is compared to the criterion, if the criterion is larger, the partition is accepted. This operation is then looped and the temperature decreases at each iteration.

## Comments
The approach presented shows good results for well chosen parameters (k and T0). By setting a temperature high enough and a cooling scheme slow enough, the simulated annealing algorithm manages to find the global minimum. When the minimum is not the expected zero energy, the algorithm converges to a statiscally optimal value.
