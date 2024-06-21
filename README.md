# Alternating algorithms for cross approximation
 
This repository contains the code to reproduce the numerical experiments of my master thesis `Alternating algorithms for cross approximation` carried out at EPFL.

## Repository description
This repository is composed of two main repositories `Matrices` and `Tensors` whicn contain respectively the code to reproduce the numerical experiments of the matrix case and the tensor case. They are structured in the similar following way :
- Every file starting by "script" is a script to reproduce some numerical experiments.
- All the other files are functions used to initialize and compute the different numerical experiments. Notice that the algorithm "Iterative Refinement" is coded several times to optimize the time depending on the experiment done.

## Reproducibility of the results
To reproduce the numerical experiments of my master thesis, you must run the files starting with "script" and set correclty the parameters. Note that you can only do an experiment for one example at a time but you may run several initializations. To enjoy the full speed up gained by using parfor loops, the "Parallel Computing Toolbox" must be installed. Remark that it is not required to install it in order to run the scripts.

## Author
Valentin Comment
