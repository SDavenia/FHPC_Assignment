# Exercise 1
This folder contains code and results for the Exercise 1. More details on the requirements can be found [here](https://github.com/Foundations-of-HPC/Foundations_of_HPC_2022/blob/main/Assignment/exercise1/Assignment_exercise1.pdf).\
The goal of the exercise is to implement a parallel version of a variant of the famous [Conway’s “Game of Life”](https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life).\
Different evolution methods are implemented:
- Static Evolution
- Ordered Evolution
- White-Black evolution

## Structure of the repository
In the `src` folder you can find the code for the different evolution methods we implemented:
- [Static evolution](https://github.com/SDavenia/FHPC_Assignment/blob/main/Exercise1/src/evolve_static.c)
- [Ordered evolution](https://github.com/SDavenia/FHPC_Assignment/blob/main/Exercise1/src/evolve_ordered.c)
- [White-Black evolution](https://github.com/SDavenia/FHPC_Assignment/blob/main/Exercise1/src/black_white.c)
- [Read, Write and random generation of the matrix](https://github.com/SDavenia/FHPC_Assignment/blob/main/Exercise1/src/read_write_parallel.c)

In the `results` folder you can find the .csv files with the times obtained by running the code. The folder is subdivided according to the evolution methods used.

## How to run this code
Clone the respository with:\
```git clone https://github.com/SDavenia/FHPC_Assignment.git```\

In the `sbatch_files` folder there are the batch files to run the code. To run one of the evolution method, enter the corresponding folder and use:\
```sbatch evolution_method.sh```
