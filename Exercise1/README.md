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
Clone the respository with:
```
git clone https://github.com/SDavenia/FHPC_Assignment.git
```

To compile use the [Makefile](https://github.com/FilippoOlivo/Foundations_of_HPC_Assignment/blob/main/excercise1/Makefile) with:
```
make
```

This will create the executable `main_parallel.exe`.

To run the code use `mpirun` with the following arguments:
- -i: to generate a new random playground
- -r: to evolve the generated playground
- -k *value*: playground size
- -e [0|1|2]: evolution method; 0 means "ordered", 1 means "static", 2 means "white-black"
- -f *string*: name of the .pgm file (the file name should be of the form *init_nnnnn.pgm* where, *nnnnn* is the playground size with 5 digits, padded with zeros.
- -n *int value*: number of steps to evolve the playground
- -s *int value*: every how many steps a snapshot of the playground is saved on a .pgm file. This will be saved in the folder image with name of the form *snapshot_nnnnn.pgm*, where *nnnnn* is the step at which the image is captured padded with zeros to obtain 5 digits.

- - - 
Let us show how an example.

1. To generate an initial playground of size 100:
```
mpirun ./main_parallel.exe -i -k 100 -f init_00100.pgm 
```
2. To evolve the playground in a static way for 10 steps and saving a snapshot every 2 steps using the playground above:
```
mpirun ./main_parallel.exe -r -k 100 -e 1 -f init_00100.pgm -n 10 -s 2
```
- - - 

Note that all the images will be stored in the `images` folder. The generated playgrounds will be stored in `images/initial_matrices`, while the generated snapshots are included in the folder `images/<evolution_method>`, where `evolution_method` is one of the proposed ones.

The `sbatch_files` folder contains the batch scripts to submit the jobs on the ORFEO cluster. Each subfolder contains the script to execute the corresponding evolution type and collect data for each scalability study considered (more detail in the report). To run one of them use:
```
sbatch batch_file.sh
```
