# Exercise 2
This folder contains code and results for Exercise 2. More details on the requirements can be found [here](https://github.com/Foundations-of-HPC/Foundations_of_HPC_2022/tree/main/Assignment/exercise2).

The aim of this exercise is to compare the performance of three math libraries which are used on HPC: openBLAS, MKL and BLIS. The latter was downloaded and compiled by us.
To be more specific the task is to compare the level 3 BLAS function called *gemm*, which performs matrix multiplication in an efficient way.

## Structure of this repository
The exercise requires to run [`gemm.c`](https://github.com/SDavenia/FHPC_Assignment/blob/main/Exercise2/gemm.c) on both the AMD Epyc and the Intel Thin nodes available on the cluster used. 
As such there are two different folders that contain measurements from that nodes:
- [`EPYC`](https://github.com/SDavenia/FHPC_Assignment/tree/main/Exercise2/EPYC)
- [`THIN`](https://github.com/SDavenia/FHPC_Assignment/tree/main/Exercise2/THIN)

Contained inside each folder there are two additional folders:
- `cores` which contains batch scripts and results related to *cores scalability*
- `size` which contains batch scripts and results related to *size scalability*

In both cases there are two batch scripts `job_close.sh` and `job_spread.sh` to run the code using different thread affinity policies, respectively *close* and *spread*.

As such the `.csv` files with the results are divided in two additional subfolders `close` and `spread`. Each of the resulting files contains the options being used in the name:
- In folder `cores`, each file is named as \<library>\_\<precision>\_\<number of cores\>.csv 
- In folder `size`, each file is named as \<library>\_\<precision>\_\<matrix size\>.csv

Where \<library> specifies the math library being used (see above) and \<precision> specifies whether double or single precision was used.

Finally, folder `parallel_initialization` contains a modified `gemm_modified.c` file, where unlike the original code the matrices were initialized using OpenMP threads. It also contains results obtained with the same structure as the original directory.

## Additional files
- There are some `analysis.ipynb` files which are used to extract and plot the results which are shown in the report.
- The folder `attempts` contains some old results and previous trials, and should be ignored.
- The `MAKEFILE` and `gemm.c` files were provided by the course teachers, and only slightly modified.

## How to run
1. Decide what scalability you want to analyze and the type of node. Then enter the corresponding folder `<node_type>/<scalability_type>`. 
2. Execute the `job_close.sh` or `job_spread.sh` batch scripts using for example:
```
sbatch job_close.sh
```
This will use the Makefile specific to the problem under investigation, which is located in the same directory.

3. All results will be written on .csv files as described above.
