# Exercise 2
This folder contains code and results for the Exercise 2. More details on the requirements can be found [here](https://github.com/Foundations-of-HPC/Foundations_of_HPC_2022/tree/main/Assignment/exercise2).

The aim of this exercise is to compare the performance of three math libraries which are used on HPC: openBLAS, MKL and BLIS. This latter was downloaded and compiled by us.
To be more specific the task is to compare the level 3 BLAS function called *gemm*, which performs matrix multiplication in an efficient way.

## Structure of this repository
The exercise requires to run *gemm* on both the AMD Epyc nodes available on the cluster and also the INTEL Thin nodes. 
As such there are two different folders that contain measurements from that nodes:
- [`EPYC`](https://github.com/SDavenia/FHPC_Assignment/tree/main/Exercise2/EPYC)
- [`THIN`](https://github.com/SDavenia/FHPC_Assignment/tree/main/Exercise2/THIN)

Inside each folder there are two additional folders:
- `cores` which contains batch files and results related to *cores scalability*
- `size` which contains batch files and results related to *size scalability*

In both cases there are two batch files `job_close.sh` and `job_spread.sh` to run the code using different thread affinity policies, respectively *close* and *spread*.

As such the `.csv` files with the results are divided in two different folders `close` and `spread`. Each of the resulting files contains final specifications in the name:
- In folder `cores`, each file is named as \<library>\_\<precision>\_\<number of cores\>.csv 
- In folder `size`, each file is named as \<library>\_\<precision>\_\<matrix size\>.csv

Where \<library> specifies the math library being used (see above) and \<precision> specifies whether double or float precision was used.

Finally, folder `parallel_initialization` contains a modified `gemm_modified.c` files where we use a parallel initialization of the matrices.

## Additional files
- There are some `analysis.ipynb` files which are used to extract and plot the results which are reported in the report.
- The folder `attempts` contains some old results and previous trials, and should be ignored.
- The `MAKEFILE` and `gemm.c` files were provided by the course teachers, and only slightly modified.

--------------------------------------------------------------------------------
As such the `EPYC` folder contains measurements from the first, while the folder `THIN` contains measurements on the latter.
Additionally the exercise required to take measurements under these conditions:
 - For a fixed matrix size increase the number of cores available (results can be found in folder *cores*)
 - For a fixed number of cores increase the matrix size (results can be found in folder *size*)

Finally we also compared different thread affinity policies, respectively CLOSE and SPREAD.
The results from each of these settings are contained in the corresponding folders.

To conclude, each of the resulting files contains final specifications in the name. 
- In folder *cores*, each file is named as <library>_<precision>_<number of cores>.csv 
- In folder *size*, each file is named as <library>_<precision>_<matrix size>.csv
Where library specifies the math library being used (see above) and precision specifies whether double or float precision was used.

As a final note on the files contained in this folder:
- There are some `.sh` files, which were used to submit the jobs to the computer cluster, which uses SLURM as a workload manager.
- The folder `attempts` contains some old results and previous trials, and should be ignored.
- There are some `analysis.ipynb` files which are used to extract and plot the results which are reported in the report.
- The `MAKEFILE` and `gemm.c` files were provided by the course teachers, and only slightly modified.
