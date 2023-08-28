# Exercise 2
The aim of this exercise is to compare the performance of three math libraries which are used on HPC. 
These are openBLAS, MKL and BLIS. This latter was downloaded and compiled by us.
To be more specific the task is to compare the level 3 BLAS function called *gemm*, which performs matrix multiplication in an efficient way.

The exercise requires to run this on both the AMD Epyc nodes available on the cluster and also the INTEL Thin nodes. 
As such the *EPYC* folder contains measurements from the first, while the folder *THIN* contains measurements on the latter.
Additionally the exercise required to take measurements under these conditions:
 - For a fixed matrix size increase the number of cores available (results can be found in folder *cores*)
 - For a fixed number of cores increase the matrix size (results can be found in folder *size*)

Finally we also compared different thread affinity policies, respectively the DEFAULT one (without specifying anything), CLOSE and SPREAD. 
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
