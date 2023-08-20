#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=gemm_first_attempt
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb 
#SBATCH --time=00:05:00 
#SBATCH --output=attempt.csv
for size in 2000 4000 6000 8000 10000 12000 14000 16000 18000 20000
do
    echo "Running using matrix of size $size"
    srun -n1 --cpu-per-task=64 ./gemm_mkl.x $size $size $size 
    srun -n1 --cpu-per-task=64 ./gemm_oblas.x $size $size $size 
done