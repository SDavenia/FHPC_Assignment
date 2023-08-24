#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=gemm_first_attempt
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb 
#SBATCH --time=01:00:00 
#SBATCH --output=output.out
module load architecture/AMD
module load mkl
module load openBLAS/0.3.23-omp

srun -n1 make cpu

for size in {2000..6000..500}
do
    output_file_name="size_$size"
    echo "Running using matrix of size $size" > "$output_file_name.txt"
    echo "MKL library" >> "$output_file_name.txt"
    srun -n1 --cpus-per-task=64 ./gemm_mkl.x $size $size $size >> "$output_file_name.txt"
    echo "BLAS library" >> "$output_file_name.txt"
    srun -n1 --cpus-per-task=64 ./gemm_oblas.x $size $size $size >> "$output_file_name.txt"
    echo ---------------------------------------------------------------------------------------------------------  
done
