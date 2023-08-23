#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=gemm_first_attempt
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=128
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --output=output.out

module load architecture/AMD
module load mkl
module load openBLAS/0.3.23-omp

export LD_LIBRARY_PATH=/u/dssc/sdaven00/myblis/lib:$LD_LIBRARY_PATH

export OMP_PLACES=cores
export OMP_PROC_BIND=close

srun -n1 make cpu # Now I have all the needed executables.

m_size=10000


for implem in 'oblas' 'mkl' 'blis'
do
    for type in 'double' 'float'
    do
        for n_threads in {2..128..2}
        do
            export OMP_NUM_THREADS=$n_threads
            export BLIS_NUM_THREADS=$n_threads 

            # Run everything with openBLAS double
            for j in 1 2 3 4 5 # Take multiple measurements
            do
                srun -n1 --cpus-per-task=64 ./gemm_"$implem"_"$type".x $m_size $m_size $m_size > output.txt #just a temporary file
                
                # Extract information using grep and regular expressions
                times=$(grep -o 'Time: [0-9.]*' output.txt| cut -d' ' -f2)
                gflops=$(grep -o 'GFLOPS: [0-9.]*' output.txt| cut -d' ' -f2)
                # Store the extracted information in a CSV file
                filename=close/"$implem"_"$type"_$n_threads.csv

                if [ ! -e $filename ]; then
                echo "n_threads,Time,GFLOPS" > $filename
                fi
                echo "$n_threads,$times,$gflops" >> $filename
            done
        done
    done
done
rm output.txt # Delete the temporary file
