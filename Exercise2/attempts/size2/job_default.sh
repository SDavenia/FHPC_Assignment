#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=gemm_first_attempt
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --output=output.out

module load architecture/AMD
module load mkl
module load openBLAS/0.3.23-omp

export LD_LIBRARY_PATH=/u/dssc/sdaven00/myblis/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=64
export BLIS_NUM_THREADS=64 

srun -n1 make cpu # Now I have all the needed executables.


for implem in 'oblas' 'mkl' 'blis'
do
    for type in 'double' 'float'
    do
        for m_size in {2000..3000..1000}
        do
            for j in 1 2 3 4 5 # Take multiple measurements
            do
                srun -n1 --cpus-per-task=64 ./gemm_"$implem"_"$type".x $m_size $m_size $m_size > output.txt #just a temporary file
                # Extract information using grep and regular expressions
                size=$(grep -o 'Size: [0-9]*' output.txt| cut -d' ' -f2)
                times=$(grep -o 'Time: [0-9.]*' output.txt| cut -d' ' -f2)
                gflops=$(grep -o 'GFLOPS: [0-9.]*' output.txt| cut -d' ' -f2)
                # Store the extracted information in a CSV file
                filename=default/"$implem"_"$type"_$size.csv

                if [ ! -e $filename ]; then
                echo "Size,Time,GFLOPS" > $filename
                fi
                echo "$size,$times,$gflops" >> $filename
            done
        done
    done
done
rm output.txt # Delete the temporary file
