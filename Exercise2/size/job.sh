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

srun -n1 make cpu # Now I have all the needed executables.

rm -rf *.csv # Remove previous csv file.

for m_size in 1000 2000
do
    # Run everything with openBLAS double
    for j in 1 2 3 4 5 # Take multiple measurements
    do
        srun -n1 --cpus-per-task=64 ./gemm_oblas_double.x $m_size $m_size $m_size > output.txt #just a temporary file
        # Extract information using grep and regular expressions
        size=$(grep -o 'Size: [0-9]*' output.txt| cut -d' ' -f2)
        times=$(grep -o 'Time: [0-9.]*' output.txt| cut -d' ' -f2)
        gflops=$(grep -o 'GFLOPS: [0-9.]*' output.txt| cut -d' ' -f2)
        # Store the extracted information in a CSV file
        filename="OBLAS_double_$size.csv"

        if [ ! -e $filename ]; then
		echo "Size, Time, GFLOPS" > $filename
        fi
        echo "$size,$times,$gflops" >> $filename
    done

    # Run everything with openBLAS float
    for j in 1 2 3 4 5 # Take multiple measurements
    do
        srun -n1 --cpus-per-task=64 ./gemm_oblas_float.x $m_size $m_size $m_size > output.txt #just a temporary file

        # Extract information using grep and regular expressions
        size=$(grep -o 'Size: [0-9]*' output.txt| cut -d' ' -f2)
        times=$(grep -o 'Time: [0-9.]*' output.txt| cut -d' ' -f2)
        gflops=$(grep -o 'GFLOPS: [0-9.]*' output.txt| cut -d' ' -f2)
        # Store the extracted information in a CSV file
        filename="OBLAS_float_$size.csv"

        if [ ! -e $filename ]; then
            echo "Size,Time,GFLOPS" > $filename
        fi
        echo "$size,$times,$gflops" >> $filename
    done

    # Now repeat for MKL library
    for j in 1 2 3 4 5
    do
        srun -n1 --cpus-per-task=64 ./gemm_mkl_double.x $m_size $m_size $m_size > output.txt #just a temporary file

        # Extract information using grep and regular expressions
        size=$(grep -o 'Size: [0-9]*' output.txt| cut -d' ' -f2)
        times=$(grep -o 'Time: [0-9.]*' output.txt| cut -d' ' -f2)
        gflops=$(grep -o 'GFLOPS: [0-9.]*' output.txt| cut -d' ' -f2)
        # Store the extracted information in a CSV file
        filename="MKL_double_$size.csv"

        if [ ! -e $filename ]; then
            echo "Size,Time,GFLOPS" > $filename
        fi
        echo "$size,$times,$gflops" >> $filename
    done

    # Now repeat for MKL library float
    for j in 1 2 3 4 5
    do
        srun -n1 --cpus-per-task=64 ./gemm_mkl_float.x $m_size $m_size $m_size > output.txt #just a temporary file

        # Extract information using grep and regular expressions
        size=$(grep -o 'Size: [0-9]*' output.txt| cut -d' ' -f2)
        times=$(grep -o 'Time: [0-9.]*' output.txt| cut -d' ' -f2)
        gflops=$(grep -o 'GFLOPS: [0-9.]*' output.txt| cut -d' ' -f2)
        # Store the extracted information in a CSV file
        filename="MKL_float_$size.csv"

        if [ ! -e $filename ]; then
            echo "Size,Time,GFLOPS" > $filename
        fi
        echo "$size,$times,$gflops" >> $filename
    done
done
