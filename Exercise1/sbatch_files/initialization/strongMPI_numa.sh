#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=SMPIN_Init
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16 
#SBATCH --mem=200gb 
#SBATCH --time=00:45:00 
#SBATCH --exclusive
#SBATCH --output=SMPIN_Init.out

module load openMPI/4.1.5/gnu

cd ../../ 
make
# cd sbatch_files/initialization

# Define MPI binding and OMP affinity
MAPBY=node
BINDTO=numa

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16

out_generate=results/initialization/strong_MPI_numa_generate_10000.csv # To write times
out_write=results/initialization/strong_MPI_numa_write_10000.csv # To write times
echo "size,processes,time" > $out_generate
echo "size,processes,time" > $out_write

for ksize in 10000 20000
do
    formatted_number=$(printf "%05d" "$ksize")
    filename="init_"$formatted_number".pgm" # To write image
    for n_processes in {1..8..1}
    do
        for j in {1..5..1}
        do
            mpirun -n $n_processes --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -i -k $ksize -f $filename > output_initialization_strong_MPI_numa.txt
            
            time_value=$(grep -o 'Generate time: [0-9.]*' output_initialization_strong_MPI_numa.txt | awk '{print $3}')
            echo "$ksize,$n_processes,$time_value" >> $out_generate
            time_value=$(grep -o 'Write time: [0-9.]*' output_initialization_strong_MPI_numa.txt | awk '{print $3}')
            echo "$ksize,$n_processes,$time_value" >> $out_write
        done
    done
    for n_processes in {10..16..2}
    do
        for j in {1..5..1}
        do
            mpirun -n $n_processes --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -i -k $ksize -f $filename > output_initialization_strong_MPI_numa.txt
            
            time_value=$(grep -o 'Generate time: [0-9.]*' output_initialization_strong_MPI_numa.txt | awk '{print $3}')
            echo "$ksize,$n_processes,$time_value" >> $out_generate
            time_value=$(grep -o 'Write time: [0-9.]*' output_initialization_strong_MPI_numa.txt | awk '{print $3}')
            echo "$ksize,$n_processes,$time_value" >> $out_write
        done
    done
done

rm output_initialization_strong_MPI_numa.txt # Remove useless temporary file


