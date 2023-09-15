#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=WMPIS_Init
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=2 
#SBATCH --cpus-per-task=64 
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=WMPIS_Init.out

module load openMPI/4.1.5/gnu

cd ../../ 
make
# cd sbatch_files/initialization

# Define MPI binding and OMP affinity
MAPBY=node
BINDTO=socket

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=64

out_generate=results/initialization/weak_MPI_generate.csv # To write times
out_write=results/initialization/weak_MPI_write.csv # To write times
echo "size,processes,time" > $out_generate
echo "size,processes,time" > $out_write

initial_size=10000

for n_processes in {1..6..1}
do  
    current_size=$(echo "scale=0; sqrt($n_processes * $initial_size^2)" | bc)
    formatted_number=$(printf "%05d" "$current_size")
    filename="init_"$formatted_number".pgm" # To write image
    for j in {1..5..1}
    do
        mpirun -n $n_processes --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -i -k $current_size -f $filename > output_initialization_weak_MPI.txt
        time_value=$(grep -o 'Generate time: [0-9.]*' output_initialization_weak_MPI.txt | awk '{print $3}')
        echo "$current_size,$n_processes,$time_value" >> $out_generate
        time_value=$(grep -o 'Write time: [0-9.]*' output_initialization_weak_MPI.txt | awk '{print $3}')
        echo "$current_size,$n_processes,$time_value" >> $out_write
    done
    

done
rm output_initialization_weak_MPI.txt # Remove useless temporary file



