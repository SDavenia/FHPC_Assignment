#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=SMPIS_Init
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2 
#SBATCH --cpus-per-task=16 
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=SMPIS_Init.out

module load openMPI/4.1.5/gnu

cd ../../ 
make
# cd sbatch_files/initialization

# Define MPI binding and OMP affinity
MAPBY=node
BINDTO=socket

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16

out_generate=results/initialization/strong_MPI_socket_generate.csv # To write times
out_write=results/initialization/strong_MPI_socket_write.csv # To write times
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
            mpirun -n $n_processes --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -i -k $ksize -f $filename > output_initialization_strong_MPI_socket.txt
            time_value=$(grep -o 'Generate time: [0-9.]*' output_initialization_strong_MPI_socket.txt | awk '{print $3}')
            echo "$ksize,$n_processes,$time_value" >> $out_generate
            time_value=$(grep -o 'Write time: [0-9.]*' output_initialization_strong_MPI_socket.txt | awk '{print $3}')
            echo "$ksize,$n_processes,$time_value" >> $out_write
        done
    done
done

rm output_initialization_strong_MPI_socket.txt # Remove useless temporary file


