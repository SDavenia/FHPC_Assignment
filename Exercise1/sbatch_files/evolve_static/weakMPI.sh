#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=WMPIS_Static
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=2 
#SBATCH --cpus-per-task=64 
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=WMPIS_Static.out

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

out_filename=results/evolve_static/weak_MPI.csv # To write times
echo "size,processes,time" > $out_filename

initial_size=10000

for n_processes in {1..6..1}
do  
    current_size=$(echo "scale=0; sqrt($n_processes * $initial_size^2)" | bc)
    formatted_number=$(printf "%05d" "$current_size")
    filename="init_"$formatted_number".pgm" # To write image
    for j in {1..5..1}
    do
        mpirun -n $n_processes --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -r -k $current_size -e 1 -f $filename -n $nsteps -s $s > output_static_strong_MPI_socket.txt
        time_value=$(grep -o 'Static time: [0-9.]*' output_static_weak_MPI.txt | awk '{print $3}')
        echo "$current_size,$n_processes,$time_value" >> $out_filename
    done
    

done
rm output_static_weak_MPI.txt # Remove useless temporary file



