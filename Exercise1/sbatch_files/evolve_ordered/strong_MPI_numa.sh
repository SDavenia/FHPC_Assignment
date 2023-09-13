#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=SMPIN_Ordered
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2 
#SBATCH --cpus-per-task=16 
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=SMPIN_Ordered.out

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

out_filename=results/evolve_ordered/strong_MPI_numa.csv # To write times
echo "size,processes,time" > $out_filename

nsteps=100
s=`expr $nsteps + 2`

for ksize in 10000
do
    formatted_number=$(printf "%05d" "$ksize")
    filename="init_"$formatted_number".pgm" # To write image
    for n_processes in 1 2 4 8 16
    do
        for j in {1..5..1}
        do
            mpirun -n $n_processes --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -r -k $ksize -e 1 -f $filename -n $nsteps -s $s > output_ordered_strong_MPI_numa.txt
            time_value=$(grep -o 'Ordered time: [0-9.]*' output_ordered_strong_MPI_numa.txt | awk '{print $3}')
            echo "$ksize,$n_processes,$time_value" >> $out_filename
        done
    done
done

rm output_ordered_strong_MPI_numa.txt # Remove useless temporary file



# mpirun -n 4 --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -r -k $ksize -e 1 -f $filename -n $nsteps -s $s