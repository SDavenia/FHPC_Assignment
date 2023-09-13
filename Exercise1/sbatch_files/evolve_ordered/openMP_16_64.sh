#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=OMP_Ordered
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2 
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=OMP_Ordered.out

module load openMPI/4.1.5/gnu

cd ../../ 
make
# cd sbatch_files/initialization

# Define MPI binding and OMP affinity
MAPBY=node
BINDTO=socket

export OMP_PLACES=cores
export OMP_PROC_BIND=close

nsteps=100
s=`expr $nsteps + 2`
# To ensure printing times are not included in the evolve times

out_filename=results/evolve_ordered/openMP.csv # To write times
echo "size,threads,time" > $out_filename

for ksize in 10000
do
    formatted_number=$(printf "%05d" "$ksize")
    filename="init_"$formatted_number".pgm" # To write image
    for n_threads in 16 32 64
    do
        export OMP_NUM_THREADS=$n_threads
        for j in {1..5..1}
        do
            mpirun -n 4 --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -r -k $ksize -e 1 -f $filename -n $nsteps -s $s  > output_ordered_openMP.txt
            time_value=$(grep -o 'Ordered time: [0-9.]*' output_ordered_openMP.txt | awk '{print $3}')
            echo "$ksize,$n_threads,$time_value" >> $out_filename
        done
    done
done

rm output_ordered_openMP.txt # Remove useless temporary file


