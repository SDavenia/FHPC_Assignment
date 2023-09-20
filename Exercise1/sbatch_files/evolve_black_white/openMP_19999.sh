#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=OMP_BlackWhite
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=128
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=OMP_BlackWhite_19999.out

module load openMPI/4.1.5/gnu

cd ../../ 
make
# cd sbatch_files/initialization

MAPBY=socket

export OMP_PLACES=cores
export OMP_PROC_BIND=close

nsteps=100
s=`expr $nsteps + 2`
# To ensure printing times are not included in the evolve times

out_filename=results/evolve_black_white/openMP_19999.csv # To write times
echo "size,threads,time" > $out_filename

for ksize in 19999
do
    filename="init_19999.pgm" # To write image
    for n_threads in 1 2 4 8 16 32 64 128
    do
        export OMP_NUM_THREADS=$n_threads
        for j in {1..5..1}
        do
            mpirun -n 1 --map-by $MAPBY ./main_parallel.exe -r -k $ksize -e 2 -f $filename -n $nsteps -s $s  > output_black_white_openMP_19999.txt
            time_value=$(grep -o 'BlackWhite time: [0-9.]*' output_black_white_openMP_19999.txt | awk '{print $3}')
            echo "$ksize,$n_threads,$time_value" >> $out_filename
        done
    done
done

rm output_black_white_openMP_19999.txt # Remove useless temporary file


