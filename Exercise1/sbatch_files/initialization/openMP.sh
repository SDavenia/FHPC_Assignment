#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=OMP_Init
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2 
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb 
#SBATCH --time=00:30:00 
#SBATCH --exclusive
#SBATCH --output=OMP_Init.out

module load openMPI/4.1.5/gnu

cd ../../ 
make

# Define MPI binding and OMP affinity
MAPBY=node
BINDTO=socket

export OMP_PLACES=cores
export OMP_PROC_BIND=close

out_write=results/initialization/openMP_write.csv
out_generate=results/initialization/openMP_generate.csv

echo "size,threads,time" > $out_write
echo "size,threads,time" > $out_generate

for ksize in 20000
do
    formatted_number=$(printf "%05d" "$ksize")
    filename="init_"$formatted_number".pgm" # To write image
    for n_threads in 1 2 4 
    do
        export OMP_NUM_THREADS=$n_threads
        for j in {1..5..1}
        do
            mpirun -n 4 --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -i -k $ksize -f $filename > output_initialization_openMP.txt
            
            time_value=$(grep -o 'Generate time: [0-9.]*' output_initialization_openMP.txt | awk '{print $3}')
            echo "$ksize,$n_threads,$time_value" >> $out_generate
            time_value=$(grep -o 'Write time: [0-9.]*' output_initialization_openMP.txt | awk '{print $3}')
            echo "$ksize,$n_threads,$time" >> $out_write
        done
    done

    for n_threads in {8..64..4}
    do
        export OMP_NUM_THREADS=$n_threads
        for j in {1..5..1}
        do
            mpirun -n 4 --map-by $MAPBY --bind-to $BINDTO ./main_parallel.exe -i -k $ksize -f $filename > output_initialization_openMP.txt
            
            time_value=$(grep -o 'Generate time: [0-9.]*' output_initialization_openMP.txt | awk '{print $3}')
            echo "$ksize,$n_threads,$time_value" >> $out_generate
            time_value=$(grep -o 'Write time: [0-9.]*' output_initialization_openMP.txt | awk '{print $3}')
            echo "$ksize,$n_threads,$time_value" >> $out_write
        done
    done
done

rm output_initialization_openMP.txt # Remove useless temporary file


