#define _GNU_SOURCE // sched_getcpu(3) is glibc-specific (see the man page)

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <sched.h> 


int main(int argc, char* argv[]){

    int k; // Matrix size
    if (argc > 1)
        k = atoi(argv[1]);
    else
        k = 5;
    //printf("Initialising matrix of size %d\n", k);

    int rank, size;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD,&rank );
    MPI_Comm_size( MPI_COMM_WORLD,&size );  

    // Defines how many rows each process has to initialise
    int rows_read = k / size; 
    if (rank < k%size) // For remainder 
        rows_read += 1;

    char image_name[] = "init.pgm";
    unsigned char ptr[rows_read * k];

    MPI_Offset disp;
    MPI_File   fh;
    MPI_File_open(  MPI_COMM_WORLD, image_name, 
                    MPI_MODE_RDONLY,
                    MPI_INFO_NULL, &fh  );
    
    if (rank >= k % size)
        rows_read += k % size;

    disp = rank * rows_read * k *sizeof(unsigned char) + 64; // 64 in this case is the hardcoded representation of the header of the images

    MPI_File_seek(fh, disp, MPI_SEEK_CUR);
    MPI_File_read_all(fh, ptr, rows_read*k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);
    MPI_Barrier(MPI_COMM_WORLD);

    FILE* prova_file;
    char nome_file[] = "prova.txt";
    prova_file = fopen(nome_file, "w");
    if(rank==0){
        fprintf(prova_file,"I am process %d\n", rank);
        for(int i = 0; i < 3*5; i++)
            fprintf(prova_file, "%u ",ptr[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    fclose(prova_file);
    prova_file = fopen(nome_file, "a");
    if(rank==1){
        fprintf(prova_file,"I am process %d\n", rank);
        for(int i = 0; i < 2*5; i++)
            fprintf(prova_file, "%u ",ptr[i]);
    }
    fclose(prova_file);

    MPI_Finalize();





    
    return 0;
}