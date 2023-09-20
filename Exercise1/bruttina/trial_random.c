#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>


int main(int argc, char* argv[]){

    int k;
    if(argc > 1)
        k = atoi(argv[1]);
    else
        k = 5;


    int rank, size;
    MPI_Init( NULL, NULL );
    MPI_Comm_rank( MPI_COMM_WORLD,&rank );
    MPI_Comm_size( MPI_COMM_WORLD,&size );

    int rows_read = k / size; 
    rows_read = (rank < k % size) ? rows_read+1 : rows_read;
    printf("I am and I have to do %d rows\n", rows_read);
    unsigned char* ptr = (unsigned char*)malloc(rows_read*k * sizeof(unsigned char));

    double Tstart_generate;
    //if(rank == 0) 
    Tstart_generate = omp_get_wtime();
    #pragma omp parallel 
    {
        int s = 10000;
        

    }
    //MPI_Barrier(MPI_COMM_WORLD);
    //if(rank == 0){
        double Time_generate = omp_get_wtime() - Tstart_generate;
        printf("I am %d and Generate time: %lf\n",rank, Time_generate);
    //}
    MPI_Finalize();

}