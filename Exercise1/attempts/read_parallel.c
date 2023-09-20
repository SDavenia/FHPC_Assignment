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
    if(rank == 0){
        FILE* image_file; 
        image_file = fopen(image_name, "r");
        // *image = NULL; //address of the first element of image
        int* xsize;
        int* ysize;
        int* max_val;

        // *xsize = *ysize = *maxval = 0; // set to 0 the value of xsize, ysize and maxval
        *xsize = 0; 
        *ysize = 0;

        char    MagicN[2]; // define a string of 2 elements
        char   *line = NULL; //define a pointer "line" to NULL
        size_t  t, n = 0;
        int counter;
        // get the Magic Number
        t = fscanf(image_file, "%2s%*c", MagicN ); // This one reads P5
        counter += t;
        t = getline( &line, &n, image_file); // Here we read all the lines starting with #, i.e. all the comments.
        counter += t;
        while((line[0]=='#')){
            t = getline( &line, &n, image_file); // Here we read all the lines starting with #, i.e. all the comments.
            counter += t;
            printf("t is %zu\n", t);
        }
        

      
    }
    MPI_Finalize();
    
    return 0;
}