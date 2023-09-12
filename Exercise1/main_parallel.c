// FINAL VERSION OF THE OPTIMIZED MAIN (parallel)
// To compile: srun mpicc -fopenmp main_parallel.c -o main_parallel.exe
// To run executable to generate playground: mpirun -n 4 --map-by socket ./main_parallel.exe -i -k 5 -f init_00005.pgm
// To run execubtable to play on playground: mpirun -n 4 --map-by socket ./main_parallel.exe -r -k 5 -e 1 -f init_00005.pgm -n 5 -s 1
/*
  This code is an update on main_parallel_old.
  The difference is that read_pgm allocates the row before and after the ones that each process has to evolve.
  double Tstart_init = omp_get_wtime();
  double Time_init = omp_get_wtime() - Tstart_init;
  printf("I am process %d and generating random matrix took %lf s\n",rank, Time_init);
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "evolve_ordered.h"
#include "evolve_static.h"
#include "read_write_parallel.h"

#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1
#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9)


int   action = 0;
int   k      = K_DFLT;
int   e      = ORDERED;
int   n      = 10000;
int   s      = 1;
char *fname  = NULL;

int main ( int argc, char **argv )
{
  // 1- Read the command line arguments
  struct timespec ts;
  int action = 0;
  char *optstring = "irk:e:f:n:s:"; //optstring is a list of characters, each representing a single character option

  int c;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch(c) {
      //c takes the current option (ex. k)
      // optarg takes the value of the current option (ex. 400 for k)
      
    case 'i':
      action = INIT; break;
      
    case 'r':
      action = RUN; break;
      
    case 'k':
      k = atoi(optarg); break;

    case 'e':
      e = atoi(optarg); break;

    case 'f':
      fname = (char*)malloc( sizeof(optarg)+1 ); //fname now is an "array of char" variable (a string)
      sprintf(fname, "%s", optarg );
      break;

    case 'n':
      n = atoi(optarg); break;

    case 's':
      s = atoi(optarg); break;

    default :
      printf("argument -%c not known\n", c ); break;
    }
  }

  // Where the initial matrices are stored 
  char *file_path = (char*)malloc( sizeof(optarg)+ 1 +30);
  strcpy(file_path,"images/initial_matrices/");
  strcat(file_path,fname);

  int rank, size;
  MPI_Init( NULL, NULL );
  MPI_Comm_rank( MPI_COMM_WORLD,&rank );
  MPI_Comm_size( MPI_COMM_WORLD,&size );

  // Defines how many rows each process has to initialise
  int rows_read = k / size; 
  rows_read = (rank < k % size) ? rows_read+1 : rows_read;

  // 2- Depending on whether initialisation or execution is required, perform it.
  if(action == INIT){
    // create initial conditions
    //double Tstart_init = omp_get_wtime();
    initialize_parallel(k,file_path, rank, size, rows_read);
    //double Time_init = omp_get_wtime() - Tstart_init;
    //printf("write time : %lf\n", Time_init);
  }else{ 
    // Read and run a playground
    unsigned char* current;
    double Tstart_read = omp_get_wtime();
    //read_pgm_parallel(&current, k, file_path, rank, size, rows_read);
    read_pgm_parallel(&current, k, file_path, rank, size, rows_read);
    double Time_read = omp_get_wtime() - Tstart_read;
    printf("I am process %d, Read time : %lf\n", rank, Time_read);

    unsigned char* next = (unsigned char*)malloc((rows_read+2)*k*sizeof(unsigned char));
    /*
    printf("Initial one:\n");
    print_image(current, k+2, k);
    printf("\n");
    */
    
    if(e == 0){ // Ordered
      printf("ORDERED EXECUTION\n");
      double Tstart_ord = omp_get_wtime();
      evolve_ordered(current, k, n, rank, size, rows_read, s);
      double Time_ord = omp_get_wtime() - Tstart_ord;
      printf("I am process %d and evolving the matrix in ordered way for %d steps took %lf s\n",rank, n, Time_ord);
    }else{ 
      printf("STATIC EXECUTION\n");
      double Tstart_static = omp_get_wtime();
      evolve_static(current, next, k, n, rank, size, rows_read, s); 
      double Time_static = omp_get_wtime() - Tstart_static;
      printf("I am process %d and evolving the matrix statically for %d steps took %lf s\n",rank, n, Time_static);
    }
    
    
    if (fname != NULL)
      free(fname);
    free(next);
    free(current);
  }
  free(file_path);
  MPI_Finalize();


  return 0;
}

