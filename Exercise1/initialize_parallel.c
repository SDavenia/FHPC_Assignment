#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <sched.h> // Needed to find out which core each thread is running on

#define _GNU_SOURCE // sched_getcpu(3) is glibc-specific (see the man page)

// #define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9)
// To compile on mac /usr/local/bin/gcc-12 -fopenmp 

/*
 The matrix will be initialised in parallel, but for now for writing it to file that will
 be done in serial, since there are additional complications for doing that.

 N.B.
 If you use a shared seed this will create seed contention issues among the threads, leading to the code not scaling at all.
 In this case the seed is initialised to the time, and then to ensure it is different among threads, each one of them is added its thread it.
*/

int main(int argc, char* argv[])
{
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
    int rows_initialize = k / size; 
    if (rank < k%size) // For remainder 
        rows_initialize += 1;

    // printf( "I am %d of %d and I have to generate %d rows\n", rank, size, rows_initialize);

    unsigned char* ptr = (unsigned char*)calloc(rows_initialize*rows_initialize, sizeof(unsigned char)); // Allocate memory for the rows you have to generate.
    
    double Tstart_init = omp_get_wtime();
    #pragma omp parallel
    {
        int my_id = omp_get_thread_num();
        int cpu_num = sched_getcpu(); // To see what core it is using
        unsigned int seed = clock();
        seed += my_id;
        printf("I am thread %d of process %d and I am running on core %d\n", my_id, rank, cpu_num);

        #pragma omp for
        for (int i = 0; i < rows_initialize*rows_initialize; i++){
            unsigned char random_num = (unsigned char) rand_r(&seed) % 2;
            //printf("random_num is %d\n", random_num);
            ptr[i] = random_num==1 ? 255 : 0;
        }

    }
    
    double Time_init = omp_get_wtime() - Tstart_init;
    // printf("I am %d and generating random matrix took %lf s\n",rank, Time_init);
    
}

