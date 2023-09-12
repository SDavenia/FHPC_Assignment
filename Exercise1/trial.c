#define _GNU_SOURCE // sched_getcpu(3) is glibc-specific (see the man page)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <sched.h> 

int main()
{
    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD,&rank );
    MPI_Comm_size( MPI_COMM_WORLD,&size );

    char name[50];
    int resultlen;
    MPI_Get_processor_name(name, &resultlen);
    
    #pragma omp parallel
    {
        int core = sched_getcpu();
        int mythread = omp_get_thread_num();
        printf("I am thread %d of process %d, I run on node %s and core %d\n", mythread, rank, name, core);
    }

    //printf("I am process %d, I run on node %s and core %d\n", rank, name, cpu);
    MPI_Finalize();

}

// Vogliamo vedere come si comporta con il map-by numa (forse quando ho più per farlo cosi map-by node --bind-to numa)
