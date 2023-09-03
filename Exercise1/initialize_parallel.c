#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9)
// To compile on mac /usr/local/bin/gcc-12 -fopenmp 

/*
 The matrix will be initialised in parallel, but for now for writing it to file that will
 be done in serial, since there are additional complications for doing that.

 N.B.
 If you use a shared seed this will create seed contention issues among the threads, leading to the code not scaling at all.
 In this case the seed is initialised to the time, and then to ensure it is different among threads, each one of them is added its thread it.
*/

void shared_random_playground(unsigned char* ptr, int k);

int main(int argc, char* argv[])
{
    int k;
    if (argc > 1)
        k = atoi(argv[1]);
    else
        k = 5;

    printf("Initialising matrix of size %d\n", k);
    unsigned char* ptr = (unsigned char*)calloc(k*k, sizeof(unsigned char)); // Allocate memory which will be used by all
    
    double Tstart_init = omp_get_wtime();
    shared_random_playground(ptr, k);
    double Time_init = omp_get_wtime() - Tstart_init;
    printf("Generating random matrix took %lf s\n", Time_init);
    
}

void shared_random_playground(unsigned char* ptr, int k){
    struct timespec ts;
    
    #pragma omp parallel
    {
        int my_id = omp_get_thread_num();
        unsigned int seed = clock();
        seed += my_id;

        #pragma omp for
        for (int i = 0; i < k*k; i++){
            unsigned char random_num = (unsigned char) rand_r(&seed) % 2;
            //printf("random_num is %d\n", random_num);
            ptr[i] = random_num==1 ? 255 : 0;
        }  

    }
}

