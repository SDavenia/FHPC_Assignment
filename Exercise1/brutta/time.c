#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9)

void main(){

    struct timespec ts;
    double Tstart = CPU_TIME ;
    double Time = CPU_TIME - Tstart;
    printf("%lf\n", Time);

}