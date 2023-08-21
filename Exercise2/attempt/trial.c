#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
    int size;
    if (argc > 1)
        size = atoi(argv[1]);
    else{
        printf("You need to pass the size\n");
        return 0;
    }
    
    printf("Size: %dx%dx%d\tTime: %d\tGFLOPS: %d\n", size, size, size, size/100, size*50/243);

    

}


/*
Using double

 Computing matrix product using gemm function via CBLAS interface

 Elapsed time 0.62182575 s


2000x2000x2000  0.062183 s      257.306810 GFLOPS
BLAS library

 This example computes real matrix C=alpha*A*B+beta*C using
 BLAS function dgemm, where A, B, and  C are matrices and
 alpha and beta are scalars

 Initializing data for matrix multiplication C=A*B for matrix
 A(2000x2000) and matrix B(2000x2000)

 Using double

 Computing matrix product using gemm function via CBLAS interface

 Elapsed time 0.68869930 s


2000x2000x2000  0.068870 s      232.322002 GFLOPS
*/