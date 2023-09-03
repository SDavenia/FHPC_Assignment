// FINAL VERSION OF THE OPTIMIZED MAIN (ONLY SERIAL FUNCTIONS)
// To compile: gcc main.c -o main.exe
// To run executable to generate playground: ./main.exe -i -k 5 -f init.pgm
// To run execubtable to play on playground: ./main.exe -r -k 5 -f init.pgm -n 3
// OPTIMIZED and with the errors fixed
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
//#include <mpi.h>
#include <omp.h>
#ifndef _OPENMP
#error "openmp support is required to compile this code"
#endif


void initialize_current(unsigned char* input, unsigned char* current, int k);

#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1
#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9)


char fname_deflt[] = "game_of_life.pgm";

int   action = 0;
int   k      = K_DFLT;
int   e      = ORDERED;
int   n      = 10000;
int   s      = 1;
char *fname  = NULL;

void print_image(unsigned char* ptr, int ncol){
    for(int i = 0; i < ncol; i++){
        for(int j = 0; j < ncol; j++){
            printf("%d ", ptr[i*ncol + j]/255);
        }
        printf("\n");
    }
}

void random_playground(int k, char *fname){
   
  unsigned char* ptr = (unsigned char*)calloc(k*k, sizeof(unsigned char)); // creates a k*k array of unsigned char
  
  #pragma omp parallel
  {
  int me = omp_get_thread_num();
  //printf("Thread (pid: %d, tid: %ld) nr %d:\n", syscall(SYS_gettid), me); // process ID, Thread ID, number of thread inside the process
  printf("Thread nr %d\n", me);

  // generate a random matrix of 0 and 255
  unsigned int seed = clock();

  #pragma omp for
  for (int i = 0; i < k*k; i++) {
    unsigned char rand_num = (unsigned char) rand_r(&seed) % 2;
    ptr[i] = rand_num==1 ? 255 : 0;
  }
  //write_pgm_image(ptr, 255, k, k,fname);
  }
  print_image(ptr, k);
  free(ptr);
}

int main ( int argc, char **argv )
{
  struct timespec ts;
  int action = 0;
  char *optstring = "irk:e:f:n:s:"; //optstring is a list of characters, each representing a single character option

  int c;
  /*
  Generally, the getopt() function is called from inside of a loop’s conditional statement.
  The loop terminates when the getopt() function returns -1.
  A switch statement is then executed with the value returned by getopt() function.
  */
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
      printf("%s \n",fname);
      break;

    case 'n':
      n = atoi(optarg); break;

    case 's':
      s = atoi(optarg); break;

    default :
      printf("argument -%c not known\n", c ); break;
    }
  }
  random_playground(k,fname);

  if(action == INIT){
    // create initial conditions
    printf("Initialize\n");
    random_playground(k,fname);
  }else{ 
    // Read and run a playground
    printf("Boh ci pensiamo dopo\n");
  }

  return 0;
}

void initialize_current(unsigned char* input, unsigned char* current, int k){

  current[0] = input[k*k-1]; // Top left of current
  // Last row of input -> First row of current
  for (int i = 0; i < k; i++){
    current[i+1] = input[(k-1)*k + i];
  }
  current[(k+2) - 1] = input[(k-1)*k]; // Top right of current


  // Initialize the inner values
  for(int i = 0; i < k; i++){
    // First column of that row
    current[(i+1)*(k+2)] = input[(i+1)*k - 1];
    for(int j = 0; j < k; j++){
      current[(i+1)*(k+2) + (j+1)] = input[i*k + j];
    }
    // Last column of that row
    current[(i+2)*(k+2)-1] = input[i*k];
  }

  current[(k+1)*(k+2)] = input[k-1]; // Bottom left of current

  // Initialize the corners and the frame rows and columns
  // First row of input -> Last row of current
  for (int i = 0; i < k; i++){
    current[(k+1)*(k+2) + i + 1]=input[i];
  }
  current[(k+2)*(k+2)-1] = input[0]; // Bottom right of current
}
