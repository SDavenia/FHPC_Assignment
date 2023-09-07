// FINAL VERSION OF THE OPTIMIZED MAIN (parallel)
// To compile: srun mpicc -fopenmp main_parallel.c -o main_parallel.exe
// To run executable to generate playground: srun ./main_parallel.exe -i -k 5 -f init.pgm
// To run execubtable to play on playground: srun ./main_parallel.exe -r -k 5 -f init.pgm -n 3

/*
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

void print_image(unsigned char* ptr, int k){
  /*
  This function is used only in intermediate steps to print to terminal a square matrix with k columns
  */
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            printf("%d ", ptr[i*k + j]/255);
        }
        printf("\n");
    }
}

// Here we report the functions which are used to generate the image, write to a pgm file, and read from a pgm file respectively
void initialize_parallel(int k, char *fname);
void write_pgm_parallel(unsigned char *ptr, int maxval, int xsize, int ysize, const char *fname, int rank, int size, int rows_initialize);
void read_pgm_parallel(unsigned char **ptr, int k, const char *image_name);
void read_pgm_parallel_frame(unsigned char **ptr, int k, const char *image_name);


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

  // 2- Depending on whether initialisation or execution is required, perform it.
  if(action == INIT){
    // create initial conditions
    //double Tstart_init = omp_get_wtime();
    initialize_parallel(k,fname);
    //double Time_init = omp_get_wtime() - Tstart_init;
    //printf("write time : %lf\n", Time_init);
  }else{ 
    // Read and run a playground
    unsigned char* input;
    read_pgm_parallel_frame(&input, k, fname);
    
    if(e == 0){ // Ordered
      printf("ORDERED EXECUTION\n");
      //evolve_dynamic_parallel(current, k, n);
    }else{ 
      printf("STATIC EXECUTION\n");    
    }
    
    
    if (fname != NULL)
      free(fname);
    free(input);
  }


  return 0;
}



void initialize_parallel(int k, char *fname){
  /*
  INPUT:
    - k: the size of the square where we want to execute GOL.
    - fname: the string containing the name of the pgm file where we want to store the image
  
  1. It generates a random matrix of size k x k. 
     Each MPI process is assigned a number of rows, and then OpenMP threads are used to randomly initialize different rows of the matrix
  
  2. It writes the random matrix to file called fname.
  */
  int rank, size;
  MPI_Init( NULL, NULL );
  MPI_Comm_rank( MPI_COMM_WORLD,&rank );
  MPI_Comm_size( MPI_COMM_WORLD,&size );  

  // Defines how many rows each process has to initialise
  int rows_initialize = (rank < k%size) ? (k / size + 1) : (k / size);
  
  MPI_File fh;
  MPI_Offset disp;

  // printf( "I am %d of %d and I have to generate %d rows\n", rank, size, rows_initialize);

  // Allocate the memory required for the rows of that processor.
  unsigned char* ptr = (unsigned char*)malloc(rows_initialize*k * sizeof(unsigned char)); // Allocate memory for the rows you have to generate.
  
  // In this parallel region the different threads generate random numbers on different sections of the matrix.
  //double Tstart_init = omp_get_wtime();
  #pragma omp parallel
  {
    int my_id = omp_get_thread_num();
    //int cpu_num = sched_getcpu(); // To see what core it is using
    unsigned int seed = clock();
    seed = seed * rank + my_id;
    // seed += my_id;
    // printf("I am thread %d of process %d and I am running on core %d\n", my_id, rank, cpu_num);

    #pragma omp for
    for (int i = 0; i < rows_initialize*k; i++){
        unsigned char random_num = (unsigned char) rand_r(&seed) % 2;
        //printf("random_num is %d\n", random_num);
        ptr[i] = (random_num==1) ? 255 : 0;
    }

  }
  // double Time_init = omp_get_wtime() - Tstart_init;
  // printf("I am process %d and generating random matrix took %lf s\n",rank, Time_init);

  // THIS SECTION HERE IS ONLY NEEDED FOR TESTING (HERE WE ONLY HAVE 3 PROCESSORS).
  /*
  char fname2[] = "prova_write.txt";
  // Open a FILE* stream
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 0){
      FILE* prova_file = fopen(fname2, "w");
      fprintf(prova_file,"I am process %d\n", rank);
      for(int i =0; i<k*rows_initialize;i++){
          fprintf(prova_file,"%u ", ptr[i]);
      }
      fclose(prova_file);
  }
  fprintf(prova_file,"\n");

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 1){
      FILE* prova_file = fopen(fname2, "a");
      fprintf(prova_file,"\nI am process %d\n", rank);
      for(int i =0; i<k*rows_initialize;i++){
          fprintf(prova_file,"%u ", ptr[i]);
      }
      fclose(prova_file);
  }
  fprintf(prova_file, "\n");
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 2){
      FILE* prova_file = fopen(fname2, "a");
      fprintf(prova_file,"\nI am process %d\n", rank);
      for(int i =0; i<k*rows_initialize;i++){
          fprintf(prova_file,"%u ", ptr[i]);
      }
      fclose(prova_file);
  }
  fprintf(prova_file, "\n");
  */
  

  // Write to file.
  write_pgm_parallel(ptr, 255, k, k, fname, rank, size, rows_initialize);

  free(ptr);
  MPI_Finalize();
}

void write_pgm_parallel( unsigned char *ptr, int maxval, int xsize, int ysize, const char *fname, int rank, int size, int rows_initialize){
  /*
  INPUT:
    - ptr: pointer to the memory location where the matrix is stored.
    - maxval: 255
    - xsize, ysize: size of the pgm image, in our case it is kxk
    - fname: name of the file where the pgm image is stored.
    - rank, size: MPI quantities
    - rows_initialize: how many rows each MPI process should write
  */
  
  MPI_File fh;
  MPI_Offset disp;

  MPI_File_delete(fname, MPI_INFO_NULL);
  MPI_File_open(  MPI_COMM_WORLD, fname, 
                MPI_MODE_CREATE | MPI_MODE_RDWR, 
                MPI_INFO_NULL, &fh  );
  MPI_File_close(&fh);

  // Write pgm header
  if (rank == 0) {
    FILE* file_stream = fopen(fname, "w");  // Open a FILE* stream
    /*int xsize, ysize, maxval;
    xsize = k;
    ysize = k;
    maxval = 255;*/

    if (file_stream != NULL) {
        int max_value = size - 1;
        fprintf(file_stream, "P5\n# generated by\n# Elena Rivaroli and Samuele D'Avenia\n%d %d\n%d\n", xsize, ysize, maxval);
        fclose(file_stream);
    } else {
        fprintf(stderr, "Failed to open the PGM file for writing.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD); // Wait for the header to have finished writing.
  // Open the file using MPI collective functions and write with each thread
  MPI_File_open(MPI_COMM_WORLD, fname, 
                MPI_MODE_APPEND | MPI_MODE_RDWR, 
                MPI_INFO_NULL, &fh  );
      
  // Decide where on file each MPI process should write.
  disp = (rank >= k % size) ? (rank * rows_initialize +k%size)* k *sizeof(unsigned char) : rank * rows_initialize * k *sizeof(unsigned char);
  
  MPI_File_seek(fh, disp, MPI_SEEK_CUR);
  MPI_File_write_all(fh, ptr, rows_initialize*k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
}

void read_pgm_parallel(unsigned char **ptr, int k, const char *image_name){
  /*
  INPUT:
    - ptr: pointer to the memory location where the matrix will be stored
    - k: matrix size
    - image_name: name of file where the matrix is stored.
  */
  int rank, size;
  MPI_Init( NULL, NULL );
  MPI_Comm_rank( MPI_COMM_WORLD,&rank );
  MPI_Comm_size( MPI_COMM_WORLD,&size );

  // Defines how many rows each process has to initialise
  int rows_read = k / size; 
  rows_read = (rank < k % size) ? rows_read+1 : rows_read;
  // printf("I am process %d and I have to read in my memory %d rows\n", rank, rows_read);

  *ptr = (unsigned char*)malloc(rows_read*k * sizeof(unsigned char));

  MPI_Offset disp;
  MPI_File   fh;
  MPI_File_open(  MPI_COMM_WORLD, image_name, 
                  MPI_MODE_RDONLY,
                  MPI_INFO_NULL, &fh  );

  
  disp = (rank >= k % size) ? (rank * rows_read + k % size) * k * sizeof(unsigned char) : rank * rows_read * k * sizeof(unsigned char);
  disp += 64;  // 64 in this case is the hardcoded representation of the header of the images
  // printf("I am process %d and I have to start writing at %d\n", rank, disp);

  MPI_File_seek(fh, disp, MPI_SEEK_CUR);
  MPI_File_read_all(fh, *ptr, rows_read*k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
  MPI_Barrier(MPI_COMM_WORLD);

  // This part is needed only for testing (ASSUMES 3 MPI processes are being used).
  
  /*
  FILE* prova_file;
  char nome_file[] = "prova_read.txt";
  if(rank==0){
    prova_file = fopen(nome_file, "w");
    fprintf(prova_file,"I am process %d\n", rank);
    for(int i = 0; i < rows_read * k; i++)
        fprintf(prova_file, "%u ",(*ptr)[i]);
    fprintf(prova_file, "\n");
    fclose(prova_file);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==1){
    prova_file = fopen(nome_file, "a");
    fprintf(prova_file,"I am process %d\n", rank);
    for(int i = 0; i < rows_read * k; i++)
        fprintf(prova_file, "%u ",(*ptr)[i]);
    fprintf(prova_file, "\n");
    fclose(prova_file);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==2){
    prova_file = fopen(nome_file, "a");
    fprintf(prova_file,"I am process %d\n", rank);
    for(int i = 0; i < rows_read * k; i++)
        fprintf(prova_file, "%u ",(*ptr)[i]);
    fprintf(prova_file, "\n");
    fclose(prova_file);
  }
  */
  MPI_Finalize();
  
}

void read_pgm_parallel_frame(unsigned char **ptr, int k, const char *image_name){
  /*
  INPUT:
    - ptr: pointer to the memory location where the matrix will be stored
    - k: matrix size
    - image_name: name of file where the matrix is stored.
  
  It initializes a memory area pointed by ptr with an additional row for the first and last process.
  To be more precise, the last row is also pasted before the first, and viceversa for the last. 
  This is needed
  */
  int rank, size;
  MPI_Init( NULL, NULL );
  MPI_Comm_rank( MPI_COMM_WORLD,&rank );
  MPI_Comm_size( MPI_COMM_WORLD,&size );

  // Defines how many rows each process has to initialise
  int rows_read = k / size; 
  rows_read = (rank < k % size) ? rows_read+1 : rows_read;
  // printf("I am process %d and I have to read in my memory %d rows\n", rank, rows_read);

  if (rank == 0 || rank == (size-1))
    *ptr = (unsigned char*)malloc((rows_read+1)*k * sizeof(unsigned char));
  else
    *ptr = (unsigned char*)malloc(rows_read*k * sizeof(unsigned char));

  MPI_Offset disp;
  MPI_File   fh;
  MPI_File_open(  MPI_COMM_WORLD, image_name, 
                  MPI_MODE_RDONLY,
                  MPI_INFO_NULL, &fh  );

  
  disp = (rank >= k % size) ? (rank * rows_read + k % size) * k * sizeof(unsigned char) : rank * rows_read * k * sizeof(unsigned char);
  disp += 64;  // 64 in this case is the hardcoded representation of the header of the images
  // printf("I am process %d and I have to start writing at %d\n", rank, disp);
  int u=0;
  if(rank == 0){
    MPI_Offset disp_temp;
    disp_temp = k*k + 64 - k;
    MPI_File_seek(fh, disp_temp, MPI_SEEK_SET);
    MPI_File_read(fh, *ptr, k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
    MPI_File_seek(fh, 0, MPI_SEEK_SET); // Needed because MPI_SEEK_CUR does offset + current, otherwise it goes after the end.
    u=k;
  }

  MPI_File_seek(fh, disp, MPI_SEEK_CUR);
  MPI_File_read_all(fh, (*ptr)+u, rows_read*k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

  if(rank == size -1){
    MPI_Offset disp_temp;
    disp_temp = 64;
    MPI_File_seek(fh, disp_temp, MPI_SEEK_SET);
    MPI_File_read(fh, (*ptr) + rows_read*k, k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
  }

  MPI_File_close(&fh);
  MPI_Barrier(MPI_COMM_WORLD);

  // This part is needed only for testing (ASSUMES 3 MPI processes are being used).
   /*
  FILE* prova_file;
  char nome_file[] = "prova_read.txt";
  if(rank==0){
    prova_file = fopen(nome_file, "w");
    fprintf(prova_file,"I am process %d\n", rank);
    for(int i = 0; i < (rows_read+1) * k; i++)
        fprintf(prova_file, "%u ",(*ptr)[i]);
    
    fprintf(prova_file, "\n");
    fclose(prova_file);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  if(rank==1){
      prova_file = fopen(nome_file, "a");
      fprintf(prova_file,"I am process %d\n", rank);
      for(int i = 0; i < rows_read * k; i++)
          fprintf(prova_file, "%u ",(*ptr)[i]);

     fprintf(prova_file, "\n");
     fclose(prova_file);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==1){
    prova_file = fopen(nome_file, "a");
    fprintf(prova_file,"I am process %d\n", rank);
    for(int i = 0; i < (rows_read+1) * k; i++)
        fprintf(prova_file, "%u ",(*ptr)[i]);
    fprintf(prova_file, "\n");
    fclose(prova_file);
  }
  */
  MPI_Finalize();

  
}