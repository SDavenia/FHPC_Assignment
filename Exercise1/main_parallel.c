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
// To initialize and write the matrix
void initialize_parallel(int k, char *fname);
void write_pgm_parallel(unsigned char *ptr, int maxval, int xsize, int ysize, const char *fname, int rank, int size, int rows_initialize);

// To read matrix from file
void read_pgm_parallel_frame(unsigned char **ptr, int k, const char *image_name, int rank, int size, int rows_read);

// To evolve static: evolve_static is a wrapper which decides whether to call OMP or MPI version.
//                   this is because MPI version does not work unless we have at least 2 processes (otherwise no messages to receive).
void evolve_static(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s);
void evolve_static_OMP(unsigned char* current, unsigned char* next, int k, int n_steps, int s);
void evolve_static_MPI(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s);

// To evolve ordered: evolve_ordered is a wrapper which decides whether to call OMP or MPI version.
//                   this is because MPI version does not work unless we have at least 2 processes (otherwise no messages to receive).
void evolve_ordered(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read, int s);
void evolve_ordered_OMP(unsigned char* current, int k, int n_steps, int s);
void evolve_ordered_MPI(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read, int s);

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

  // 2- Depending on whether initialisation or execution is required, perform it.
  if(action == INIT){
    // create initial conditions
    //double Tstart_init = omp_get_wtime();
    initialize_parallel(k,file_path);
    //double Time_init = omp_get_wtime() - Tstart_init;
    //printf("write time : %lf\n", Time_init);
  }else{ 
    int rank, size;
    MPI_Init( NULL, NULL );
    MPI_Comm_rank( MPI_COMM_WORLD,&rank );
    MPI_Comm_size( MPI_COMM_WORLD,&size );

    // Defines how many rows each process has to initialise
    int rows_read = k / size; 
    rows_read = (rank < k % size) ? rows_read+1 : rows_read;
    
    // Read and run a playground
    unsigned char* current;
    read_pgm_parallel_frame(&current, k, file_path, rank, size, rows_read);

    unsigned char* next = (unsigned char*)malloc((rows_read+2)*k*sizeof(unsigned char));

    
    if(e == 0){ // Ordered
      printf("ORDERED EXECUTION\n");
      double Tstart_ord = omp_get_wtime();
      evolve_ordered(current, k, n, rank, size, rows_read, s);
      double Time_ord = omp_get_wtime() - Tstart_ord;
      printf("I am process %d and evolving the matrix statically for %d steps took %lf s\n",rank, n, Time_ord);
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

    MPI_Finalize();
  }
  free(file_path);


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

void read_pgm_parallel_frame(unsigned char **ptr, int k, const char *image_name, int rank, int size, int rows_read){
  /*
  INPUT:
    - ptr: pointer to the memory location where the matrix will be stored
    - k: matrix size
    - image_name: name of file where the matrix is stored.
  
  Each process gets allocated an amount of memory which corresponds to the rows it has to evolve plus the one above and below
    which are needed for correct update in parallel.
  */
  
  // printf("I am process %d and I have to read in my memory %d rows\n", rank, rows_read);

  // Added +2 since each process will also need the row above and below the rows it has to update.
  *ptr = (unsigned char*)malloc((rows_read+2)*k * sizeof(unsigned char));

  MPI_Offset disp;
  MPI_File   fh;
  MPI_File_open(  MPI_COMM_WORLD, image_name, 
                  MPI_MODE_RDONLY,
                  MPI_INFO_NULL, &fh  );

  // disp is the starting point for the rows each process has to evolve in the file
  disp = (rank >= k % size) ? (rank * rows_read + k % size) * k * sizeof(unsigned char) : rank * rows_read * k * sizeof(unsigned char);
  disp += 64;  // 64 in this case is the hardcoded representation of the header of the images
  // printf("I am process %d and I have to start writing at %d\n", rank, disp);
  
  // process 0 has to read the LAST row in the file
  //   all the others have to read the one before their "working" rows.
  if(rank==0)
    MPI_File_seek(fh, 64+k*k-k, MPI_SEEK_SET);
  else 
    MPI_File_seek(fh, disp-k, MPI_SEEK_SET);
  // Read into ptr the leftmost row
  MPI_File_read_all(fh, *ptr, k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

  // Now we read the working rows
  if(rank==0)
    MPI_File_seek(fh, 64, MPI_SEEK_SET);
  MPI_File_read_all(fh, (*ptr)+k, rows_read*k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

  // Finally we read the last row (symmetric case for the last processor as P0)
  if(rank==size-1)
    MPI_File_seek(fh, 64, MPI_SEEK_SET);
  MPI_File_read_all(fh, (*ptr)+k+rows_read*k, k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
  MPI_Barrier(MPI_COMM_WORLD);

  // This part is needed only for testing (ASSUMES 3 MPI processes are being used).
  

  /*
  FILE* prova_file;
  char nome_file[] = "prova_read.txt";
  if(rank==0){
    prova_file = fopen(nome_file, "w");
    fprintf(prova_file,"I am process %d\n", rank);
    for(int i = 0; i < (rows_read+2) * k; i++)
        fprintf(prova_file, "%u ",(*ptr)[i]);
    
    fprintf(prova_file, "\n");
    fclose(prova_file);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  if(rank==1){
      prova_file = fopen(nome_file, "a");
      fprintf(prova_file,"I am process %d\n", rank);
      for(int i = 0; i < (rows_read+2) * k; i++)
          fprintf(prova_file, "%u ",(*ptr)[i]);

     fprintf(prova_file, "\n");
     fclose(prova_file);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==2){
    prova_file = fopen(nome_file, "a");
    fprintf(prova_file,"I am process %d\n", rank);
    for(int i = 0; i < (rows_read+2) * k; i++)
        fprintf(prova_file, "%u ",(*ptr)[i]);
    fprintf(prova_file, "\n");
    fclose(prova_file);
  }
  */
    
}

void evolve_static(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s){
  if(size == 1)
    evolve_static_OMP(current, next, k, n_steps, s);
  else
    evolve_static_MPI(current, next,k , n_steps, rank, size, rows_read, s);
}

void evolve_static_MPI(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s){
  MPI_Request request[2]; // Array of MPI_Request, necessary to ensure that the non blocking receive is complete

  for(int n_step=0; n_step < n_steps; n_step++){
      int nthreads;
      #pragma omp parallel
      {
        int myid = omp_get_thread_num();
        #pragma omp master
          nthreads = omp_get_num_threads();
      
        #pragma omp for
        for(int i=1;i<rows_read+1;i++){
          //printf("I am thread %d doing row %d\n", myid, i);
          for(int j=0; j<k;j++){
              int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                  current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                  current[(i-1)*k + j] + current[(i+1)*k + j];
              next[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
          }
        }
        // Here there is an implicit barrier of the for loop
        /* QUESTI VANNO SOSTITUITI CON UN MESSAGGIO
        #pragma omp single nowait
        for (int i = 0; i < k; i++){
            next[i] = next[(k)*k + i];
        }

        #pragma omp single nowait
        for (int i = 0; i < k; i++){
            next[(k+1)*k+i]=next[k+i];
        }
        */
      }
      // Here there is an implcit barrier for the end of the parallel region.

      if(rank == 0){ // TAG IS ALWAYS THE RANK OF THE SENDING PROCESS + n_step
        // Send message to last process
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, size-1, rank + n_step + 1, MPI_COMM_WORLD, &request[0]);
        // Send message to the next process
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
        
        // Blocking receive message
        // Lower row receive
        MPI_Recv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Upper row receive from final process
        MPI_Recv(next, k, MPI_UNSIGNED_CHAR, size-1, size-1 + n_step + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      }else if(rank == size-1){
        // Send message to process before
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
        // Send message to process 0
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, rank + n_step +1 , MPI_COMM_WORLD, &request[1]);
        
        // BLOCKING RECEIVE
        // Lower row receive
        MPI_Recv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, 0+n_step+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Upper row receive
        MPI_Recv(next, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }else{
        // Send message to process before
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
        // Send message to process after
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
        
        // int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
        //   MPI_Comm comm, MPI_Status *status)
        // BLOCKING RECEIVE
        // Upper row receive
        MPI_Recv(next, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Lower row receive
        MPI_Recv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1+n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      MPI_Waitall(2, request, MPI_STATUS_IGNORE); // To ensure they all do not alter the buffer by moving to next iteration.

      unsigned char* tmp;
      tmp = next;
      next = current;
      current=tmp;

      /*
      // Questa parte serve solo a controllare
      if(n_step == 0){
        FILE* prova_file;
        char nome_file[] = "prova_read.txt";
        if(n_step==0){
          if(rank==0){
            prova_file = fopen(nome_file, "w");
            fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
            for(int i = 0; i < (rows_read+2) * k; i++)
                fprintf(prova_file, "%u ",current[i]);
            
            fprintf(prova_file, "\n");
            fclose(prova_file);
          }
        }else{
          if(rank==0){
            prova_file = fopen(nome_file, "w");
            fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
            for(int i = 0; i < (rows_read+2) * k; i++)
                fprintf(prova_file, "%u ",current[i]);
            
            fprintf(prova_file, "\n");
            fclose(prova_file);
          } 
        }

        MPI_Barrier(MPI_COMM_WORLD);
        
        if(rank==1){
            prova_file = fopen(nome_file, "a");
            fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
            for(int i = 0; i < (rows_read+2) * k; i++)
                fprintf(prova_file, "%u ",current[i]);

          fprintf(prova_file, "\n");
          fclose(prova_file);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==2){
          prova_file = fopen(nome_file, "a");
          fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
          for(int i = 0; i < (rows_read+2) * k; i++)
              fprintf(prova_file, "%u ",current[i]);
          fprintf(prova_file, "\n");
          fclose(prova_file);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==3){
          prova_file = fopen(nome_file, "a");
          fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
          for(int i = 0; i < (rows_read+2) * k; i++)
              fprintf(prova_file, "%u ",current[i]);
          fprintf(prova_file, "\n");
          fclose(prova_file);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==4){
          prova_file = fopen(nome_file, "a");
          fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
          for(int i = 0; i < (rows_read+2) * k; i++)
              fprintf(prova_file, "%u ",current[i]);
          fprintf(prova_file, "\n");
          fclose(prova_file);
        }
      }
      */

      //printf("Step %d:\n", n_step+1);
      //print_image(current, k+2,k);
    if(n_step % s == 0){
      char file_path[45] = "images/evolve_static/"; // Sufficiently large
      char filename[20];
      
      snprintf(filename, 20, "snapshot_%05d.pgm", n_step);
      strcat(file_path, filename);
      write_pgm_parallel(current+k, 255, k, k, file_path, rank, size, rows_read);
    }
  }

}

void evolve_static_OMP(unsigned char* current, unsigned char* next, int k, int n_steps, int s){
    for(int n_step=0; n_step < n_steps; n_step++){
      int nthreads;
      #pragma omp parallel
      {
        int myid = omp_get_thread_num();
        #pragma omp master
          nthreads = omp_get_num_threads();
      
        #pragma omp for
        for(int i=1;i<k+1;i++){
          //printf("I am thread %d doing row %d\n", myid, i);
          for(int j=0; j<k;j++){
              int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                  current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                  current[(i-1)*k + j] + current[(i+1)*k + j];
              next[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
          }
        }
        // Here there is an implicit barrier of the for loop

        #pragma omp single nowait
        for (int i = 0; i < k; i++){
            next[i] = next[(k)*k + i];
        }

        #pragma omp single nowait
        for (int i = 0; i < k; i++){
            next[(k+1)*k+i]=next[k+i];
        }
      }
      // Here there is an implicit barrier for the end of the parallel region.

      unsigned char* tmp;
      tmp = next;
      next = current;
      current=tmp;

      //printf("Step %d:\n", n_step+1);
      //print_image(current, k+2,k);
      if(n_step % s == 0){
      char file_path[45] = "images/evolve_static/"; // Sufficiently large
      char filename[20];
      
      snprintf(filename, 20, "snapshot_%05d.pgm", n_step);
      strcat(file_path, filename);
      write_pgm_parallel(current+k, 255, k, k, file_path, 0, 1, k+2);
      }
    
    }

}

void evolve_ordered(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read, int s){
  if(size == 1)
    evolve_ordered_OMP(current, k, n_steps, s);
  else
    evolve_ordered_MPI(current, k, n_steps, rank, size, rows_read, s);
}

void evolve_ordered_OMP(unsigned char* current, int k, int n_steps, int s){
  for(int n_step=0; n_step < n_steps; n_step++){
    /*
    FILE* prova_file;
    char nome_file[] = "prova_file.txt";
    prova_file=fopen(nome_file, "w");
    */
    int nthreads;
    #pragma omp parallel
    {
      int myid = omp_get_thread_num();
      #pragma omp master
        nthreads = omp_get_num_threads();
        
      #pragma omp for ordered
      for(int i=1;i<k;i++){
        for(int j=0; j<k;j++){
          #pragma omp ordered
          {
            // fprintf(prova_file,"I am thread %d and I am updating element (%d,%d)\n", myid,i,j);
            int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
            current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
            current[(i-1)*k + j] + current[(i+1)*k + j];
            current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
          }
        }
      }

      // Update bottom row of frame: since we are only copying no need for ordered
      #pragma omp for
      for(int j=0;j<k;j++){
        current[(k+1)*k+j] = current[k+j];
      }

      // Update last inner row of current and upper row of frame
      #pragma omp for ordered
      for(int j = 0; j< k; j++){
        #pragma omp ordered
        {
          int i = k;
          int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                current[(i-1)*k + j] + current[(i+1)*k + j];
          current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
          // Frame
          current[j] = current[i*k+j];
        }
      }
    }
    if(n_step % s == 0){
      char file_path[45] = "images/evolve_ordered/"; // Sufficiently large
      char filename[20];
      
      snprintf(filename, 20, "snapshot_%05d.pgm", n_step);
      strcat(file_path, filename);
      write_pgm_parallel(current+k, 255, k, k, file_path, 0, 1, k+2);
    }
    /*
    if(n_step == 50){
      printf("Step %d of OMP\n", n_step);
      print_image(current, k+2, k);
    }
    */
  }
}

void evolve_ordered_MPI(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read, int s){
  MPI_Request request[2];
  for(int n_step=0; n_step < n_steps; n_step++){

    // All processes wait until they have the current status of the matrix to update.

    if(rank == size -1){   // Lower row receive
      MPI_Recv(current + k + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, 0+n_step+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(rank != 0){ // Upper row receive
      MPI_Recv(current, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    

    // Evolve parallel using openMP
    int nthreads;
    #pragma omp parallel
    {
      int myid = omp_get_thread_num();
      #pragma omp master
        nthreads = omp_get_num_threads();
        
      #pragma omp for ordered
      for(int i=1;i<rows_read+1;i++){
        for(int j=0; j<k;j++){
          #pragma omp ordered
          {
            // fprintf(prova_file,"I am thread %d and I am updating element (%d,%d)\n", myid,i,j);
            int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
            current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
            current[(i-1)*k + j] + current[(i+1)*k + j];
            current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
          }
        }
      }
    } // Implicit barrier at the end of parallel region

    if(rank == 0){ // TAG IS ALWAYS THE RANK OF THE SENDING PROCESS + n_step
      // Send message to the next process
      MPI_Isend(current + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
      // Send message to last process
      MPI_Isend(current+k, k, MPI_UNSIGNED_CHAR, size-1, rank + n_step + 1, MPI_COMM_WORLD, &request[0]);
      
      // Blocking receive message
      // Lower row receive
      MPI_Recv(current + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // Upper row receive
      MPI_Recv(current, k, MPI_UNSIGNED_CHAR, size-1, size-1 + n_step + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }else if(rank == size-1){
      // Send message to process 0
      MPI_Isend(current + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, rank + n_step +1 , MPI_COMM_WORLD, &request[1]);
      // Send message to process before
      MPI_Isend(current+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
    }else{
      // Send message to process after
      MPI_Isend(current + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
      // Send message to process before
      MPI_Isend(current+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
      
      // BLOCKING RECEIVE
      // Lower row receive
      MPI_Recv(current + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1+n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    /*
    if(n_step == 50){
      FILE* prova_file;
      char nome_file[] = "prova_read3.txt";
      if(n_step==0){
        if(rank==0){
          prova_file = fopen(nome_file, "w");
          fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
          for(int i = 0; i < (rows_read+2) * k; i++)
              fprintf(prova_file, "%u ",current[i]);
          
          fprintf(prova_file, "\n");
          fclose(prova_file);
        }
      }else{
        if(rank==0){
          prova_file = fopen(nome_file, "w");
          fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
          for(int i = 0; i < (rows_read+2) * k; i++)
              fprintf(prova_file, "%u ",current[i]);
          
          fprintf(prova_file, "\n");
          fclose(prova_file);
        } 
      }

      MPI_Barrier(MPI_COMM_WORLD);
      
      if(rank==1){
          prova_file = fopen(nome_file, "a");
          fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
          for(int i = 0; i < (rows_read+2) * k; i++)
              fprintf(prova_file, "%u ",current[i]);

        fprintf(prova_file, "\n");
        fclose(prova_file);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if(rank==2){
        prova_file = fopen(nome_file, "a");
        fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
        for(int i = 0; i < (rows_read+2) * k; i++)
            fprintf(prova_file, "%u ",current[i]);
        fprintf(prova_file, "\n");
        fclose(prova_file);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if(rank==3){
        prova_file = fopen(nome_file, "a");
        fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
        for(int i = 0; i < (rows_read+2) * k; i++)
            fprintf(prova_file, "%u ",current[i]);
        fprintf(prova_file, "\n");
        fclose(prova_file);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if(rank==4){
        prova_file = fopen(nome_file, "a");
        fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
        for(int i = 0; i < (rows_read+2) * k; i++)
            fprintf(prova_file, "%u ",current[i]);
        fprintf(prova_file, "\n");
        fclose(prova_file);
      }
    }
    */
    
    // MPI_Waitall(2, request, MPI_STATUS_IGNORE);
    //MPI_Barrier(MPI_COMM_WORLD);
    if(n_step % s == 0){
      char file_path[45] = "images/evolve_ordered/"; // Sufficiently large
      char filename[20];
      
      snprintf(filename, 20, "snapshot_%05d.pgm", n_step);
      strcat(file_path, filename);
      write_pgm_parallel(current+k, 255, k, k, file_path, rank, size, rows_read);
    }
  }
}

