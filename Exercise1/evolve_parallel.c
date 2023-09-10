#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
/*
  double Tstart_init = omp_get_wtime();
  double Time_init = omp_get_wtime() - Tstart_init;
  printf("I am process %d and generating random matrix took %lf s\n",rank, Time_init);
*/

/*
THIS IS THE CORRECT VERSION FOR THE SERIAL CODE, with frame only above and below.
Use this as reference to check whether the Conway evolution was correct.
*/

void read_pgm_parallel_frame(unsigned char **ptr, int k, const char *image_name, int rank, int size, int rows_read);
void initialize_current(unsigned char* input, unsigned char* current, int k);
void evolve_static_OMP(unsigned char* current, unsigned char* next, int k, int n_steps);
void evolve_static_MPI(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read);
void evolve_ordered_OMP(unsigned char* current, int k, int n_steps);
void evolve_ordered_MPI(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read);

void print_image(unsigned char* ptr, int nrow, int ncol){
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            printf("%d ", ptr[i*ncol + j]/255);
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]){
    printf("Run\n");
    int xsize;
    int ysize;
    int maxval;
    char fname[] = "init.pgm";
    int k;
    if(argc > 1)
      k = atoi(argv[1]);
    else
      k = 5;
    
    unsigned char* current;

    // MPI STUFF
    int rank, size;
    MPI_Init( NULL, NULL );
    MPI_Comm_rank( MPI_COMM_WORLD,&rank );
    MPI_Comm_size( MPI_COMM_WORLD,&size );
  
    int rows_read = k / size; 
    rows_read = (rank < k % size) ? rows_read+1 : rows_read;

    read_pgm_parallel_frame(&current, k, fname, rank, size, rows_read);
    
    // unsigned char* next = (unsigned char*)malloc((rows_read+2)*k*sizeof(unsigned char));
    // evolve_static_MPI(current, next, k, 100, rank, size, rows_read);
    /*
    double Tstart_stat = omp_get_wtime();
    double Time_stat = omp_get_wtime() - Tstart_stat;
    printf("Time %lf\n", Time_stat);
    */
    //double Tstart_ord = omp_get_wtime();
    //evolve_ordered_OMP(current, k, 51);
    evolve_ordered_MPI(current, k, 51, rank, size, rows_read);
    //double Time_ord = omp_get_wtime() - Tstart_ord;
    //printf("Time %lf\n", Time_ord);

    free(current);
    //free(next);
    MPI_Finalize();

    return 0;
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
  /*
  int rank, size;
  MPI_Init( NULL, NULL );
  MPI_Comm_rank( MPI_COMM_WORLD,&rank );
  MPI_Comm_size( MPI_COMM_WORLD,&size );
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

void evolve_static_OMP(unsigned char* current, unsigned char* next, int k, int n_steps){
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
    }

}

void evolve_static_MPI(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read){
  
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
      }
      // Here there is an implcit barrier for the end of the parallel region.

      if(rank == 0){ // TAG IS ALWAYS THE RANK OF THE SENDING PROCESS + n_step
        // Send message to last process
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, size-1, rank + n_step + 1, MPI_COMM_WORLD, &request[0]);
        // Send message to the next process
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
        
        // Blocking receive message
        // Upper row receive
        MPI_Recv(next, k, MPI_UNSIGNED_CHAR, size-1, size-1 + n_step + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Lower row receive
        MPI_Recv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      }else if(rank == size-1){
        // Send message to process before
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
        // Send message to process 0
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, rank + n_step +1 , MPI_COMM_WORLD, &request[1]);
        
        // int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
        //   MPI_Comm comm, MPI_Status *status)
        // BLOCKING RECEIVE
        // Upper row receive
        MPI_Recv(next, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Lower row receive
        MPI_Recv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, 0+n_step+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
      MPI_Waitall(2, request, MPI_STATUS_IGNORE);

      unsigned char* tmp;
      tmp = next;
      next = current;
      current=tmp;

      // Questa parte serve solo a controllare
      if(n_step == 50){
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

      //printf("Step %d:\n", n_step+1);
      //print_image(current, k+2,k);
  }

}

void evolve_ordered_OMP(unsigned char* current, int k, int n_steps){
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
    if(n_step == 50){
      printf("Step %d of OMP\n", n_step);
      print_image(current, k+2, k);
    }
  }
}

void evolve_ordered_MPI(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read){
  MPI_Request request[2];
  for(int n_step=0; n_step < n_steps; n_step++){

    // All processes wait until they have the current status of the matrix to update.
    if(rank != 0){ // Upper row receive
      MPI_Recv(current, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
      
    if(rank == size -1){   // Lower row receive
      MPI_Recv(current + k + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, 0+n_step+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
    }
    if(rank == 0){ // TAG IS ALWAYS THE RANK OF THE SENDING PROCESS + n_step
      // Send message to last process
      MPI_Isend(current+k, k, MPI_UNSIGNED_CHAR, size-1, rank + n_step + 1, MPI_COMM_WORLD, &request[0]);
      // Send message to the next process
      MPI_Isend(current + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
      
      // Blocking receive message
      // Upper row receive
      MPI_Recv(current, k, MPI_UNSIGNED_CHAR, size-1, size-1 + n_step + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // Lower row receive
      MPI_Recv(current + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }else if(rank == size-1){
      // Send message to process before
      MPI_Isend(current+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
      // Send message to process 0
      MPI_Isend(current + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, rank + n_step +1 , MPI_COMM_WORLD, &request[1]);
    }else{
      // Send message to process before
      MPI_Isend(current+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
      // Send message to process after
      MPI_Isend(current + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
      
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
    MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Waitall(2, request, MPI_STATUS_IGNORE);
  }
}

