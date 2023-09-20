#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>


void print_image(unsigned char* ptr, int nrow, int ncol){
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            printf("%d ", ptr[i*ncol + j]/255);
        }
        printf("\n");
    }
}

void read_pgm_parallel_messages_Recv(unsigned char **ptr, int k, const char *image_name, int rank, int size, int rows_read){
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
  disp += 72;  // 64 in this case is the hardcoded representation of the header of the images
  // printf("I am process %d and I have to start writing at %d\n", rank, disp);
  
  // process 0 has to read the LAST row in the file
  //   all the others have to read the one before their "working" rows.
  /*if(rank==0)
    MPI_File_seek(fh, 64+k*k-k, MPI_SEEK_SET);
  else 
    MPI_File_seek(fh, disp-k, MPI_SEEK_SET);
    
  // Read into ptr the leftmost row
  MPI_File_read_all(fh, *ptr, k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
  */
 
  // Now we read the working rows
  /*if(rank==0)
    MPI_File_seek(fh, 64, MPI_SEEK_SET);
    */
  MPI_File_seek(fh, disp, MPI_SEEK_SET);
  MPI_File_read_all(fh, (*ptr)+k, rows_read*k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Request request[2];


  if(rank == 0){ // TAG IS ALWAYS THE RANK OF THE SENDING PROCESS + 10
    // Send message to last process
    MPI_Isend((*ptr)+k, k, MPI_UNSIGNED_CHAR, size-1, rank + 10 + 1, MPI_COMM_WORLD, &request[0]);
    // Send message to the next process
    MPI_Isend((*ptr) + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + 10, MPI_COMM_WORLD, &request[1]);
    
    // Blocking receive message
    // Lower row receive
    MPI_Recv((*ptr) + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1 + 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Upper row receive from final process
    MPI_Recv((*ptr), k, MPI_UNSIGNED_CHAR, size-1, size-1 + 10 + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  }else if(rank == size-1){
    // Send message to process before
    MPI_Isend((*ptr)+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + 10, MPI_COMM_WORLD, &request[0]);
    // Send message to process 0
    MPI_Isend((*ptr) + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, rank + 10 +1 , MPI_COMM_WORLD, &request[1]);
    
    // BLOCKING RECEIVE
    // Lower row receive
    MPI_Recv((*ptr) + k + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, 0+10+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Upper row receive
    MPI_Recv((*ptr), k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }else{
    // Send message to process before
    MPI_Isend((*ptr)+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + 10, MPI_COMM_WORLD, &request[0]);
    // Send message to process after
    MPI_Isend((*ptr) + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + 10, MPI_COMM_WORLD, &request[1]);
    
    // int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
    //   MPI_Comm comm, MPI_Status *status)
    // BLOCKING RECEIVE
    // Upper row receive
    MPI_Recv((*ptr), k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Lower row receive
    MPI_Recv((*ptr) + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1+10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  MPI_Waitall(2, request, MPI_STATUS_IGNORE); // To ensure they all do not alter the buffer by moving to next iteration.




/*
  // Finally we read the last row (symmetric case for the last processor as P0)
  if(rank==size-1)
    MPI_File_seek(fh, 64, MPI_SEEK_SET);
  MPI_File_read_all(fh, (*ptr)+k+rows_read*k, k, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
*/
  

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
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==3){
    prova_file = fopen(nome_file, "a");
    fprintf(prova_file,"I am process %d\n", rank);
    for(int i = 0; i < (rows_read+2) * k; i++)
      fprintf(prova_file, "%u ",(*ptr)[i]);
    fprintf(prova_file, "\n");
    fclose(prova_file);
  }
  */
  
    
}

