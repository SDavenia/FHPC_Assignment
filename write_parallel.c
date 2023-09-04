/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Writing of a 1D buffer in parallel.

  Write a 1D buffer in parallel using 4 different versions.
  For a contiguous data type, use:

  VERSION == 1 employs shared file pointer
  VERSION == 2 employs individual file pointer with offset computed
               by the MPI_File_seek()
  VERSION == 3 defines a file view with offset depending on the process
               rank
  VERSION == 4 similar to VERSION == 3, but using a contiguous MPI
               datatype

  A non-contiguous version is handled with VERSION == 5, which defines
  a MPI vector datatype and a file view.
   
  \author A. Mignone (mignone@to.infn.it)
  \date   March 1, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define NELEM    3

#define VERSION  3
int main(int argc, char **argv)
{
  int i, rank, size;
  double buf[NELEM];
  char fname[] = "arr1D.bin";
  MPI_File   fh;
  MPI_Offset disp;

/* --------------------------------------------------------
   0. Initialize the MPI execution environment
   -------------------------------------------------------- */

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

/* --------------------------------------------------------
   1. Initialize array
   -------------------------------------------------------- */

  for (i = 0; i < NELEM; i++) buf[i] = rank;

/* --------------------------------------------------------
   2. Delete, re-open file and write
   -------------------------------------------------------- */

  MPI_File_delete(fname, MPI_INFO_NULL);
  MPI_File_open(  MPI_COMM_WORLD, fname, 
                  MPI_MODE_CREATE | MPI_MODE_RDWR, 
                  MPI_INFO_NULL, &fh  );

#if VERSION == 1     // Contiguous data, shared file pointer
  MPI_File_write_ordered(fh, buf, NELEM, MPI_DOUBLE, MPI_STATUS_IGNORE);  
#elif VERSION == 2   // Contiguous data, individual file pointer
  disp = rank*NELEM*sizeof(double);   // In bytes
  MPI_File_seek(fh, disp, MPI_SEEK_SET);
  MPI_File_write(fh, buf, NELEM, MPI_DOUBLE, MPI_STATUS_IGNORE);
#elif VERSION == 3   // Contiguous data, file view
  disp = rank*NELEM*sizeof(double);
  MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
  MPI_File_write(fh, buf, NELEM, MPI_DOUBLE, MPI_STATUS_IGNORE);  
#elif VERSION == 4   // Contiguous data, file view with MPI datatype
  MPI_Datatype cntg_type;

  MPI_Type_contiguous(NELEM, MPI_DOUBLE, &cntg_type);
  MPI_Type_commit(&cntg_type);

  disp = rank*NELEM*sizeof(double);

  MPI_File_set_view(fh, disp, MPI_BYTE, cntg_type, "native", MPI_INFO_NULL);
  MPI_File_write(fh, buf, 1, cntg_type, MPI_STATUS_IGNORE);  
  MPI_Type_free(&cntg_type);
#elif VERSION == 5   // Non-contiguous data, file view, vector type
  MPI_Datatype vec_type;

  for (i = 0; i < NELEM; i++) buf[i] = rank + 0.1*i;

  MPI_Type_vector(NELEM, 1, size, MPI_DOUBLE, &vec_type);
  MPI_Type_commit(&vec_type);

  disp = rank*sizeof(double);
  MPI_File_set_view(fh, disp, MPI_DOUBLE, vec_type, "native", MPI_INFO_NULL);
  MPI_File_write(fh, buf, NELEM, MPI_DOUBLE, MPI_STATUS_IGNORE);  
  MPI_Type_free(&vec_type);
#endif

  MPI_File_close(&fh);
  MPI_Finalize();
  return 0;
}
