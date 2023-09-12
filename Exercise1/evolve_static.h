// To evolve static: evolve_static is a wrapper which decides whether to call OMP or MPI version.
//                   this is because MPI version does not work unless we have at least 2 processes (otherwise no messages to receive).
void evolve_static(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s);
void evolve_static_OMP(unsigned char* current, unsigned char* next, int k, int n_steps, int s);
void evolve_static_MPI(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s);
void evolve_static_MPI_blocking(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s);
