// To evolve ordered: evolve_ordered is a wrapper which decides whether to call OMP or MPI version.
//                   this is because MPI version does not work unless we have at least 2 processes (otherwise no messages to receive).
void evolve_ordered(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read, int s);
void evolve_ordered_OMP(unsigned char* current, int k, int n_steps, int s);
void evolve_ordered_MPI(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read, int s);
