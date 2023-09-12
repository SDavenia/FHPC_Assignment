
// Here we report the functions which are used to generate the image, write to a pgm file, and read from a pgm file respectively
// To initialize and write the matrix
void initialize_parallel(int k, char *fname, int rank, int size, int rows_initialize);
void write_pgm_parallel(unsigned char *ptr, int maxval, int xsize, int ysize, const char *fname, int rank, int size, int rows_initialize);

// To read matrix from file
void read_pgm_parallel(unsigned char **ptr, int k, const char *image_name, int rank, int size, int rows_read);
