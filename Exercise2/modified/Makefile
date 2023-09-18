### MKL libraries
###
###
MKL= -L${MKLROOT}/lib/intel64  -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

### OpenBLAS libraries 
OPENBLASROOT=${OPENBLAS_ROOT}

### BLIS library
BLISROOT=/u/dssc/sdaven00/myblis

#cpu: gemm_mkl_double.x gemm_oblas_double.x gemm_mkl_float.x gemm_oblas_float.x gemm_blis_double.x gemm_blis_float.x
cpu: gemm_mkl_double.x #gemm_oblas_double.x gemm_blis_double.x


gemm_mkl_double.x: gemm_modified.c
	gcc -fopenmp -DUSE_DOUBLE -DMKL $^ -m64 -I${MKLROOT}/include $(MKL)  -o $@

#gemm_oblas_double.x: gemm_modified.c
#	gcc -DUSE_DOUBLE -DOPENBLAS $^ -m64 -I${OPENBLASROOT}/include -L/${OPENBLASROOT}/lib -lopenblas -lpthread -o $@ -fopenmp

#gemm_blis_double.x: gemm_modified.c
#	gcc -fopenmp -DUSE_DOUBLE  -DBLIS $^ -m64 -I${BLISROOT}/include/blis -L/${BLISROOT}/lib -o $@ -lpthread  -lblis -fopenmp -lm

clean:
	rm -rf *.x
	rm -rf *.txt
	rm -rf *.out
