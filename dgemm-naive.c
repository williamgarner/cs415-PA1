/* 
    Please include compiler name below (you may also include any other modules you would like to be loaded)

COMPILER= gnu

    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
 
CC = cc
OPT = -O3
CFLAGS = -Wall -std=gnu99 $(OPT)
MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

*/

#include <string.h>

#define C(row,col) C[row*n + col]
#define A(row, col) A[row*n + col]
#define B(row, bCol) B[row*n + bCol]
#define likely(x) __builtin_expect(!!(x), 1)

const char* dgemm_desc = "Naive, three-loop dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */    
void square_dgemm (int n, double* A, double* B, double* C)
{
//	memset(C, 0, sizeof(C));
//	for(int row = 0; row < n; ++row)
//	{
//		for(int col = 0; col < n; ++col)
//		{
//			for(int bCol = 0; bCol < n; ++bCol)
//			{
//				C(row, col) += A(row, col) * B(row, bCol);
//			}
//		}
//	}


	for (int i = 0; i < n; i++)
		for (int k = 0; k < n; k++)
			for (int j = 0; j < n; j++)
				if (likely(k))
					C[i*n + j] += A[i*n + k] * B[k*n + j];
				else
					C[i*n + j] = A[i*n + k] * B[k*n + j];


  /* For each row i of A */
//  for (int i = 0; i < n; ++i)
//    /* For each column j of B */
//    for (int j = 0; j < n; ++j)
//    {
//      /* Compute C(i,j) */
//      double cij = C[i+j*n];
//      for( int k = 0; k < n; k++ )
//      	cij += A[i+k*n] * B[k+j*n];
//      C[i+j*n] = cij;
//    }
}
