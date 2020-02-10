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
const char* dgemm_desc = "Naive, three-loop dgemm.";

//#include <stdio.h>

#define min(a,b) (((a)<(b))?(a):(b))
#define BLOCK_SIZE 50

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */    
void square_dgemm (const int n, double*  A, double* B, double* restrict C)
{

//	printf("THIS IS THE SIZE: %d\n", n);
//			double T[n*n];
//			for(int blockStart = 0; blockStart < n; blockStart += BLOCK_SIZE)
//			{
//				int blockEnd = min(blockStart + BLOCK_SIZE, n);
//
//
////					int temp;
////					scanf("%d", &temp);
//				for(int i = blockStart; i < blockEnd; ++i)
//				{
//					for(int j = blockStart; j < blockEnd; ++j)
//					{
////							printf("%d, %d\n", i, j);
//
//						T[i*n + j] = A[j*n + i];
//					}
//				}
//			}
//
//
//
//				for (int j = 0; j < n; ++j)
//					for (int i = 0; i < n; ++i)
//						#pragma vector unaligned
//						for( int k = 0; k < n; ++k)
//								C[i+ j*n] += T[k+ i*n] * B[k+ j*n];


	double T[n*n];
	for(int i = 0; i < n; ++i)
#pragma vector unaligned
			for(int j = 0; j < n; ++j)
				T[i*n + j] = A[j*n + i];



	for( int k = 0; k < n; ++k)
		for (int i = 0; i < n; ++i)
#pragma vector unaligned
			for (int j = 0; j < n; ++j)



					C[i+ j*n] += T[k+ i*n] * B[k+ j*n];

}
