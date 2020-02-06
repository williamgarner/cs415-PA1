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

const char* dgemm_desc = "Simple blocked dgemm.";

#include <stdio.h>

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 41
#endif

#define min(a,b) (((a)<(b))?(a):(b))

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, double* A, double* B, double* C)
{
  /* For each row i of A */
  for (int i = 0; i < M; ++i)
    /* For each column j of B */
    for (int j = 0; j < N; ++j)
    {
      /* Compute C(i,j) */
      double cij = C[i+j*lda];
      for (int k = 0; k < K; ++k)
	cij += A[i+k*lda] * B[k+j*lda];
      C[i+j*lda] = cij;
    }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
//void square_dgemm_old (int lda, double* A, double* B, double* C)
//{
//  /* For each block-row of A */
//  for (int i = 0; i < lda; i += BLOCK_SIZE)
//    /* For each block-column of B */
//    for (int j = 0; j < lda; j += BLOCK_SIZE)
//      /* Accumulate block dgemms into block of C */
//      for (int k = 0; k < lda; k += BLOCK_SIZE)
//      {
//	/* Correct block dimensions if block "goes off edge of" the matrix */
//	int M = min (BLOCK_SIZE, lda-i);
//	int N = min (BLOCK_SIZE, lda-j);
//	int K = min (BLOCK_SIZE, lda-k);
//
//	/* Perform individual block dgemm */
//	do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
//      }
//}


void doBlock(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax, double* A, double* B, double* C)
{
//	printf("iMax: %d", iMax);
//	printf("jMax: %d", jMax);
//	printf("kMax: %d", kMax);
	for(int i = iMin; i < iMax; i++)
		for(int j = jMin; j < jMax; j++)
			for(int k = kMin; k < kMax; k++)
			{
//				printf("C:\n");
//				printf("\ti: %d\n", i);
//				printf("\tj: %d\n", j);
//				printf("\tk: %d\n", k);
//				printf("i + j * jMax: %d\n", i + j * jMax);
//				printf("i + k * kMax: %d\n", i + k * kMax);
//				printf("k + j * jMax: %d\n", k + j * jMax);


				C[i + j * jMax] += A[i + k * kMax] * B[k + j * jMax];

			}
}



#define LEVEL_1_BLOCK 60/sizeof(double)

void level1Block(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax, double* A, double* B, double* C, int lda)
{
	for(int iStart = iMin, iEnd = iMax; iStart < iMax; iStart = iEnd, iEnd = min(iEnd + LEVEL_1_BLOCK, iMax))
		for(int jStart = jMin, jEnd = jMax; jStart < jMax; jStart = jEnd, jEnd = min(jEnd + LEVEL_1_BLOCK, jMax))
			for(int kStart = kMin, kEnd = kMax; kStart < kMax; kStart = kEnd, kEnd = min(kEnd + LEVEL_1_BLOCK, kMax))
				doBlock(iStart, iEnd, jStart, jEnd, kStart, kEnd, A, B, C);
//				do_block(lda, iMax, jMax, kMax, A, B, C);
}

#define LEVEL_2_BLOCK 170/sizeof(double)

void level2Block(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax, double* A, double* B, double* C, int lda)
{
	for(int iStart = iMin, iEnd = iMax; iStart < iMax; iStart = iEnd, iEnd = min(iEnd + LEVEL_2_BLOCK, iMax))
		for(int jStart = jMin, jEnd = jMax; jStart < jMax; jStart = jEnd, jEnd = min(jEnd + LEVEL_2_BLOCK, jMax))
			for(int kStart = kMin, kEnd = kMax; kStart < kMax; kStart = kEnd, kEnd = min(kEnd + LEVEL_2_BLOCK, kMax))
				level1Block(iStart, iEnd, jStart, jEnd, kStart, kEnd, A, B, C, lda);
}

#define LEVEL_3_BLOCK 2019/sizeof(double)

void square_dgemm (int lda, double* A, double* B, double* C)
{
	for(int iStart = 0, iEnd = min(LEVEL_3_BLOCK, lda); iStart < lda; iStart = iEnd, iEnd = min(iEnd + LEVEL_3_BLOCK, lda))
		for(int jStart = 0, jEnd = min(LEVEL_3_BLOCK, lda); jStart < lda; jStart = jEnd, jEnd = min(jEnd + LEVEL_3_BLOCK, lda))
			for(int kStart = 0, kEnd = min(LEVEL_3_BLOCK, lda); kStart < lda; kStart = kEnd, kEnd = min(kEnd + LEVEL_3_BLOCK, lda))
				level2Block(iStart, iEnd, jStart, jEnd, kStart, kEnd, A, B, C, lda);

}
