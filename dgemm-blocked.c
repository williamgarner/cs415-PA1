///*
//    Please include compiler name below (you may also include any other modules you would like to be loaded)
//
//COMPILER= gnu
//
//    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
//
//CC = cc
//OPT = -O3
//CFLAGS = -Wall -std=gnu99 $(OPT)
//MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
//LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
//
//*/
//
//const char* dgemm_desc = "Simple blocked dgemm.";
//
//#include <stdio.h>
//
//#if !defined(BLOCK_SIZE)
//#define BLOCK_SIZE 41
//#endif
//
//#define min(a,b) (((a)<(b))?(a):(b))
//
///* This auxiliary subroutine performs a smaller dgemm operation
// *  C := C + A * B
// * where C is M-by-N, A is M-by-K, and B is K-by-N. */
//static void do_block (int lda, int M, int N, int K, double* A, double* B, double* C)
//{
//  /* For each row i of A */
//  for (int i = 0; i < M; ++i)
//    /* For each column j of B */
//    for (int j = 0; j < N; ++j)
//    {
//      /* Compute C(i,j) */
//      double cij = C[i+j*lda];
//      for (int k = 0; k < K; ++k)
//	cij += A[i+k*lda] * B[k+j*lda];
//      C[i+j*lda] = cij;
//    }
//}
//
///* This routine performs a dgemm operation
// *  C := C + A * B
// * where A, B, and C are lda-by-lda matrices stored in column-major format.
// * On exit, A and B maintain their input values. */
////void square_dgemm_old (int lda, double* A, double* B, double* C)
////{
////  /* For each block-row of A */
////  for (int i = 0; i < lda; i += BLOCK_SIZE)
////    /* For each block-column of B */
////    for (int j = 0; j < lda; j += BLOCK_SIZE)
////      /* Accumulate block dgemms into block of C */
////      for (int k = 0; k < lda; k += BLOCK_SIZE)
////      {
////	/* Correct block dimensions if block "goes off edge of" the matrix */
////	int M = min (BLOCK_SIZE, lda-i);
////	int N = min (BLOCK_SIZE, lda-j);
////	int K = min (BLOCK_SIZE, lda-k);
////
////	/* Perform individual block dgemm */
////	do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
////      }
////}
//
//
//void doBlock(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax, double* A, double* B, double* C)
//{
////	printf("iMax: %d", iMax);
////	printf("jMax: %d", jMax);
////	printf("kMax: %d", kMax);
//	for(int i = iMin; i < iMax; i++)
//		for(int j = jMin; j < jMax; j++)
//			for(int k = kMin; k < kMax; k++)
//			{
////				printf("C:\n");
////				printf("\ti: %d\n", i);
////				printf("\tj: %d\n", j);
////				printf("\tk: %d\n", k);
////				printf("i + j * jMax: %d\n", i + j * jMax);
////				printf("i + k * kMax: %d\n", i + k * kMax);
////				printf("k + j * jMax: %d\n", k + j * jMax);
//
//
////				C[i + j * jMax] += A[i + k * kMax] * B[k + j * jMax];
//
//			}
//}
//
//
//
//#define LEVEL_1_BLOCK 60/sizeof(double)
//
//void level1Block(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax, double* A, double* B, double* C, int lda)
//{
//	for(int iStart = iMin, iEnd = iMax; iStart < iMax; iStart = iEnd, iEnd = min(iEnd + LEVEL_1_BLOCK, iMax))
//		for(int jStart = jMin, jEnd = jMax; jStart < jMax; jStart = jEnd, jEnd = min(jEnd + LEVEL_1_BLOCK, jMax))
//			for(int kStart = kMin, kEnd = kMax; kStart < kMax; kStart = kEnd, kEnd = min(kEnd + LEVEL_1_BLOCK, kMax))
////				doBlock(iStart, iEnd, jStart, jEnd, kStart, kEnd, A, B, C, lda);
//				do_block(lda, iMax, jMax, kMax, A, B, C);
//}
//
//#define LEVEL_2_BLOCK 170/sizeof(double)
//
//void level2Block(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax, double* A, double* B, double* C, int lda)
//{
//	for(int iStart = iMin, iEnd = iMax; iStart < iMax; iStart = iEnd, iEnd = min(iEnd + LEVEL_2_BLOCK, iMax))
//		for(int jStart = jMin, jEnd = jMax; jStart < jMax; jStart = jEnd, jEnd = min(jEnd + LEVEL_2_BLOCK, jMax))
//			for(int kStart = kMin, kEnd = kMax; kStart < kMax; kStart = kEnd, kEnd = min(kEnd + LEVEL_2_BLOCK, kMax))
//				level1Block(iStart, iEnd, jStart, jEnd, kStart, kEnd, A, B, C, lda);
//}
//
//#define LEVEL_3_BLOCK 2019/sizeof(double)
//
//void square_dgemm (int lda, double* A, double* B, double* C)
//{
//	for(int iStart = 0, iEnd = min(LEVEL_3_BLOCK, lda); iStart < lda; iStart = iEnd, iEnd = min(iEnd + LEVEL_3_BLOCK, lda))
//		for(int jStart = 0, jEnd = min(LEVEL_3_BLOCK, lda); jStart < lda; jStart = jEnd, jEnd = min(jEnd + LEVEL_3_BLOCK, lda))
//			for(int kStart = 0, kEnd = min(LEVEL_3_BLOCK, lda); kStart < lda; kStart = kEnd, kEnd = min(kEnd + LEVEL_3_BLOCK, lda))
//				level2Block(iStart, iEnd, jStart, jEnd, kStart, kEnd, A, B, C, lda);
//
////	for(int i = 0; i < lda * lda; i++)
////		for(int j = 0; j < lda; ++j)
////			for(int k = 0; k < lda; ++k)
////			{
////				C[i] = -1;
//////				printf("%lf ", C[i]);
////			}
//
//}
//



#include <stdio.h>
#include <mmintrin.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <emmintrin.h>

const char* dgemm_desc = "our optimized dgemm";

#if !defined(BLOCK_L1)
#define BLOCK_L1 256
#endif

#if !defined(BLOCK_L2)
#define BLOCK_L2 512
#endif

#define A(i,j) A[(j)*lda + (i)]
#define B(i,j) B[(j)*lda + (i)]
#define C(i,j) C[(j)*lda + (i)]
#define min(a,b) (((a)<(b))?(a):(b))

typedef __m128d sse_reg;

/* Perform SSE arithmetic to fill a 4x4 block of matrix C */
static void do_4x4 (int lda, int K, double* a, double* b, double* c)
{
	sse_reg
			a0x_1x, a2x_3x,
//			bx0, bx1, bx2, bx3,
			c00_10, c20_30,
			c01_11, c21_31,
			c02_12, c22_32,
			c03_13, c23_33;


	double* c01_11_pntr = c + lda;
	double* c02_12_pntr = c01_11_pntr + lda;
	double* c03_13_pntr = c02_12_pntr + lda;

	c00_10 = _mm_loadu_pd(c);
	c20_30 = _mm_loadu_pd(c+2);
	c01_11 = _mm_loadu_pd(c01_11_pntr);
	c21_31 = _mm_loadu_pd(c01_11_pntr+2);
	c02_12 = _mm_loadu_pd(c02_12_pntr);
	c22_32 = _mm_loadu_pd(c02_12_pntr+2);
	c03_13 = _mm_loadu_pd(c03_13_pntr);
	c23_33 = _mm_loadu_pd(c03_13_pntr+2);

	for (int x = 0; x < K; ++x) {
		a0x_1x = _mm_load_pd(a);
		a2x_3x = _mm_load_pd(a+2);
		a += 4;

		__m128d bx0 = _mm_loaddup_pd(b++);
		__m128d bx1 = _mm_loaddup_pd(b++);
		__m128d bx2 = _mm_loaddup_pd(b++);
		__m128d bx3 = _mm_loaddup_pd(b++);

		c00_10 = _mm_add_pd(c00_10, _mm_mul_pd(a0x_1x, bx0));
		c20_30 = _mm_add_pd(c20_30, _mm_mul_pd(a2x_3x, bx0));
		c01_11 = _mm_add_pd(c01_11, _mm_mul_pd(a0x_1x, bx1));
		c21_31 = _mm_add_pd(c21_31, _mm_mul_pd(a2x_3x, bx1));
		c02_12 = _mm_add_pd(c02_12, _mm_mul_pd(a0x_1x, bx2));
		c22_32 = _mm_add_pd(c22_32, _mm_mul_pd(a2x_3x, bx2));
		c03_13 = _mm_add_pd(c03_13, _mm_mul_pd(a0x_1x, bx3));
		c23_33 = _mm_add_pd(c23_33, _mm_mul_pd(a2x_3x, bx3));
	}

	/* Fill 4x4 block of C with results */
	_mm_storeu_pd(c, c00_10);
	_mm_storeu_pd((c+2), c20_30);
	_mm_storeu_pd(c01_11_pntr, c01_11);
	_mm_storeu_pd((c01_11_pntr+2), c21_31);
	_mm_storeu_pd(c02_12_pntr, c02_12);
	_mm_storeu_pd((c02_12_pntr+2), c22_32);
	_mm_storeu_pd(c03_13_pntr, c03_13);
	_mm_storeu_pd((c03_13_pntr+2), c23_33);
}

/* Store A so that we stride through it continuously */
static void store_a (int lda, const int K, double* a_src, double* a_dest) {
	/* For each 4xK block-row of A */
	for (int w = 0; w < K; ++w) {
		*a_dest++ = *a_src;
		*a_dest++ = *(a_src+1);
		*a_dest++ = *(a_src+2);
		*a_dest++ = *(a_src+3);
		a_src += lda;
	}
}

/* Store B so that we stride through it continuously */
static void store_b (int lda, const int K, double* b_src, double* b_dest) {
	double *b_pntr0, *b_pntr1, *b_pntr2, *b_pntr3;
	b_pntr0 = b_src;
	b_pntr1 = b_pntr0 + lda;
	b_pntr2 = b_pntr1 + lda;
	b_pntr3 = b_pntr2 + lda;
	/* For each Kx4 block-column of B */
	for (int w = 0; w < K; ++w) {
		*b_dest++ = *b_pntr0++;
		*b_dest++ = *b_pntr1++;
		*b_dest++ = *b_pntr2++;
		*b_dest++ = *b_pntr3++;
	}
}

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, double* A, double* B, double* C)
{
	double A_blocked[M*K], B_blocked[K*N];
	double *ablk_pntr, *bblk_pntr, *c;
	const int Nmax = N-3;
	int Mmax = M-3;
	int remainder = M%4;
	int i, j, p;
	/* For each column j of B */
	for (j = 0 ; j < Nmax; j += 4) {
		bblk_pntr = &B_blocked[j*K];
		store_b(lda, K, B + j*lda, &B_blocked[j*K]);
		/* For each row i of A */
		for (i = 0; i < Mmax; i += 4) {
			ablk_pntr = &A_blocked[i*K];
			if (j == 0) store_a(lda, K, A + i, &A_blocked[i*K]);
			c = C + i + j*lda;
			do_4x4(lda, K, ablk_pntr, bblk_pntr, c);
		}
	}
	/* If we have a remainder, handle it now. */
	if (remainder != 0) {
		/* For each row i of A */
		for ( ; i < M; ++i)
			/* For each column p of B */
			for (p=0; p < N; ++p)
			{
				/* Compute C(i,j) */
				double cip = C(i,p);
				for (int k = 0; k < K; ++k)
					cip += A(i,k) * B(k,p);
				C(i,p) = cip;
			}
	}
	if (N%4 != 0) {
		Mmax = M - remainder;
		/* For each column j of B */
		for ( ; j < N; ++j)
			/* For each row i of A */
			for (i=0; i < Mmax; ++i)
			{
				/* Compute C(i,j) */
				double cij = C(i,j);
				for (int k = 0; k < K; ++k)
					cij += A(i,k) * B(k,j);
				C(i,j) = cij;
			}
	}
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm (int lda, double* A, double* B, double* C)
{
	/* Accumulate block dgemms into block of C */
	for (int t = 0; t < lda; t += BLOCK_L2) {
		/* For each L2-sized block-column of B */
		for (int s = 0; s < lda; s += BLOCK_L2) {
			// int s_cntr = s * rounded_dim;
			/* For each L2-sized block-row of A */
			for (int r = 0; r < lda; r += BLOCK_L2) {
				// int r_cntr = r * BLOCK_L2;
				/* Correct block dimensions if block "goes off edge of" the matrix */
				int end_i = r + min(BLOCK_L2, lda-r);
				int end_j = s + min(BLOCK_L2, lda-s);
				int end_k = t + min(BLOCK_L2, lda-t);
				/* Accumulate block dgemms into block of C */
				for (int k = t; k < end_k; k += BLOCK_L1) {
					// int K = min (BLOCK_L1, lda-k);
					/* For each L1-sized block-column of B */
					for (int j = s; j < end_j; j += BLOCK_L1) {
						/* For each L1-sized block-row of A */
						for (int i = r; i < end_i; i += BLOCK_L1) {
							int M = min(BLOCK_L1, end_i-i);
							int N = min(BLOCK_L1, end_j-j);
							int K = min(BLOCK_L1, end_k-k);
							do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
						}
					}
				}
			}
		}
	}
}