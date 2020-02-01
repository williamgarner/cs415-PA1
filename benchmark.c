//#include <stdlib.h> // For: exit, drand48, malloc, free, NULL, EXIT_FAILURE
//#include <stdio.h>  // For: perror
//#include <string.h> // For: memset
//
//#include <float.h>  // For: DBL_EPSILON
//#include <math.h>   // For: fabs
//
//#ifdef GETTIMEOFDAY
//#include <sys/time.h> // For struct timeval, gettimeofday
//#else
//#include <time.h> // For struct timespec, clock_gettime, CLOCK_MONOTONIC
//#endif
//
//#define MAX_SPEED 42.9  // defining Bridges Max Gflops/s per core with peak TurboBoost Frequency
//
///* reference_dgemm wraps a call to the BLAS-3 routine DGEMM, via the standard FORTRAN interface - hence the reference semantics. */
//#define DGEMM dgemm_
//extern void DGEMM (char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
//void reference_dgemm (int N, double ALPHA, double* A, double* B, double* C)
//{
//  char TRANSA = 'N';
//  char TRANSB = 'N';
//  int M = N;
//  int K = N;
//  double BETA = 1.;
//  int LDA = N;
//  int LDB = N;
//  int LDC = N;
//  DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
//}
//
///* Your function must have the following signature: */
//extern const char* dgemm_desc;
//extern void square_dgemm (int, double*, double*, double*);
//
//double wall_time ()
//{
//#ifdef GETTIMEOFDAY
//  struct timeval t;
//  gettimeofday (&t, NULL);
//  return 1.*t.tv_sec + 1.e-6*t.tv_usec;
//#else
//  struct timespec t;
//  clock_gettime (CLOCK_MONOTONIC, &t);
//  return 1.*t.tv_sec + 1.e-9*t.tv_nsec;
//#endif
//}
//
//void die (const char* message)
//{
//  perror (message);
//  exit (EXIT_FAILURE);
//}
//
//void fill (double* p, int n)
//{
//  for (int i = 0; i < n; ++i)
//    p[i] = 2 * drand48() - 1; // Uniformly distributed over [-1, 1]
//}
//
//void absolute_value (double *p, int n)
//{
//  for (int i = 0; i < n; ++i)
//    p[i] = fabs (p[i]);
//}
//
///* The benchmarking program */
//int main (int argc, char **argv)
//{
//  printf ("Description:\t%s\n\n", dgemm_desc);
//
//  /* Test sizes should highlight performance dips at multiples of certain powers-of-two */
//
//  int test_sizes[] =
//
//  /* Multiples-of-32, +/- 1. Currently commented. */
//  /* {31,32,33,63,64,65,95,96,97,127,128,129,159,160,161,191,192,193,223,224,225,255,256,257,287,288,289,319,320,321,351,352,353,383,384,385,415,416,417,447,448,449,479,480,481,511,512,513,543,544,545,575,576,577,607,608,609,639,640,641,671,672,673,703,704,705,735,736,737,767,768,769,799,800,801,831,832,833,863,864,865,895,896,897,927,928,929,959,960,961,991,992,993,1023,1024,1025}; */
//
//  /* A representative subset of the first list. Currently uncommented. */
//  { 31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
//    319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767, 768, 769 };
//
//  int nsizes = sizeof(test_sizes)/sizeof(test_sizes[0]);
//
//  /* assume last size is also the largest size */
//  int nmax = test_sizes[nsizes-1];
//
//  /* allocate memory for all problems */
//  double* buf = NULL;
//  buf = (double*) malloc (3 * nmax * nmax * sizeof(double));
//  if (buf == NULL) die ("failed to allocate largest problem size");
//
//  double Mflops_s[nsizes],per[nsizes],aveper,grade;
//
//  /* For each test size */
//  for (int isize = 0; isize < sizeof(test_sizes)/sizeof(test_sizes[0]); ++isize)
//  {
//    /* Create and fill 3 random matrices A,B,C*/
//    int n = test_sizes[isize];
//
//    double* A = buf + 0;
//    double* B = A + nmax*nmax;
//    double* C = B + nmax*nmax;
//
//    fill (A, n*n);
//    fill (B, n*n);
//    fill (C, n*n);
//
//    /* Measure performance (in Gflops/s). */
//
//    /* Time a "sufficiently long" sequence of calls to reduce noise */
//    double Gflops_s, seconds = -1.0;
//    double timeout = 0.1; // "sufficiently long" := at least 1/10 second.
//    for (int n_iterations = 1; seconds < timeout; n_iterations *= 2)
//    {
//      /* Warm-up */
//      square_dgemm (n, A, B, C);
//
//      /* Benchmark n_iterations runs of square_dgemm */
//      seconds = -wall_time();
//      for (int it = 0; it < n_iterations; ++it)
//	square_dgemm (n, A, B, C);
//      seconds += wall_time();
//
//      /*  compute Gflop/s rate */
//      Gflops_s = 2.e-9 * n_iterations * n * n * n / seconds;
//    }
//
//    /* Storing Mflop rate and calculating percentage of peak */
//    Mflops_s[isize] = Gflops_s*1000;
//    per[isize] = Gflops_s*100/MAX_SPEED;
//
//    printf ("Size: %d\tMflop/s: %8g\tPercentage:%6.2lf\n", n, Mflops_s[isize],per[isize]);
//
//    /* Ensure that error does not exceed the theoretical error bound. */
//
//    /* C := A * B, computed with square_dgemm */
//    memset (C, 0, n * n * sizeof(double));
//    square_dgemm (n, A, B, C);
//
//    /* Do not explicitly check that A and B were unmodified on square_dgemm exit
//     *  - if they were, the following will most likely detect it:
//     * C := C - A * B, computed with reference_dgemm */
//    reference_dgemm(n, -1., A, B, C);
//
//    /* A := |A|, B := |B|, C := |C| */
//    absolute_value (A, n * n);
//    absolute_value (B, n * n);
//    absolute_value (C, n * n);
//
//    /* C := |C| - 3 * e_mach * n * |A| * |B|, computed with reference_dgemm */
//    reference_dgemm (n, -3.*DBL_EPSILON*n, A, B, C);
//
//    /* If any element in C is positive, then something went wrong in square_dgemm */
//    for (int i = 0; i < n * n; ++i)
//      if (C[i] > 0)
//	  {
//      	printf("C[i]: %lf\n", C[i]);
//      	printf("i: %d\n", i);
//      	printf("\n");
//	  fflush(stdout);
//	  die("*** FAILURE *** Error in matrix multiply exceeds componentwise error bounds.\n" );
//	  }
//
//  }
//
//  /* Calculating average percentage of peak reached by algorithm */
//  aveper=0;
//  for (int i=0; i<nsizes;i++)
//    aveper+= per[i];
//  aveper/=nsizes*1.0;
//
//  /* Assigning grade based on average percentage reached (50% gets 75; 80% gets 100; rest distributed proportionally) */
//  if (aveper > 80) grade = 100.0;
//  else if (aveper > 50) grade = (aveper-50)*0.25*100.0/30.0 + 75.0;
//  else
//     grade = aveper * 2 * 0.75;
//
//  /* Printing average percentage and grade to screen */
//  printf("Average percentage of Peak = %g\nGrade = %g\n",aveper,grade);
//
//  free (buf);
//
//  return 0;
//}





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
			bx0, bx1, bx2, bx3,
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

		bx0 = _mm_loaddup_pd(b++);
		bx1 = _mm_loaddup_pd(b++);
		bx2 = _mm_loaddup_pd(b++);
		bx3 = _mm_loaddup_pd(b++);

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