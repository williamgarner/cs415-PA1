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

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */    
void square_dgemm (const int n, double*  A, double* B, double* restrict C)
{

//	if(n < 191)
//	{
//		#pragma loop_count min(31), max(769), avg(345)
//		#pragma block_loop factor(1024) level(1)
//		for (int j = 0; j < n; ++j)
//			for( int k = 0; k < n; k++ )
//				#pragma vector unaligned
//				for (int i = 0; i < n; ++i)
//					C[i+ j*n] += A[i+ k*n] * B[k+ j*n];
//	}
//	else
//	{
			double T[n*n];
//			#pragma loop_count min(31), max(769), avg(345)
			#pragma block_loop level(2) factor(16)
				for(int i = 0; i < n; ++i)
					#pragma vector unaligned
					for(int j = 0; j < n; ++j)
							T[i*n + j] = A[j*n + i];


//			#pragma loop_count min(31), max(769), avg(345)
//			#pragma block_loop factor(256)  level(2)
//				#pragma block_loop level(1:2)
				for (int j = 0; j < n; ++j)
					for (int i = 0; i < n; ++i)
						#pragma vector unaligned
						for( int k = 0; k < n; ++k)
								C[i+ j*n] += T[k+ i*n] * B[k+ j*n];
//	}
}
