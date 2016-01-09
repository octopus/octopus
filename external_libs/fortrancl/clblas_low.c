#include <config.h>

#ifdef HAVE_CLBLAS

#include <clBLAS.h>

void FC_FUNC_(clblasgetversion_low, CLBLASGETVERSION_LOW)(int * major, int * minor, int * patch, int * status){
  cl_uint cl_major, cl_minor, cl_patch;

  *status = clblasGetVersion(&cl_major, &cl_minor, &cl_patch);
  *major = cl_major;
  *minor = cl_minor;
  *patch = cl_patch;
}

void FC_FUNC_(clblassetup_low, CLBLASSETUP_LOW)(int * status){
  *status = clblasSetup();
}

void FC_FUNC_(clblasteardown_low, CLBLASTEARDOWN_LOW)(){
  clblasTeardown();
}

void FC_FUNC_(clblasdtrsmex_low, CLBLASDTRSMEX_LOW)(int * order,
							  int * side,
							  int * uplo,
							  int * transA,
							  int * diag,
							  cl_long * M,
							  cl_long * N,
							  double * alpha,
							  const cl_mem * A,
							  size_t * offA,
							  size_t * lda,
							  cl_mem * B,
							  size_t * offB,
							  size_t * ldb, 
							  cl_command_queue * CommandQueue, 
							  int * status){


  *status = clblasDtrsm((clblasOrder) *order, (clblasSide) *side, (clblasUplo) *uplo, 
			(clblasTranspose) *transA, (clblasDiag) *diag,
			(size_t) *M, (size_t) *N, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, 
			*B, (size_t) *offB, (size_t) *ldb, 
			1, CommandQueue, 0, NULL, NULL);
}


void FC_FUNC_(clblasztrsmex_low, CLBLASZTRSMEX_LOW)(int * order,
							  int * side,
							  int * uplo,
							  int * transA,
							  int * diag,
							  cl_long * M,
							  cl_long * N,
							  DoubleComplex * alpha,
							  const cl_mem * A,
							  size_t * offA,
							  size_t * lda,
							  cl_mem * B,
							  size_t * offB,
							  size_t * ldb, 
							  cl_command_queue * CommandQueue, 
							  int * status){


  *status = clblasZtrsm((clblasOrder) *order, (clblasSide) *side, (clblasUplo) *uplo, 
			(clblasTranspose) *transA, (clblasDiag) *diag,
			(size_t) *M, (size_t) *N, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, 
			*B, (size_t) *offB, (size_t) *ldb, 
			1, CommandQueue, 0, NULL, NULL);
}

void FC_FUNC_(clblasdgemmex_low, CLBLASDGEMMEX_LOW)(int * order,
							  int * transA, 
							  int * transB, 
							  cl_long * M,
							  cl_long * N,
							  cl_long * K,
							  double * alpha,
							  const cl_mem * A,
							  cl_long * offA,
							  cl_long * lda,
							  const cl_mem * B,
							  cl_long * offB,
							  cl_long * ldb, 
							  double * beta, 
							  cl_mem * C, 
							  cl_long * offC, 
							  cl_long * ldc, 
							  cl_command_queue * CommandQueue,
							  int * status){

  *status = clblasDgemm((clblasOrder) *order, (clblasTranspose) *transA, (clblasTranspose) *transB, 
			(size_t) *M, (size_t) *N, (size_t) *K, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, 
			*B, (size_t) *offB, (size_t) *ldb, *beta, 
			*C, (size_t) *offC, (size_t) *ldc, 
			1, CommandQueue, 0, NULL, NULL);
}

void FC_FUNC_(clblaszgemmex_low, CLBLASDGEMMEX_LOW)(int * order,
							  int * transA, 
							  int * transB, 
							  cl_long * M,
							  cl_long * N,
							  cl_long * K,
							  DoubleComplex * alpha,
							  const cl_mem * A,
							  cl_long * offA,
							  cl_long * lda,
							  const cl_mem * B,
							  cl_long * offB,
							  cl_long * ldb, 
							  DoubleComplex * beta, 
							  cl_mem * C, 
							  cl_long * offC, 
							  cl_long * ldc, 
							  cl_command_queue * CommandQueue,
							  int * status){

  *status = clblasZgemm((clblasOrder) *order, (clblasTranspose) *transA, (clblasTranspose) *transB, 
			(size_t) *M, (size_t) *N, (size_t) *K, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, 
			*B, (size_t) *offB, (size_t) *ldb, *beta, 
			*C, (size_t) *offC, (size_t) *ldc, 
			1, CommandQueue, 0, NULL, NULL);
}

void FC_FUNC_(clblasdsyrkex_low, CLBLASDSYRKEX_LOW)(int * order,
							  int * uplo, 
							  int * transA, 
							  cl_long * N,
							  cl_long * K,
							  double * alpha,
							  const cl_mem * A,
							  cl_long * offA,
							  cl_long * lda,
							  double * beta, 
							  cl_mem * C, 
							  cl_long * offC, 
							  cl_long * ldc, 
							  cl_command_queue * CommandQueue,
							  int * status){

  *status = clblasDsyrk((clblasOrder) *order, (clblasUplo) *uplo, (clblasTranspose) *transA, 
			(size_t) *N, (size_t) *K, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, *beta, 
			*C, (size_t) *offC, (size_t) *ldc, 
			1, CommandQueue, 0, NULL, NULL);
}

void FC_FUNC_(clblaszherkex_low, CLBLASZHERKEX_LOW)(int * order,
							  int * uplo, 
							  int * transA, 
							  cl_long * N,
							  cl_long * K,
							  double * alpha,
							  const cl_mem * A,
							  cl_long * offA,
							  cl_long * lda,
							  double * beta, 
							  cl_mem * C, 
							  cl_long * offC, 
							  cl_long * ldc, 
							  cl_command_queue * CommandQueue,
							  int * status){

  *status = clblasZherk((clblasOrder) *order, (clblasUplo) *uplo, (clblasTranspose) *transA, 
			(size_t) *N, (size_t) *K, *alpha, 
			*A, (size_t) *offA, (size_t) *lda, *beta, 
			*C, (size_t) *offC, (size_t) *ldc, 
			1, CommandQueue, 0, NULL, NULL);
}


clblasStatus FC_FUNC_(clblasddot_low, CLBLASDDOT_LOW)(cl_long * N,
						      cl_mem * dotProduct,
						      cl_long * offDP,
						      const cl_mem *X,
						      cl_long * offx,
						      int *incx,
						      const cl_mem * Y,
						      cl_long * offy,
						      int * incy,
						      cl_mem * scratchBuff,
						      cl_command_queue * CommandQueue,
						      int * status){

  
  *status = clblasDdot((size_t) *N, *dotProduct, (size_t) *offDP, *X, (size_t) *offx, *incx,
		       *Y, (size_t) *offy, *incy, *scratchBuff, 1, CommandQueue, 0, NULL, NULL);

}


clblasStatus FC_FUNC_(clblaszdotc_low, CLBLASZDOTC_LOW)(cl_long * N,
							cl_mem * dotProduct,
							cl_long * offDP,
							const cl_mem *X,
							cl_long * offx,
							int *incx,
							const cl_mem * Y,
							cl_long * offy,
							int * incy,
							cl_mem * scratchBuff,
							cl_command_queue * CommandQueue,
							int * status){

  
  *status = clblasZdotc((size_t) *N, *dotProduct, (size_t) *offDP, *X, (size_t) *offx, *incx,
			*Y, (size_t) *offy, *incy, *scratchBuff, 1, CommandQueue, 0, NULL, NULL);

}


clblasStatus FC_FUNC_(clblasdnrm2_low, CLBLASDNRM2_LOW)(cl_long * N,
							cl_mem * NRM2,
							cl_long * offNRM2,
							const cl_mem *X,
							cl_long * offx,
							int *incx,
							cl_mem * scratchBuff,
							cl_command_queue * CommandQueue,
							int * status){

  
  *status = clblasDnrm2((size_t) *N, *NRM2, (size_t) *offNRM2, *X, (size_t) *offx, *incx,
			*scratchBuff, 1, CommandQueue, 0, NULL, NULL);

}

clblasStatus FC_FUNC_(clblasdznrm2_low, CLBLASDZNRM2_LOW)(cl_long * N,
							  cl_mem * NRM2,
							  cl_long * offNRM2,
							  const cl_mem *X,
							  cl_long * offx,
							  int *incx,
							  cl_mem * scratchBuff,
							  cl_command_queue * CommandQueue,
							  int * status){
  
  
  *status = clblasDznrm2((size_t) *N, *NRM2, (size_t) *offNRM2, *X, (size_t) *offx, *incx,
			 *scratchBuff, 1, CommandQueue, 0, NULL, NULL);
  
}

#endif
