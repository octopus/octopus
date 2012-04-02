#include <config.h>

#ifdef HAVE_CLAMDBLAS

#include <clAmdBlas.h>

void FC_FUNC_(clamdblasgetversion_low, CLAMDBLASGETVERSION_LOW)(int * major, int * minor, int * patch, int * status){
  cl_uint cl_major, cl_minor, cl_patch;

  *status = clAmdBlasGetVersion(&cl_major, &cl_minor, &cl_patch);
  *major = cl_major;
  *minor = cl_minor;
  *patch = cl_patch;
}

void FC_FUNC_(clamdblassetup_low, CLAMDBLASSETUP_LOW)(int * status){
  *status = clAmdBlasSetup();
}

void FC_FUNC_(clamdblasteardown_low, CLAMDBLASTEARDOWN_LOW)(){
  clAmdBlasTeardown();
}

void FC_FUNC_(clamdblasdtrsmex_low, CLAMDBLASDTRSMEX_LOW)(int * order,
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


  *status = clAmdBlasDtrsmEx((clAmdBlasOrder) *order, (clAmdBlasSide) *side, (clAmdBlasUplo) *uplo, 
			     (clAmdBlasTranspose) *transA, (clAmdBlasDiag) *diag,
			     (size_t) *M, (size_t) *N, *alpha, 
			     *A, (size_t) *offA, (size_t) *lda, 
			     *B, (size_t) *offB, (size_t) *ldb, 
			     1, CommandQueue, 0, NULL, NULL);
}


void FC_FUNC_(clamdblasztrsmex_low, CLAMDBLASZTRSMEX_LOW)(int * order,
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


  *status = clAmdBlasZtrsmEx((clAmdBlasOrder) *order, (clAmdBlasSide) *side, (clAmdBlasUplo) *uplo, 
			     (clAmdBlasTranspose) *transA, (clAmdBlasDiag) *diag,
			     (size_t) *M, (size_t) *N, *alpha, 
			     *A, (size_t) *offA, (size_t) *lda, 
			     *B, (size_t) *offB, (size_t) *ldb, 
			     1, CommandQueue, 0, NULL, NULL);
}

void FC_FUNC_(clamdblasdgemmex_low, CLAMDBLASDGEMMEX_LOW)(int * order,
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

  *status = clAmdBlasDgemmEx((clAmdBlasOrder) *order, (clAmdBlasTranspose) *transA, (clAmdBlasTranspose) *transB, 
			     (size_t) *M, (size_t) *N, (size_t) *K, *alpha, 
			     *A, (size_t) *offA, (size_t) *lda, 
			     *B, (size_t) *offB, (size_t) *ldb, *beta, 
			     *C, (size_t) *offC, (size_t) *ldc, 
			     1, CommandQueue, 0, NULL, NULL);
}

void FC_FUNC_(clamdblaszgemmex_low, CLAMDBLASDGEMMEX_LOW)(int * order,
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

  *status = clAmdBlasZgemmEx((clAmdBlasOrder) *order, (clAmdBlasTranspose) *transA, (clAmdBlasTranspose) *transB, 
			     (size_t) *M, (size_t) *N, (size_t) *K, *alpha, 
			     *A, (size_t) *offA, (size_t) *lda, 
			     *B, (size_t) *offB, (size_t) *ldb, *beta, 
			     *C, (size_t) *offC, (size_t) *ldc, 
			     1, CommandQueue, 0, NULL, NULL);
}

#endif
