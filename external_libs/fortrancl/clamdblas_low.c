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
							  cl_command_queue * commandQueue, 
							  int * status){


  *status = clAmdBlasDtrsmEx((clAmdBlasOrder) *order, (clAmdBlasSide) *side, (clAmdBlasUplo) *uplo, 
			     (clAmdBlasTranspose) *transA, (clAmdBlasDiag) *diag,
			     (size_t) *M, (size_t) *N, *alpha, 
			     *A, (size_t) *offA, (size_t) *lda, 
			     *B, (size_t) *offB, (size_t) *ldb, 
			     1, commandQueue, 
			     0, NULL, NULL);
}

#endif
