module clAmdBlas
  use cl

  implicit none

  private

  public ::                  &
    clAmdBlasGetVersion,     &
    clAmdBlasSetup,          &
    clAmdBlasTeardown,       &
    clAmdBlasDtrsmEx

  integer, public, parameter ::     &
    clAmdBlasRowMajor        = 0,   &
    clAmdBlasColumnMajor     = 1

  integer, public, parameter ::     &
    clAmdBlasNoTrans         = 0,   &
    clAmdBlasTrans           = 1,   &
    clAmdBlasConjTrans       = 2

  integer, public, parameter ::     &
    clAmdBlasUpper           = 0,   &
    clAmdBlasLower           = 1

  integer, public, parameter ::     &
    clAmdBlasUnit            = 0,   &
    clAmdBlasNonUnit         = 1

  integer, public, parameter ::     &
    clAmdBlasLeft            = 0,   &
    clAmdBlasRight           = 1

  integer, public, parameter ::                                  &
    clAmdBlasSuccess               = CL_SUCCESS,                 &
    clAmdBlasInvalidValue          = CL_INVALID_VALUE,           &
    clAmdBlasInvalidCommandQueue   = CL_INVALID_COMMAND_QUEUE,   &
    clAmdBlasInvalidContext        = CL_INVALID_CONTEXT,         &
    clAmdBlasInvalidMemObject      = CL_INVALID_MEM_OBJECT,      &
    clAmdBlasInvalidDevice         = CL_INVALID_DEVICE,          &
    clAmdBlasInvalidEventWaitList  = CL_INVALID_EVENT_WAIT_LIST, &
    clAmdBlasOutOfResources        = CL_OUT_OF_RESOURCES,        &
    clAmdBlasOutOfHostMemory       = CL_OUT_OF_HOST_MEMORY,      &
    clAmdBlasInvalidOperation      = CL_INVALID_OPERATION,       &
    clAmdBlasCompilerNotAvailable  = CL_COMPILER_NOT_AVAILABLE,  &
    clAmdBlasBuildProgramFailure   = CL_BUILD_PROGRAM_FAILURE

  integer, public, parameter ::             &
    clAmdBlasNotImplemented        = -1024, &
    clAmdBlasNotInitialized        = -1023, &
    clAmdBlasInvalidMatA           = -1022, &
    clAmdBlasInvalidMatB           = -1021, &
    clAmdBlasInvalidMatC           = -1020, &
    clAmdBlasInvalidVecX           = -1019, &
    clAmdBlasInvalidVecY           = -1018, &
    clAmdBlasInvalidDim            = -1017, &
    clAmdBlasInvalidLeadDimA       = -1016, &
    clAmdBlasInvalidLeadDimB       = -1015, &
    clAmdBlasInvalidLeadDimC       = -1014, &
    clAmdBlasInvalidIncX           = -1013, &
    clAmdBlasInvalidIncY           = -1012, &
    clAmdBlasInsufficientMemMatA   = -1011, &
    clAmdBlasInsufficientMemMatB   = -1010, &
    clAmdBlasInsufficientMemMatC   = -1009, &
    clAmdBlasInsufficientMemVecX   = -1008, &
    clAmdBlasInsufficientMemVecY   = -1007

  ! SUPPORT FUNCTIONS

  interface clAmdBlasGetVersion 
    subroutine clamdblasgetversion_low(major, minor, patch, status)
      implicit none

      integer, intent(out) :: major
      integer, intent(out) :: minor
      integer, intent(out) :: patch
      integer, intent(out) :: status
    end subroutine clamdblasgetversion_low
  end interface clAmdBlasGetVersion

  interface clAmdBlasSetup
    subroutine clamdblassetup_low(status)
      implicit none

      integer, intent(out) :: status
    end subroutine clamdblassetup_low
  end interface clAmdBlasSetup

  interface clAmdBlasTeardown
    subroutine clamdblasteardown_low()
    end subroutine clamdblasteardown_low
  end interface clAmdBlasTeardown

  interface clAmdBlasDtrsmEx
    subroutine clamdblasdtrsmex_low(order, side, uplo, transA, diag, M, N, alpha, A, offA, lda, B, offB, ldb, commandQueue, status)
      use cl

      implicit none

      integer,                intent(in)    :: order
      integer,                intent(in)    :: side
      integer,                intent(in)    :: uplo
      integer,                intent(in)    :: transA
      integer,                intent(in)    :: diag
      integer(8),             intent(in)    :: M
      integer(8),             intent(in)    :: N
      real(8),                intent(in)    :: alpha
      type(cl_mem),           intent(inout) :: A
      integer(8),             intent(in)    :: offA
      integer(8),             intent(in)    :: lda
      type(cl_mem),           intent(inout) :: B
      integer(8),             intent(in)    :: offB
      integer(8),             intent(in)    :: ldb 
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(in)    :: status
    end subroutine clamdblasdtrsmex_low
  end interface clAmdBlasDtrsmEx
  
end module clAmdBlas
