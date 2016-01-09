module clblas
  use cl

  implicit none

  private

  public ::                  &
    clblasGetVersion,     &
    clblasSetup,          &
    clblasTeardown,       &
    clblasDtrsmEx,        &
    clblasZtrsmEx,        &
    clblasDgemmEx,        &
    clblasZgemmEx,        &
    clblasDsyrkEx,        &
    clblasZherkEx,        &
    clblasDdot,           &
    clblasZdotc,          &
    clblasDnrm2,          &
    clblasDznrm2

  integer, public, parameter ::     &
    clblasRowMajor        = 0,   &
    clblasColumnMajor     = 1

  integer, public, parameter ::     &
    clblasNoTrans         = 0,   &
    clblasTrans           = 1,   &
    clblasConjTrans       = 2

  integer, public, parameter ::     &
    clblasUpper           = 0,   &
    clblasLower           = 1

  integer, public, parameter ::     &
    clblasUnit            = 0,   &
    clblasNonUnit         = 1

  integer, public, parameter ::     &
    clblasLeft            = 0,   &
    clblasRight           = 1

  integer, public, parameter ::                                  &
    clblasSuccess               = CL_SUCCESS,                 &
    clblasInvalidValue          = CL_INVALID_VALUE,           &
    clblasInvalidCommandQueue   = CL_INVALID_COMMAND_QUEUE,   &
    clblasInvalidContext        = CL_INVALID_CONTEXT,         &
    clblasInvalidMemObject      = CL_INVALID_MEM_OBJECT,      &
    clblasInvalidDevice         = CL_INVALID_DEVICE,          &
    clblasInvalidEventWaitList  = CL_INVALID_EVENT_WAIT_LIST, &
    clblasOutOfResources        = CL_OUT_OF_RESOURCES,        &
    clblasOutOfHostMemory       = CL_OUT_OF_HOST_MEMORY,      &
    clblasInvalidOperation      = CL_INVALID_OPERATION,       &
    clblasCompilerNotAvailable  = CL_COMPILER_NOT_AVAILABLE,  &
    clblasBuildProgramFailure   = CL_BUILD_PROGRAM_FAILURE

  integer, public, parameter ::             &
    clblasNotImplemented        = -1024, &
    clblasNotInitialized        = -1023, &
    clblasInvalidMatA           = -1022, &
    clblasInvalidMatB           = -1021, &
    clblasInvalidMatC           = -1020, &
    clblasInvalidVecX           = -1019, &
    clblasInvalidVecY           = -1018, &
    clblasInvalidDim            = -1017, &
    clblasInvalidLeadDimA       = -1016, &
    clblasInvalidLeadDimB       = -1015, &
    clblasInvalidLeadDimC       = -1014, &
    clblasInvalidIncX           = -1013, &
    clblasInvalidIncY           = -1012, &
    clblasInsufficientMemMatA   = -1011, &
    clblasInsufficientMemMatB   = -1010, &
    clblasInsufficientMemMatC   = -1009, &
    clblasInsufficientMemVecX   = -1008, &
    clblasInsufficientMemVecY   = -1007

  ! SUPPORT FUNCTIONS

  interface clblasGetVersion 
    subroutine clblasgetversion_low(major, minor, patch, status)
      implicit none

      integer, intent(out) :: major
      integer, intent(out) :: minor
      integer, intent(out) :: patch
      integer, intent(out) :: status
    end subroutine clblasgetversion_low
  end interface clblasGetVersion

  interface clblasSetup
    subroutine clblassetup_low(status)
      implicit none

      integer, intent(out) :: status
    end subroutine clblassetup_low
  end interface clblasSetup

  interface clblasTeardown
    subroutine clblasteardown_low()
    end subroutine clblasteardown_low
  end interface clblasTeardown

  ! -------------------------------------------------

  interface clblasDtrsmEx
    subroutine clblasdtrsmex_low(order, side, uplo, transA, diag, M, N, alpha, A, offA, lda, B, offB, ldb, commandQueue, status)
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
      integer,                intent(out)   :: status
    end subroutine clblasdtrsmex_low
  end interface clblasDtrsmEx
  
  ! -------------------------------------------------

  interface clblasZtrsmEx
    subroutine clblasztrsmex_low(order, side, uplo, transA, diag, M, N, alpha, A, offA, lda, B, offB, ldb, commandQueue, status)
      use cl

      implicit none

      integer,                intent(in)    :: order
      integer,                intent(in)    :: side
      integer,                intent(in)    :: uplo
      integer,                intent(in)    :: transA
      integer,                intent(in)    :: diag
      integer(8),             intent(in)    :: M
      integer(8),             intent(in)    :: N
      complex(8),             intent(in)    :: alpha
      type(cl_mem),           intent(inout) :: A
      integer(8),             intent(in)    :: offA
      integer(8),             intent(in)    :: lda
      type(cl_mem),           intent(inout) :: B
      integer(8),             intent(in)    :: offB
      integer(8),             intent(in)    :: ldb 
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(out)   :: status
    end subroutine clblasztrsmex_low
  end interface clblasZtrsmEx

  ! -------------------------------------------------

  interface clblasDgemmEx
    subroutine clblasDgemmEx_low(order, transA, transB, M, N, K, &
      alpha, A, offA, lda, B, offB, ldb, beta, C, offC, ldc, commandQueue, status)
      use cl

      implicit none

      integer,                intent(in)    :: order
      integer,                intent(in)    :: transA
      integer,                intent(in)    :: transB
      integer(8),             intent(in)    :: M
      integer(8),             intent(in)    :: N
      integer(8),             intent(in)    :: K
      real(8),                intent(in)    :: alpha
      type(cl_mem),           intent(in)    :: A
      integer(8),             intent(in)    :: offA
      integer(8),             intent(in)    :: lda
      type(cl_mem),           intent(in)    :: B
      integer(8),             intent(in)    :: offB
      integer(8),             intent(in)    :: ldb
      real(8),                intent(in)    :: beta
      type(cl_mem),           intent(inout) :: C
      integer(8),             intent(in)    :: offC
      integer(8),             intent(in)    :: ldc
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(out)   :: status
    end subroutine clblasDgemmEx_low
  end interface clblasDgemmEx

  ! -------------------------------------------------

  interface clblasZgemmEx
    subroutine clblasZgemmEx_low(order, transA, transB, M, N, K, &
      alpha, A, offA, lda, B, offB, ldb, beta, C, offC, ldc, commandQueue, status)
      use cl

      implicit none

      integer,                intent(in)    :: order
      integer,                intent(in)    :: transA
      integer,                intent(in)    :: transB
      integer(8),             intent(in)    :: M
      integer(8),             intent(in)    :: N
      integer(8),             intent(in)    :: K
      complex(8),             intent(in)    :: alpha
      type(cl_mem),           intent(in)    :: A
      integer(8),             intent(in)    :: offA
      integer(8),             intent(in)    :: lda
      type(cl_mem),           intent(in)    :: B
      integer(8),             intent(in)    :: offB
      integer(8),             intent(in)    :: ldb
      complex(8),             intent(in)    :: beta
      type(cl_mem),           intent(inout) :: C
      integer(8),             intent(in)    :: offC
      integer(8),             intent(in)    :: ldc
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(out)   :: status
    end subroutine clblasZgemmEx_low
  end interface clblasZgemmEx

  ! -------------------------------------------------

  interface clblasDsyrkEx
    subroutine clblasDsyrkEx_low(order, uplo, transA, N, K, &
      alpha, A, offA, lda, beta, C, offC, ldc, commandQueue, status)
      use cl

      implicit none

      integer,                intent(in)    :: order
      integer,                intent(in)    :: uplo
      integer,                intent(in)    :: transA
      integer(8),             intent(in)    :: N
      integer(8),             intent(in)    :: K
      real(8),                intent(in)    :: alpha
      type(cl_mem),           intent(in)    :: A
      integer(8),             intent(in)    :: offA
      integer(8),             intent(in)    :: lda
      real(8),                intent(in)    :: beta
      type(cl_mem),           intent(inout) :: C
      integer(8),             intent(in)    :: offC
      integer(8),             intent(in)    :: ldc
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(out)   :: status
    end subroutine clblasDsyrkEx_low
  end interface clblasDsyrkEx
  ! -------------------------------------------------

  interface clblasZherkEx
    subroutine clblasZherkEx_low(order, uplo, transA, N, K, &
      alpha, A, offA, lda, beta, C, offC, ldc, commandQueue, status)
      use cl

      implicit none

      integer,                intent(in)    :: order
      integer,                intent(in)    :: uplo
      integer,                intent(in)    :: transA
      integer(8),             intent(in)    :: N
      integer(8),             intent(in)    :: K
      real(8),                intent(in)    :: alpha
      type(cl_mem),           intent(in)    :: A
      integer(8),             intent(in)    :: offA
      integer(8),             intent(in)    :: lda
      real(8),                intent(in)    :: beta
      type(cl_mem),           intent(inout) :: C
      integer(8),             intent(in)    :: offC
      integer(8),             intent(in)    :: ldc
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(out)   :: status
    end subroutine clblasZherkEx_low
  end interface clblasZherkEx


  ! --------------------------------------------------

  interface clblasDdot
    subroutine clblasDdot_low(N, dotProduct, offDP, X, offx, incx, Y, offY, incy, scratchBuff, CommandQueue, status)
      use cl

      implicit none

      integer(8),             intent(in)    :: N
      type(cl_mem),           intent(inout) :: dotProduct
      integer(8),             intent(in)    :: offDP
      type(cl_mem),           intent(in)    :: X
      integer(8),             intent(in)    :: offX
      integer,                intent(in)    :: incx
      type(cl_mem),           intent(in)    :: Y
      integer(8),             intent(in)    :: offy
      integer,                intent(in)    :: incy
      type(cl_mem),           intent(inout) :: scratchBuff
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(out)   :: status
    end subroutine clblasDdot_low
  end interface clblasDdot


  ! --------------------------------------------------

  interface clblasZdotc
    subroutine clblasZdotc_low(N, dotProduct, offDP, X, offx, incx, Y, offY, incy, scratchBuff, CommandQueue, status)
      use cl

      implicit none

      integer(8),             intent(in)    :: N
      type(cl_mem),           intent(inout) :: dotProduct
      integer(8),             intent(in)    :: offDP
      type(cl_mem),           intent(in)    :: X
      integer(8),             intent(in)    :: offX
      integer,                intent(in)    :: incx
      type(cl_mem),           intent(in)    :: Y
      integer(8),             intent(in)    :: offy
      integer,                intent(in)    :: incy
      type(cl_mem),           intent(inout) :: scratchBuff
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(out)   :: status
    end subroutine clblasZdotc_low
  end interface clblasZdotc

  ! --------------------------------------------------

  interface clblasDnrm2
    subroutine clblasDnrm2_low(N, NRM2, offNRM2, X, offx, incx, scratchBuff, CommandQueue, status)
      use cl

      implicit none

      integer(8),             intent(in)    :: N
      type(cl_mem),           intent(inout) :: NRM2
      integer(8),             intent(in)    :: offNRM2
      type(cl_mem),           intent(in)    :: X
      integer(8),             intent(in)    :: offX
      integer,                intent(in)    :: incx
      type(cl_mem),           intent(inout) :: scratchBuff
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(out)   :: status
    end subroutine clblasDnrm2_low
  end interface clblasDnrm2

  ! --------------------------------------------------

  interface clblasDznrm2
    subroutine clblasDznrm2_low(N, NRM2, offNRM2, X, offx, incx, scratchBuff, CommandQueue, status)
      use cl

      implicit none

      integer(8),             intent(in)    :: N
      type(cl_mem),           intent(inout) :: NRM2
      integer(8),             intent(in)    :: offNRM2
      type(cl_mem),           intent(in)    :: X
      integer(8),             intent(in)    :: offX
      integer,                intent(in)    :: incx
      type(cl_mem),           intent(inout) :: scratchBuff
      type(cl_command_queue), intent(inout) :: CommandQueue 
      integer,                intent(out)   :: status
    end subroutine clblasDznrm2_low
  end interface clblasDznrm2

end module clblas
