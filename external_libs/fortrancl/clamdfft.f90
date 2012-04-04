module clAmdFft
  use cl

  implicit none

  private

  integer, parameter ::                                                           &
    CLFFT_INVALID_GLOBAL_WORK_SIZE         = CL_INVALID_GLOBAL_WORK_SIZE,         &
    CLFFT_INVALID_MIP_LEVEL                = CL_INVALID_MIP_LEVEL,                &
    CLFFT_INVALID_BUFFER_SIZE              = CL_INVALID_BUFFER_SIZE,              &
    CLFFT_INVALID_GL_OBJECT                = CL_INVALID_GL_OBJECT,                &
    CLFFT_INVALID_OPERATION                = CL_INVALID_OPERATION,                &
    CLFFT_INVALID_EVENT                    = CL_INVALID_EVENT,                    &
    CLFFT_INVALID_EVENT_WAIT_LIST          = CL_INVALID_EVENT_WAIT_LIST,          &
    CLFFT_INVALID_GLOBAL_OFFSET            = CL_INVALID_GLOBAL_OFFSET,            &
    CLFFT_INVALID_WORK_ITEM_SIZE           = CL_INVALID_WORK_ITEM_SIZE,           &
    CLFFT_INVALID_WORK_GROUP_SIZE          = CL_INVALID_WORK_GROUP_SIZE,          &
    CLFFT_INVALID_WORK_DIMENSION           = CL_INVALID_WORK_DIMENSION,           &
    CLFFT_INVALID_KERNEL_ARGS              = CL_INVALID_KERNEL_ARGS,              &
    CLFFT_INVALID_ARG_SIZE                 = CL_INVALID_ARG_SIZE,                 &
    CLFFT_INVALID_ARG_VALUE                = CL_INVALID_ARG_VALUE,                &
    CLFFT_INVALID_ARG_INDEX                = CL_INVALID_ARG_INDEX,                &
    CLFFT_INVALID_KERNEL                   = CL_INVALID_KERNEL,                   &
    CLFFT_INVALID_KERNEL_DEFINITION        = CL_INVALID_KERNEL_DEFINITION,        &
    CLFFT_INVALID_KERNEL_NAME              = CL_INVALID_KERNEL_NAME,              &
    CLFFT_INVALID_PROGRAM_EXECUTABLE       = CL_INVALID_PROGRAM_EXECUTABLE,       &
    CLFFT_INVALID_PROGRAM                  = CL_INVALID_PROGRAM,                  &
    CLFFT_INVALID_BUILD_OPTIONS            = CL_INVALID_BUILD_OPTIONS,            &
    CLFFT_INVALID_BINARY                   = CL_INVALID_BINARY,                   &
    CLFFT_INVALID_SAMPLER                  = CL_INVALID_SAMPLER

  integer, parameter ::                                                           &
    CLFFT_INVALID_IMAGE_SIZE               = CL_INVALID_IMAGE_SIZE,               &
    CLFFT_INVALID_IMAGE_FORMAT_DESCRIPTOR  = CL_INVALID_IMAGE_FORMAT_DESCRIPTOR,  &
    CLFFT_INVALID_MEM_OBJECT               = CL_INVALID_MEM_OBJECT,               &
    CLFFT_INVALID_HOST_PTR                 = CL_INVALID_HOST_PTR,                 &
    CLFFT_INVALID_COMMAND_QUEUE            = CL_INVALID_COMMAND_QUEUE,            &
    CLFFT_INVALID_QUEUE_PROPERTIES         = CL_INVALID_QUEUE_PROPERTIES,         &
    CLFFT_INVALID_CONTEXT                  = CL_INVALID_CONTEXT,                  &
    CLFFT_INVALID_DEVICE                   = CL_INVALID_DEVICE,                   &
    CLFFT_INVALID_PLATFORM                 = CL_INVALID_PLATFORM,                 &
    CLFFT_INVALID_DEVICE_TYPE              = CL_INVALID_DEVICE_TYPE,              &
    CLFFT_INVALID_VALUE                    = CL_INVALID_VALUE,                    &
    CLFFT_MAP_FAILURE                      = CL_MAP_FAILURE,                      &
    CLFFT_BUILD_PROGRAM_FAILURE            = CL_BUILD_PROGRAM_FAILURE,            &
    CLFFT_IMAGE_FORMAT_NOT_SUPPORTED       = CL_IMAGE_FORMAT_NOT_SUPPORTED,       &
    CLFFT_IMAGE_FORMAT_MISMATCH            = CL_IMAGE_FORMAT_MISMATCH,            &
    CLFFT_MEM_COPY_OVERLAP                 = CL_MEM_COPY_OVERLAP,                 &
    CLFFT_PROFILING_INFO_NOT_AVAILABLE     = CL_PROFILING_INFO_NOT_AVAILABLE,     &
    CLFFT_OUT_OF_HOST_MEMORY               = CL_OUT_OF_HOST_MEMORY,               &
    CLFFT_OUT_OF_RESOURCES                 = CL_OUT_OF_RESOURCES,                 &
    CLFFT_MEM_OBJECT_ALLOCATION_FAILURE    = CL_MEM_OBJECT_ALLOCATION_FAILURE,    &
    CLFFT_COMPILER_NOT_AVAILABLE           = CL_COMPILER_NOT_AVAILABLE,           &
    CLFFT_DEVICE_NOT_AVAILABLE             = CL_DEVICE_NOT_AVAILABLE,             &
    CLFFT_DEVICE_NOT_FOUND                 = CL_DEVICE_NOT_FOUND,                 &
    CLFFT_SUCCESS                          = CL_SUCCESS

  integer, parameter ::                                                           &
    CLFFT_BUGCHECK                         = 4*1024    ,                          &
    CLFFT_NOTIMPLEMENTED                   = 4*1024 + 1,                          &
    CLFFT_FILE_NOT_FOUND                   = 4*1024 + 2,                          &
    CLFFT_FILE_CREATE_FAILURE              = 4*1024 + 3,                          &
    CLFFT_VERSION_MISMATCH                 = 4*1024 + 4,                          &   
    CLFFT_INVALID_PLAN                     = 4*1024 + 5,                          &  
    CLFFT_DEVICE_NO_DOUBLE                 = 4*1024 + 6,                          &     
    CLFFT_ENDSTATUS                        = 4*1024 + 7

  integer, parameter ::                &
    CLFFT_1D                      = 1, &
    CLFFT_2D                      = 2, &
    CLFFT_3D                      = 3, &
    ENDDIMENSION                  = 4

  integer, parameter ::                &
    CLFFT_COMPLEX_INTERLEAVED     = 1, &
    CLFFT_COMPLEX_PLANAR          = 2, &
    CLFFT_HERMITIAN_INTERLEAVED   = 4, &
    CLFFT_HERMITIAN_PLANAR        = 5, &
    CLFFT_REAL                    = 6, &
    ENDLAYOUT                     = 7

  integer, parameter ::                &
    CLFFT_SINGLE                  = 1, &
    CLFFT_DOUBLE                  = 2, &
    CLFFT_SINGLE_FAST             = 3, &
    CLFFT_DOUBLE_FAST             = 4, &
    ENDPRECISION                  = 5

  integer, parameter ::                &
    CLFFT_INPLACE                 = 1, &
    CLFFT_OUTOFPLACE              = 2, &
    ENDPLACE                      = 3

  interface clAmdFftGetVersion 
    subroutine clamdfftgetversion_low(major, minor, patch, status)
      implicit none

      integer, intent(out) :: major
      integer, intent(out) :: minor
      integer, intent(out) :: patch
      integer, intent(out) :: status
    end subroutine clamdfftgetversion_low
  end interface clAmdFftGetVersion

  interface clAmdFftTeardown
    subroutine clamdfftteardown_low()
    end subroutine clamdfftteardown_low
  end interface clAmdFftTeardown

end module clAmdFft
