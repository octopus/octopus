!! Copyright (C) 2011 J. Alberdi, P. Garcia Risueño
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: pfft.F90 7027 2011-02-08 22:47:20Z pablogr $

#include "global.h"

#define PFFT_MAX 10
#define FFT_NULL -1

#if defined(SINGLE_PRECISION)
#  define PDFFT(x) dpfftf_ ## x
#else
#  define PDFFT(x) dpfft_ ## x
#endif

module pfft_m
  use c_pointer_m
  use datasets_m
  use global_m
  use lalg_basic_m
  use loct_math_m
  use messages_m
  use parser_m
  use profiling_m
  use varinfo_m
  use mpi_m

  implicit none

#ifdef HAVE_PFFT
  include "pfft.f"
  include "fftw3.f"
#endif

  private ::        &
       decompose,   &
       prime
  public ::         &
    pfft_t,         &
    pfft_all_init,  &
    pfft_all_end,   &
    pfft_init,      &
    pfft_forward_3d,&
    pfft_backward_3d

  ! fftw constants. this is just a copy from file fftw3.f,
  ! distributed with fftw package. 
  
  integer, public, parameter ::          &
       pfftw_forward            =    -1, &
       pfftw_backward           =     1, &
       pfft_measure             =     0 

  type pfft_t
    integer     :: slot         !< in which slot do we have this fft
    integer     :: n(3)         !< size of the fft
    integer     :: is_real      !< is the fft real or complex. PFFT only works with complex
    integer(ptrdiff_t_kind)    :: comm_cart_2d !< 2 dimensional catersian processor grid
    type(c_ptr) :: planf        !< the plan for forward transforms
    type(c_ptr) :: planb        !< the plan for backward transforms
  end type pfft_t

  integer      :: pfft_refs(PFFT_MAX)
  type(pfft_t) :: pfft_array(PFFT_MAX)
  logical      :: fft_optimize
  integer      :: pfft_prepare_plan

contains

  ! ---------------------------------------------------------
  !> initialize the table
  subroutine pfft_all_init()
    integer :: ii

#ifdef HAVE_PFFT
    PUSH_SUB(fft_all_init)

    !%Variable PFFTOptimize
    !%Type logical
    !%Default yes
    !%Section Mesh::FFTs
    !%Description
    !% Should <tt>octopus</tt> optimize the FFT dimensions? 
    !% This means that the cubic mesh to which FFTs are applied is not taken to be as small
    !% as possible: some points may be added to each direction in order to get a "good number"
    !% for the performance of the FFT algorithm.
    !% In some cases, namely when using
    !% the split-operator, or Suzuki-Trotter propagators, this option should be turned off.
    !%End
    call parse_logical(datasets_check('PFFTOptimize'), .true., fft_optimize)
    do ii = 1, PFFT_MAX
      pfft_refs(ii) = FFT_NULL
    end do
 
    !%Variable PFFTPreparePlan
    !%Type integer
    !%Default pfft_measure
    !%Section Mesh::FFTs
    !%Description
    !% The FFTs are performed in octopus with the help of the PFFT and FFTW package (http://www.fftw.org).
    !% Before doing the actual computations, this package prepares a "plan", which means that 
    !% the precise numerical strategy to be followed to compute the FFT is machine/compiler-dependent,
    !% and therefore the software attempts to figure out which is this precise strategy (see the
    !% FFTW documentation for details). This plan preparation, which has to be done for each particular
    !% FFT shape, can be done exhaustively and carefully (slow), or merely estimated. Since this is
    !% a rather critical numerical step, by default it is done carefully, which implies a longer initial
    !% initialization, but faster subsequent computations. You can change this behaviour by changing
    !% this <tt>PFFTPreparePlan</tt> variable, and in this way you can force PFFT to do a fast guess or
    !% estimation of which is the best way to perform the FFT.
    !%Option pfft_measure 0
    !% This is the default, and implies a longer initialization, but involves a more careful analysis
    !% of the strategy to follow, and therefore more efficient FFTs.
    !%Option pfft_estimate 64
    !% This is the "fast initialization" scheme, in which the plan is merely guessed from "reasonable"
    !% assumptions.
    !%End
    call parse_integer(datasets_check('PFFTPreparePlan'), pfft_measure, pfft_prepare_plan)
    if(.not. varinfo_valid_option('PFFTPreparePlan', pfft_prepare_plan)) call input_error('PFFTPreparePlan')

    POP_SUB(fft_all_init)
#endif

  end subroutine pfft_all_init

  ! ---------------------------------------------------------
  !> delete all plans
  subroutine pfft_all_end()
    integer :: ii

#ifdef HAVE_PFFT
    PUSH_SUB(fft_all_end)

    do ii = 1, PFFT_MAX
      if(pfft_refs(ii) /= FFT_NULL) then
        call PDFFT(destroy_plan) (pfft_array(ii)%planf)
        call PDFFT(destroy_plan) (pfft_array(ii)%planb)
        pfft_refs(ii) = FFT_NULL
      end if
    end do

    call PDFFT(cleanup)
    POP_SUB(fft_all_end)
#endif

  end subroutine pfft_all_end

  ! ---------------------------------------------------------
  subroutine pfft_init(nn, dim, is_real, pfft, optimize)
    integer,           intent(inout) :: nn(1:3)
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: is_real
    type(pfft_t),      intent(out)   :: pfft
    logical, optional, intent(in)    :: optimize
   
#ifdef HAVE_PFFT 
    FLOAT, allocatable :: rin(:, :, :)
    CMPLX, allocatable :: cin(:, :, :), cout(:, :, :)
    
    integer :: ii, jj, fft_dim, idir, ierror, process_column_size, process_row_size
    logical :: optimize_
    character(len=100) :: str_tmp
    
    integer(ptrdiff_t_kind) :: alloc_local
    integer(ptrdiff_t_kind) :: local_ni(3), local_i_start(3)
    integer(ptrdiff_t_kind) :: local_no(3), local_o_start(3)
    !    integer comm_cart_2d
    !integer(8) fftplan_forw, fftplan_back
    CMPLX, allocatable ::  data_in(:)
    CMPLX, allocatable ::  data_out(:)
    
    PUSH_SUB(pfft_init)
    ! First, figure out the dimensionality of the FFT. PFFT works efficiently with 3D
    fft_dim = 0
    do ii = 1, dim
      if(nn(ii) <= 1) exit
      fft_dim = fft_dim + 1
    end do

    if(fft_dim < 3) then
      message(1) = "Internal error in pfft_init: apparently, a 1x1x1 FFT is required."
      call write_fatal(1)
    end if

    if(fft_dim > 3) then
      message(1) = "FFT for dimension greater than 3 not implemented."
      call write_fatal(1)
    end if

    optimize_ = .true.
    if(present(optimize)) optimize_ = optimize


    ! OLD: I leave it here because maybe I revert to this method later
    ! optimize dimensions in non-periodic directions
    !    do i = sb%periodic_dim + 1, sb%dim
    !      if(nn(i) /= 1 .and. fft_optimize) &
    !           call loct_fft_optimize(nn(i), 7, 1) ! always ask for an odd number
    !    end do
    ! NEW
    ! optimize dimensions only for finite sys
    optimize_ = optimize_ .and. fft_optimize
    if(optimize_) then
      do ii = 1, fft_dim
        call loct_fft_optimize(nn(ii), 7, 1)
      end do
    end if

    ! find out if fft has already been allocated
    jj = 0
    do ii = PFFT_MAX, 1, -1
      if(pfft_refs(ii) /= FFT_NULL) then
        if(all(nn(1:dim) == pfft_array(ii)%n(1:dim)) .and. is_real == pfft_array(ii)%is_real) then
          pfft = pfft_array(ii)              ! return a copy
          pfft_refs(ii) = pfft_refs(ii) + 1  ! increment the ref count
          POP_SUB(pfft_init)
          return
        end if
      else
        jj = ii
      end if
    end do

    if(jj == 0) then
      message(1) = "Not enough slots for FFTs."
      message(2) = "Please increase PFFT_MAX in pfft.F90 and recompile."
      call write_fatal(2)
    end if

    ! jj now contains an empty slot
    pfft_refs(jj)          = 1
    pfft_array(jj)%slot    = jj
    pfft_array(jj)%n(1:3)  = nn(1:3)
    pfft_array(jj)%is_real = is_real
    
    ! With PFFT only could be 3D complex

    SAFE_ALLOCATE( cin(1:nn(1), 1:nn(2), 1:nn(3)))
    SAFE_ALLOCATE(cout(1:nn(1), 1:nn(2), 1:nn(3)))

    ! Create two-dimensional process grid of
    call decompose(mpi_world%size, process_column_size, process_row_size)
  
    call dpfft_init()

    call PDFFT(create_procmesh_2d) (ierror,MPI_COMM_WORLD,process_column_size,&
         process_row_size,pfft%comm_cart_2d)
    if (ierror .ne. 0) then
      message(1) = "The number of rows and columns in PFFT processor grid is not equal to "
      message(2) = "the number of processor in the MPI communicator."
      message(3) = "Please check it."
      call write_fatal(3)
    end if

    call PDFFT(local_size_3d) (alloc_local, nn, pfft%comm_cart_2d, PFFT_TRANSPOSED, &
         local_ni, local_i_start, local_no, local_o_start)
    
    !     Allocate memory
    SAFE_ALLOCATE(data_in(alloc_local))
    SAFE_ALLOCATE(data_out(alloc_local))

    ! Create the plan, with the processor grid 
    call PDFFT(plan_dft_3d) (pfft_array(jj)%planf, nn, data_in, data_out, &
         pfft%comm_cart_2d, &
         FFTW_FORWARD, PFFT_TRANSPOSED + PFFT_FORWARD, FFTW_MEASURE)
    call PDFFT(plan_dft_3d) (pfft_array(jj)%planb, nn, data_out, data_in, pfft%comm_cart_2d, &
         & FFTW_BACKWARD, PFFT_TRANSPOSED + PFFT_BACKWARD, FFTW_MEASURE) 
    SAFE_DEALLOCATE_A(data_in)
    SAFE_DEALLOCATE_A(data_out)
    
    pfft = pfft_array(jj)
    
    write(message(1), '(a)') "Info: PFFT allocated with size ("
    do idir = 1, dim
      write(str_tmp, '(i7,a)') nn(idir)
      if(idir == dim) then
        message(1) = trim(message(1)) // trim(str_tmp) // ") in slot "
      else
        message(1) = trim(message(1)) // trim(str_tmp) // ","
      endif
    enddo
    write(str_tmp, '(i2)') jj
    message(1) = trim(message(1)) // trim(str_tmp)
    call write_info(1)
    
    POP_SUB(pfft_init)
#endif
        
  end subroutine pfft_init

  ! ---------------------------------------------------------
  !> These routines simply call pfft
  !! first the complex to complex versions
  subroutine pfft_forward_3d(pfft)
    type(pfft_t), intent(in)  :: pfft

#ifdef HAVE_PFFT
    PUSH_SUB(pfft_forward_3d)

    call PDFFT(execute) (pfft%planf)

    POP_SUB(pfft_forward_3d)
#endif

  end subroutine pfft_forward_3d 

  subroutine pfft_backward_3d(pfft)
    type(pfft_t), intent(in) :: pfft
    
#ifdef HAVE_PFFT
    PUSH_SUB(dpfft_backward)
    
    call PDFFT(execute) (pfft%planb)

    ! multiply by 1/(N1*N2*N2)
    !call lalg_scal(fft%n(1), fft%n(2), fft%n(3), &
    !  M_ONE / (fft%n(1)*fft%n(2)*fft%mn(3)), rr)

    POP_SUB(pfft_backward)
#endif

  end subroutine pfft_backward_3d
  
  !> This function decomposes a given number of processors in a
  !! two dimension processor grid.
  !! @author Miquel Huix 
  subroutine decompose(n_proc, dim1, dim2)
    integer, intent(in)  :: n_proc !< Number of processors
    integer, intent(out) :: dim1   !< First out dimension
    integer, intent(out) :: dim2   !< Second out dimension
    integer :: np, i
    
    if(n_proc < 1) then
      message(1) = "Internal error in pfft_init: "
      message(2) = "Error creating the decomposition of the processor grid"  
      message(3) = "The number of processors could not be negative"
      call write_fatal(3)
    end if
    
    dim1 = 1
    dim2 = 1
    np = n_proc
    i = n_proc-1
    
    if (prime(n_proc)) then
      dim1 = n_proc
    else
      do
        if (i <= 1) exit
        if(mod(np,i) /= 0.or.(.not.prime(i))) then
          i=i-1
          cycle
        end if
        np=np/i
        if (dim1 <= dim2) then
          dim1 = dim1*i
        else
          dim2 = dim2*i
        end if
      end do
    end if
    
    if (dim1 * dim2 /= n_proc) then
      message(1) = "Internal error in pfft_init: "
      message(2) = "Error creating the decomposition of the processor grid"  
      message(3) = "The multiplication of the two dimensions have to be equal to the number of processors"
      call write_fatal(3)
    end if
    
    write(message(1),'(a)') "Info: PFFT processor grid"
    write(message(2),'(a, i9)') " No. of processors                = ",n_proc
    write(message(3),'(a, i9)') " No. of columns in the proc. grid = ",dim1
    write(message(4),'(a, i9)') " No. of rows    in the proc. grid = ",dim2
    call write_info(4)

  end subroutine decompose

  logical function prime(n) result(is_prime)
    integer, intent(in) :: n
    
    integer :: i, root

    if (n < 1) then
      message(1) = "Internal error in pfft_init: "
      message(2) = "Error calculating the negative prime number"
      call write_fatal(2)
    end if
    if (n == 1) then
      is_prime = .false.
    else
      root = sqrt(real(n))
      do i = 2, root
        if (mod(n,i) == 0) then
          is_prime = .false.
          return
        end if
      end do
      is_prime = .true.
    end if
  end function prime

  
end module pfft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
