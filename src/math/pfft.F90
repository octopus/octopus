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
#define PFFT_NULL -1

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
  use mpi_m
  use parser_m
  use profiling_m
  use varinfo_m

  implicit none

#ifdef HAVE_PFFT
  include "pfft.f"
  include "fftw3.f"

  private ::        &
       decompose,   &
       prime
  public ::         &
    pfft_t,         &
    pfft_all_init,  &
    pfft_all_end,   &
    pfft_end,       &
    pfft_init,      &
    pfft_forward_3d, &
    pfft_backward_3d

  ! fftw constants. this is just a copy from file fftw3.f,
  ! distributed with fftw package. 
  
  integer, public, parameter ::          &
       pfftw_forward            =    -1, &
       pfftw_backward           =     1, &
       pfft_measure             =     0 
  
  type pfft_t
    integer                 :: slot         !< in which slot do we have this fft
    integer(ptrdiff_t_kind) :: n(3)         !< size of the fft
    integer                 :: is_real      !< is the fft real or complex. PFFT only works with complex. Real = 0, Complex = 1
    integer(ptrdiff_t_kind) :: comm_cart_2d !< 2-dimensional Cartesian processor grid
    integer(ptrdiff_t_kind) :: planf        !< the plan for forward transforms
    integer(ptrdiff_t_kind) :: planb        !< the plan for backward transforms
    integer(ptrdiff_t_kind) :: alloc_local

    integer, allocatable :: begin_indexes(:) !< where does each process start
    integer, allocatable :: block_sizes(:)   !< size of the block that is going to be used in the gatherv

    integer(ptrdiff_t_kind) :: local_ni(3)      !< local input points
    integer(ptrdiff_t_kind) :: local_i_start(3) !< local input start index
    integer(ptrdiff_t_kind) :: local_no(3)      !< local output points
    integer(ptrdiff_t_kind) :: local_o_start(3) !< local output start index

    CMPLX, pointer :: data_in(:)                   !< input data
    CMPLX, pointer :: data_out(:)                  !< output data
    FLOAT, pointer :: global_data_in(:)            !< ALL the input data
    integer, allocatable :: get_local_index(:,:,:) !< Mapping vector between local and global pfft indexing
  end type pfft_t

  integer      :: pfft_refs(PFFT_MAX)
  type(pfft_t) :: pfft_array(PFFT_MAX)
  logical      :: fft_optimize
  logical      :: pfft_output
  integer      :: pfft_prepare_plan

contains

  ! ---------------------------------------------------------
  !> initialize the table
  subroutine pfft_all_init()
    integer :: ii

    PUSH_SUB(pfft_all_init)

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
      pfft_refs(ii) = PFFT_NULL
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
    
    !%Variable PFFTOutput
    !%Type logical
    !%Default no
    !%Section Mesh::FFTs
    !%Description
    !% Should <tt>octopus</tt> print the FFT values? 
    !%End
    call parse_logical(datasets_check('PFFTOutput'), .false., pfft_output)
    
    do ii = 1, PFFT_MAX
      pfft_refs(ii) = PFFT_NULL
    end do
    POP_SUB(pfft_all_init)
  end subroutine pfft_all_init

  ! ---------------------------------------------------------
  !> delete all plans
  subroutine pfft_all_end()
    integer :: ii

    PUSH_SUB(pfft_all_end)

    do ii = 1, PFFT_MAX
      if(pfft_refs(ii) /= PFFT_NULL) then
        call PDFFT(destroy_plan) (pfft_array(ii)%planf)
        call PDFFT(destroy_plan) (pfft_array(ii)%planb)
        SAFE_DEALLOCATE_P(pfft_array(ii)%data_out)
        SAFE_DEALLOCATE_P(pfft_array(ii)%data_in)
        SAFE_DEALLOCATE_P(pfft_array(ii)%global_data_in)
        SAFE_DEALLOCATE_A(pfft_array(ii)%begin_indexes)
        SAFE_DEALLOCATE_A(pfft_array(ii)%block_sizes)
        SAFE_DEALLOCATE_A(pfft_array(ii)%get_local_index)
        pfft_refs(ii) = PFFT_NULL
      end if
    end do

    call PDFFT(cleanup)
    POP_SUB(pfft_all_end)
  end subroutine pfft_all_end

  ! ---------------------------------------------------------
  subroutine pfft_end(pfft)
    type(pfft_t), intent(inout) :: pfft

    integer :: ii

    PUSH_SUB(pfft_end)
    
    ii = pfft%slot
    if(pfft_refs(ii) == PFFT_NULL) then
      message(1) = "Trying to deallocate PFFT that has not been allocated."
      call messages_warning(1)
    else
      if(pfft_refs(ii) > 1) then
        pfft_refs(ii) = pfft_refs(ii) - 1
      else
        pfft_refs(ii) = PFFT_NULL
        call PDFFT(destroy_plan) (pfft_array(ii)%planf)
        call PDFFT(destroy_plan) (pfft_array(ii)%planb)
        SAFE_DEALLOCATE_P(pfft_array(ii)%data_out)
        SAFE_DEALLOCATE_P(pfft_array(ii)%data_in)
        SAFE_DEALLOCATE_P(pfft_array(ii)%global_data_in)
        SAFE_DEALLOCATE_A(pfft_array(ii)%begin_indexes)
        SAFE_DEALLOCATE_A(pfft_array(ii)%block_sizes)
        SAFE_DEALLOCATE_A(pfft_array(ii)%get_local_index)
        write(message(1), '(a,i4)') "Info: PFFT deallocated from slot ", ii
        call messages_info(1)
      end if
    end if

    POP_SUB(pfft_end)
  end subroutine pfft_end

  ! ---------------------------------------------------------
  subroutine pfft_init(nn, dim, is_real, pfft, optimize)
    integer,           intent(inout) :: nn(1:3)
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: is_real
    type(pfft_t),      intent(out)   :: pfft
    logical, optional, intent(in)    :: optimize
   
    FLOAT, allocatable :: rin(:, :, :)
    CMPLX, allocatable :: cin(:, :, :), cout(:, :, :)
    
    integer :: ii, jj, fft_dim, idir, ierror, process_column_size, process_row_size
    logical :: optimize_
    character(len=100) :: str_tmp
  
    integer(ptrdiff_t_kind) :: nn_t_kind(1:3)
    
    PUSH_SUB(pfft_init)

    ! First, figure out the dimensionality of the FFT. PFFT works efficiently with 3D
    fft_dim = 0
    do ii = 1, dim
      if(nn(ii) <= 1) exit
      fft_dim = fft_dim + 1
    end do

    if(fft_dim < 3) then
      message(1) = "Internal error in pfft_init: apparently, a 1x1x1 FFT is required."

      call messages_fatal(1)
    end if

    if(fft_dim > 3) then
      message(1) = "FFT for dimension greater than 3 not implemented."
      call messages_fatal(1)
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
      if(pfft_refs(ii) /= PFFT_NULL) then
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
      call messages_fatal(2)
    end if

    ! jj now contains an empty slot
    pfft_refs(jj)          = 1
    pfft_array(jj)%slot    = jj

    !! Change the nn from 32 bits to 64 (if needed)  
    pfft_array(jj)%n(1:3)  = nn(1:3)
    pfft_array(jj)%is_real = is_real
    
    ! With PFFT only could be 3D complex
    SAFE_ALLOCATE( cin(1:pfft_array(jj)%n(1), 1:pfft_array(jj)%n(2), 1:pfft_array(jj)%n(3)))
    SAFE_ALLOCATE(cout(1:pfft_array(jj)%n(1), 1:pfft_array(jj)%n(2), 1:pfft_array(jj)%n(3)))
    
    ! Create two-dimensional process grid of
    call decompose(mpi_world%size, process_column_size, process_row_size)
  
    call dpfft_init()

    call PDFFT(create_procmesh_2d) (ierror,MPI_COMM_WORLD,process_column_size,&
         process_row_size,pfft_array(jj)%comm_cart_2d)
    if (ierror .ne. 0) then
      message(1) = "The number of rows and columns in PFFT processor grid is not equal to "
      message(2) = "the number of processor in the MPI communicator."
      message(3) = "Please check it."
      call messages_fatal(3)
    end if
    
    call PDFFT(local_size_3d) (pfft_array(jj)%alloc_local, pfft_array(jj)%n, pfft_array(jj)%comm_cart_2d, PFFT_TRANSPOSED, &
         pfft_array(jj)%local_ni, pfft_array(jj)%local_i_start, pfft_array(jj)%local_no, pfft_array(jj)%local_o_start)

    !     Allocate memory
    SAFE_ALLOCATE(pfft_array(jj)%data_in(pfft_array(jj)%alloc_local))
    SAFE_ALLOCATE(pfft_array(jj)%data_out(pfft_array(jj)%alloc_local))
    
    ! Collect the indexes of all processes
    call pfft_do_mapping(pfft_array(jj))

    ! Create the plan, with the processor grid 
    call PDFFT(plan_dft_3d) (pfft_array(jj)%planf, pfft_array(jj)%n, & 
         pfft_array(jj)%data_in, pfft_array(jj)%data_out, pfft_array(jj)%comm_cart_2d, &
         FFTW_FORWARD, PFFT_TRANSPOSED + PFFT_FORWARD, FFTW_MEASURE)
    call PDFFT(plan_dft_3d) (pfft_array(jj)%planb, pfft_array(jj)%n, &
         pfft_array(jj)%data_out, pfft_array(jj)%data_in, pfft_array(jj)%comm_cart_2d, &
         FFTW_BACKWARD, PFFT_TRANSPOSED + PFFT_BACKWARD, FFTW_MEASURE) 
       
    write(message(1), '(a)') "Info: PFFT allocated with size ("
    do idir = 1, dim
      write(str_tmp, '(i7,a)') pfft_array(jj)%n(idir)
      if(idir == dim) then
        message(1) = trim(message(1)) // trim(str_tmp) // ") in slot "
      else
        message(1) = trim(message(1)) // trim(str_tmp) // ","
      endif
    enddo
    write(str_tmp, '(i2)') jj
    message(1) = trim(message(1)) // trim(str_tmp)
    call messages_info(1)

    pfft = pfft_array(jj)

    POP_SUB(pfft_init)       
  end subroutine pfft_init

  ! ---------------------------------------------------------
  !> These routines simply call pfft
  !! first the complex to complex versions
  subroutine pfft_forward_3d(pfft)
    type(pfft_t), intent(inout)  :: pfft

    type(profile_t), save :: prof_fw

    PUSH_SUB(pfft_forward_3d)
    
    call profiling_in(prof_fw,"PFFT_FW")
    call PDFFT(execute) (pfft%planf)
    call profiling_out(prof_fw)
    
    POP_SUB(pfft_forward_3d)
  end subroutine pfft_forward_3d 

  ! ---------------------------------------------------------
  subroutine pfft_backward_3d(pfft)
    type(pfft_t), intent(inout) :: pfft
    integer :: index, ii, jj, kk
    type(profile_t), save :: prof_bw,prof_gtv,prof_g,prof_t
    character(len=256) ::  tmp_file !< kentzeko!!!

    PUSH_SUB(pfft_backward_3d)
    
    call profiling_in(prof_bw,"PFFT_BW")
    call PDFFT(execute) (pfft%planb)
    call profiling_out(prof_bw)
    
    call profiling_in(prof_t,"PFFT_TRANS")
    !aling the data
    index = 1
    do kk = pfft%local_i_start(3), pfft%local_i_start(3)+pfft%local_ni(3)-1
      do jj = pfft%local_i_start(2), pfft%local_i_start(2)+pfft%local_ni(2)-1
        do ii = pfft%local_i_start(1), pfft%local_i_start(1)+pfft%local_ni(1)-1
          pfft%global_data_in(pfft%get_local_index(ii,jj,kk)) = real(pfft%data_in(index))
          index = index + 1 
        end do
      end do
    end do
    close(13)
    call profiling_out(prof_t)

    !collect the data in all processes
    call profiling_in(prof_g,"PFFT_GATV")
    call MPI_Allgatherv ( &
         pfft%global_data_in(pfft%begin_indexes(mpi_world%rank+1)), &
         pfft%block_sizes(mpi_world%rank+1), MPI_FLOAT, &
         pfft%global_data_in(1), &
         pfft%block_sizes(1),pfft%begin_indexes - 1, &
         MPI_FLOAT, &
         mpi_world%comm, mpi_err )
    if (mpi_err /= 0) then
      write(message(1),'(a)')"MPI_Allgatherv failed in pfft.F90"
      call messages_fatal(1)
    end if
    call profiling_out(prof_g) 

    POP_SUB(pfft_backward_3d)
  end subroutine pfft_backward_3d
  
  ! ---------------------------------------------------------
  !> do the mapping between global and local points of PFFT arrays
  subroutine pfft_do_mapping(pfft)
    type(pfft_t), intent(inout) :: pfft
    
    integer :: tmp_local(6), position, process,ii, jj, kk, index
    integer, allocatable ::local_sizes(:)
    type(profile_t), save ::  prof_gt, prof_a
    call profiling_in(prof_gt,"PFFT_GAT")
    !!BEGIN:gather the local information into a unique vector.
    !!do a gather in 3d of all the box, into a loop
    tmp_local(1) = pfft%local_i_start(1)
    tmp_local(2) = pfft%local_i_start(2) 
    tmp_local(3) = pfft%local_i_start(3) 
    tmp_local(4) = pfft%local_ni(1)
    tmp_local(5) = pfft%local_ni(2)
    tmp_local(6) = pfft%local_ni(3) 
    SAFE_ALLOCATE(local_sizes(6*mpi_world%size))
    call MPI_Allgather(tmp_local,6,MPI_INTEGER, &
         local_sizes,6,MPI_INTEGER,&
         mpi_world%comm,mpi_err)
    call profiling_out(prof_gt)
    
    call profiling_in(prof_a,"PFFT_ALLOC")
    SAFE_ALLOCATE(pfft%begin_indexes(mpi_world%size))
    SAFE_ALLOCATE(pfft%block_sizes(mpi_world%size))
    SAFE_ALLOCATE(pfft%global_data_in(pfft%n(1)*pfft%n(2)*pfft%n(3)))
    SAFE_ALLOCATE(pfft%get_local_index(pfft%n(1),pfft%n(2),pfft%n(3)))

    do process=1,mpi_world%size
      position = ((process-1)*6)+1
      if (position == 1) then
        pfft%begin_indexes(1) = 1
        pfft%block_sizes(1)  = local_sizes(4)*local_sizes(5)*local_sizes(6)
      else
        ! calculate the begin index and size of each process
        pfft%begin_indexes(process) =  pfft%begin_indexes(process-1) +  pfft%block_sizes(process-1)
        pfft%block_sizes(process) = local_sizes(position+3)*local_sizes(position+4)*local_sizes(position+5)
      end if

      !save the mapping between the global x,y,z and the local index
      index = 0
      do kk = local_sizes(position+2), local_sizes(position+2)+local_sizes(position+5)-1
        do jj = local_sizes(position+1), local_sizes(position+1)+local_sizes(position+4)-1
          do ii = local_sizes(position), local_sizes(position)+local_sizes(position+3)-1
            pfft%get_local_index(ii,jj,kk) = index + pfft%begin_indexes(process)
            index = index + 1
          end do
        end do
      end do
    end do    

    call profiling_out(prof_a)
    
  end subroutine pfft_do_mapping
  ! ---------------------------------------------------------
  !> This function decomposes a given number of processors into a
  !! two-dimensional processor grid.
  !! @author Miquel Huix 
  subroutine decompose(n_proc, dim1, dim2)
    integer, intent(in)  :: n_proc !< Number of processors
    integer, intent(out) :: dim1   !< First out dimension
    integer, intent(out) :: dim2   !< Second out dimension
    integer :: np, i
    
    PUSH_SUB(decompose)

    if(n_proc < 1) then
      message(1) = "Internal error in pfft_init: "
      message(2) = "Error creating the decomposition of the processor grid"  
      message(3) = "The number of processors could not be negative"
      call messages_fatal(3)
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
      call messages_fatal(3)
    end if
    
    write(message(1),'(a)') "Info: PFFT processor grid"
    write(message(2),'(a, i9)') " No. of processors                = ",n_proc
    write(message(3),'(a, i9)') " No. of columns in the proc. grid = ",dim1
    write(message(4),'(a, i9)') " No. of rows    in the proc. grid = ",dim2
    write(message(5),'(a, i9)') " The size of integer is = ",ptrdiff_t_kind
    call messages_info(5)

    POP_SUB(decompose)
  end subroutine decompose

  ! ---------------------------------------------------------
  logical function prime(n) result(is_prime)
    integer, intent(in) :: n
    
    integer :: i, root

    PUSH_SUB(prime)

    if (n < 1) then
      message(1) = "Internal error in pfft_init: "
      message(2) = "Error calculating the negative prime number"
      call messages_fatal(2)
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

    POP_SUB(prime)
  end function prime
#endif
  
end module pfft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
