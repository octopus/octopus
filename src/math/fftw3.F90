!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id$

#if defined(SINGLE_PRECISION)
#  define DFFTW(x) sfftw_ ## x
#else
#  define DFFTW(x) dfftw_ ## x
#endif

module fft_m
  use global_m
  use messages_m
  use varinfo_m
  use datasets_m
  use loct_math_m
  use parser_m
  use lalg_basic_m
  use c_pointer_m
  use profiling_m

  implicit none

  private
  public ::        &
    fft_t,         &
    fft_all_init,  &
    fft_all_end,   &
    fft_init,      &
    fft_end,       &
    fft_copy,      &
    pad_feq,       &
    dfft_forward,  &
    zfft_forward,  &
    dfft_backward, &
    zfft_backward, &
    dfft_forward1, &
    dfft_backward1

  ! global constants
  integer, public, parameter ::                &
    fft_real    = 0,                   &
    fft_complex = 1

  ! fftw constants. this is just a copy from file fftw3.f,
  ! distributed with fftw package.
  integer, public, parameter ::                &
    fftw_r2hc                =      0, &
    fftw_hc2r                =      1, &
    fftw_dht                 =      2, &
    fftw_redft00             =      3, &
    fftw_redft01             =      4, &
    fftw_redft10             =      5, &
    fftw_redft11             =      6, &
    fftw_rodft00             =      7, &
    fftw_rodft01             =      8, &
    fftw_rodft10             =      9, &
    fftw_rodft11             =     10, &
    fftw_forward             =     -1, &
    fftw_backward            =      1, &
    fftw_measure             =      0, &
    fftw_destroy_input       =      1, &
    fftw_unaligned           =      2, &
    fftw_conserve_memory     =      4, &
    fftw_exhaustive          =      8, &
    fftw_preserve_input      =     16, &
    fftw_patient             =     32, &
    fftw_estimate            =     64, &
    fftw_estimate_patient    =    128, &
    fftw_believe_pcost       =    256, &
    fftw_dft_r2hc_icky       =    512, &
    fftw_nonthreaded_icky    =   1024, &
    fftw_no_buffering        =   2048, &
    fftw_no_indirect_op      =   4096, &
    fftw_allow_large_generic =   8192, &
    fftw_no_rank_splits      =  16384, &
    fftw_no_vrank_splits     =  32768, &
    fftw_no_vrecurse         =  65536, &
    fftw_no_simd             = 131072

  type fft_t
    integer   :: slot                ! in which slot do we have this fft

    integer   :: n(3)          ! size of the fft
    integer   :: is_real             ! is the fft real or complex
    type(c_ptr) :: planf ! the plan for forward transforms
    type(c_ptr) :: planb ! the plan for backward transforms
  end type fft_t

  integer :: fft_refs(FFT_MAX)
  type(fft_t) :: fft_array(FFT_MAX)
  logical :: fft_optimize
  integer :: fft_prepare_plan

contains

  ! ---------------------------------------------------------
  ! initialize the table
  subroutine fft_all_init()
    integer :: i

    call push_sub('fftw3.fft_all_init')

    !%Variable FFTOptimize
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
    call parse_logical(datasets_check('FFTOptimize'), .true., fft_optimize)
    do i = 1, FFT_MAX
      fft_refs(i) = NULL
    end do

    !%Variable FFTPreparePlan
    !%Type integer
    !%Default 0
    !%Section Mesh::FFTs
    !%Description
    !% The FFTs are performed in octopus with the help of the FFTW package (http://www.fftw.org).
    !% Before doing the actual computations, this package prepares a "plan", which means that 
    !% the precise numerical strategy to be followed to compute the FFT is machine/compiler-dependent,
    !% and therefore the software attempts to figure out which is this precise strategy (see the
    !% FFTW documentation for details). This plan preparation, which has to be done for each particular
    !% FFT shape, can be done exhaustively and carefully (slow), or merely estimated. Since this is
    !% a rather critical numerical step, by default it is done carefully, which implies a longer initial
    !% initialization, but faster subsequent computations. You can change this behaviour by changing
    !% this <tt>FFTPreparePlan</tt> variable, and in this way you can force FFTW to do a fast guess or
    !% estimation of which is the best way to perform the FFT.
    !%Option fftw_measure 0
    !% This is the default, and implies a longer initialization, but involves a more careful analysis
    !% of the strategy to follow, and therefore more efficient FFTs.
    !%Option fftw_estimate 64
    !% This is the "fast initialization" scheme, in which the plan is merely guessed from "reasonable"
    !% assumptions.
    !%End
    call parse_integer(datasets_check('FFTPreparePlan'), fftw_measure, fft_prepare_plan)
    if(.not.varinfo_valid_option('FFTPreparePlan', fft_prepare_plan)) call input_error('FFTPreparePlan')

    call pop_sub()
  end subroutine fft_all_init


  ! ---------------------------------------------------------
  ! delete all plans
  subroutine fft_all_end()
    integer :: i

    call push_sub('fftw3.fft_all_end')

    do i = 1, FFT_MAX
      if(fft_refs(i).ne.NULL) then
        call DFFTW(destroy_plan) (fft_array(i)%planf)
        call DFFTW(destroy_plan) (fft_array(i)%planb)
        fft_refs(i) = NULL
      end if
    end do

    call DFFTW(cleanup)
    call pop_sub()
  end subroutine fft_all_end


  ! ---------------------------------------------------------
  subroutine fft_init(n, is_real, fft, optimize)
    integer,           intent(inout) :: n(3)
    integer,           intent(in)    :: is_real
    type(fft_t),       intent(out)   :: fft
    logical, optional, intent(in)    :: optimize

    FLOAT, allocatable :: rin(:, :, :)
    CMPLX, allocatable :: cin(:, :, :), cout(:, :, :)

    integer :: i, j, dim
    logical :: optimize_

    call push_sub('fftw3.fft_init')

    ! First, figure out the dimensionality of the FFT.
    dim = 0
    do i = 1, 3
      if(n(i) <= 1) exit
      dim = dim + 1
    end do

    if(dim.eq.0) then
      message(1) = "Internal error in fft_init: apparently, a 1x1x1 FFT is required."
      call write_fatal(1)
    end if

    optimize_ = .true.
    if(present(optimize)) optimize_ = optimize


    ! OLD: I leave it here because maybe I revert to this method later
    ! optimize dimensions in non-periodic directions
    !    do i = sb%periodic_dim + 1, sb%dim
    !      if(n(i) /= 1 .and. fft_optimize) &
    !           call loct_fft_optimize(n(i), 7, 1) ! always ask for an odd number
    !    end do
    ! NEW
    ! optimize dimensions only for finite sys
    optimize_ = optimize_ .and. fft_optimize
    if(optimize_) then
      do i = 1, dim
        call loct_fft_optimize(n(i), 7, 1)
      end do
    end if

    ! find out if fft has already been allocated
    j = 0
    do i = FFT_MAX, 1, -1
      if(fft_refs(i).ne.NULL) then
        if(n(1) == fft_array(i)%n(1).and.n(2) == fft_array(i)%n(2).and.n(3) == fft_array(i)%n(3).and. &
          is_real == fft_array(i)%is_real) then
          fft = fft_array(i)             ! return a copy
          fft_refs(i) = fft_refs(i) + 1  ! increment the ref count
          return
        end if
      else
        j = i
      end if
    end do

    if(j == 0) then
      message(1) = "Not enough slots for FFTs"
      message(2) = "Please increase FFT_MAX in fft.F90 and recompile"
      call write_fatal(2)
    end if

    ! j now contains an empty slot
    fft_refs(j)          = 1
    fft_array(j)%slot    = j
    fft_array(j)%n       = n
    fft_array(j)%is_real = is_real
    if(is_real == fft_real) then
      SAFE_ALLOCATE(rin(1:n(1), 1:n(2), 1:n(3)))
      SAFE_ALLOCATE(cout(1:n(1)/2+1, 1:n(2), 1:n(3)))

      select case(dim)
      case(3)
        call DFFTW(plan_dft_r2c_3d) (fft_array(j)%planf, n(1), n(2), n(3), rin, cout, fft_prepare_plan+fftw_unaligned)
        call DFFTW(plan_dft_c2r_3d) (fft_array(j)%planb, n(1), n(2), n(3), cout, rin, fft_prepare_plan+fftw_unaligned)
      case(2)
        call DFFTW(plan_dft_r2c_2d) (fft_array(j)%planf, n(1), n(2), rin, cout, fft_prepare_plan+fftw_unaligned)
        call DFFTW(plan_dft_c2r_2d) (fft_array(j)%planb, n(1), n(2), cout, rin, fft_prepare_plan+fftw_unaligned)
      case(1)
        call DFFTW(plan_dft_r2c_1d) (fft_array(j)%planf, n(1), rin, cout, fft_prepare_plan+fftw_unaligned)
        call DFFTW(plan_dft_c2r_1d) (fft_array(j)%planb, n(1), cout, rin, fft_prepare_plan+fftw_unaligned)
      end select
      SAFE_DEALLOCATE_A(rin)
      SAFE_DEALLOCATE_A(cout)
    else
      SAFE_ALLOCATE( cin(1:n(1), 1:n(2), 1:n(3)))
      SAFE_ALLOCATE(cout(1:n(1), 1:n(2), 1:n(3)))
      select case(dim)
      case(3)
        call DFFTW(plan_dft_3d) (fft_array(j)%planf, n(1), n(2), n(3), cin, cout, fftw_forward,  fft_prepare_plan)
        call DFFTW(plan_dft_3d) (fft_array(j)%planb, n(1), n(2), n(3), cin, cout, fftw_backward, fft_prepare_plan)
      case(2)
        call DFFTW(plan_dft_2d) (fft_array(j)%planf, n(1), n(2), cin, cout, fftw_forward,  fft_prepare_plan)
        call DFFTW(plan_dft_2d) (fft_array(j)%planb, n(1), n(2), cin, cout, fftw_backward, fft_prepare_plan)
      case(1)
        call DFFTW(plan_dft_1d) (fft_array(j)%planf, n(1), cin, cout, fftw_forward,  fft_prepare_plan)
        call DFFTW(plan_dft_1d) (fft_array(j)%planb, n(1), cin, cout, fftw_backward, fft_prepare_plan)
      end select
      SAFE_DEALLOCATE_A(cin)
      SAFE_DEALLOCATE_A(cout)
    end if
    fft = fft_array(j)

    write(message(1), '(a,i7,a,i7,a,i7,a,i2)') "Info: FFT allocated with size (", &
      n(1), ",", n(2), ",", n(3), ") in slot ", j
    call write_info(1)

    call pop_sub()
  end subroutine fft_init


  ! ---------------------------------------------------------
  subroutine fft_copy(fft_i, fft_o)
    type(fft_t), intent( in) :: fft_i
    type(fft_t), intent(out) :: fft_o

    call push_sub('fftw3.fft_copy')

    ASSERT(fft_i%slot>=1.and.fft_i%slot<=FFT_MAX)
    ASSERT(fft_refs(fft_i%slot) > 0)

    fft_o = fft_i
    fft_refs(fft_i%slot) = fft_refs(fft_i%slot) + 1

    call pop_sub()
  end subroutine fft_copy


  ! ---------------------------------------------------------
  subroutine fft_end(fft)
    type(fft_t), intent(inout) :: fft

    integer :: i

    call push_sub('fftw3.fft_end')

    i = fft%slot
    if(fft_refs(i)==NULL) then
      message(1) = "Trying to deallocate FFT that has not been allocated"
      call write_warning(1)
    else
      if(fft_refs(i) > 1) then
        fft_refs(i) = fft_refs(i) - 1
      else
        fft_refs(i) = NULL
        call DFFTW(destroy_plan) (fft_array(i)%planf)
        call DFFTW(destroy_plan) (fft_array(i)%planb)
        write(message(1), '(a,i4,a,i4,a,i4,a,i2)') "Info: FFT deallocated from slot ", i
        call write_info(1)
      end if
    end if

    call pop_sub()
  end subroutine fft_end


  ! ---------------------------------------------------------
  subroutine fft_getdim_real(fft, d)
    type(fft_t), intent(in) :: fft
    integer,    intent(out) :: d(3)

    call push_sub('fftw3.fft_getdim_real')
    d = fft%n

    call pop_sub()
  end subroutine fft_getdim_real


  ! ---------------------------------------------------------
  subroutine fft_getdim_complex(fft, d)
    type(fft_t), intent(in) :: fft
    integer,    intent(out) :: d(3)

    call push_sub('fftw3.fft_getdim_complex')

    d = fft%n
    if(fft%is_real == fft_real)  d(1) = d(1)/2 + 1

    call pop_sub()
  end subroutine fft_getdim_complex


  ! ---------------------------------------------------------
  ! these routines simply call fftw
  ! first the real to complex versions
  subroutine dfft_forward(fft, r, c)
    type(fft_t), intent(in)  :: fft
    FLOAT,       intent(in)  :: r(:,:,:)
    CMPLX,       intent(out) :: c(:,:,:)

    call push_sub('fftw3.dfft_forward')

    call DFFTW(execute_dft_r2c) (fft%planf, r, c)

    call pop_sub()

  end subroutine dfft_forward


  ! ---------------------------------------------------------
  ! these routines simply call fftw
  ! first the real to complex versions
  subroutine dfft_forward1(fft, r, c)
    type(fft_t), intent(in)  :: fft
    FLOAT,       intent(in)  :: r(:)
    CMPLX,       intent(out) :: c(:)

    call push_sub('fftw3.dfft_forward1')

    call DFFTW(execute_dft_r2c) (fft%planf, r, c)

    call pop_sub()

  end subroutine dfft_forward1


  ! ---------------------------------------------------------
  subroutine dfft_backward(fft, c, r)
    type(fft_t), intent(in) :: fft
    CMPLX, intent(in)  :: c(fft%n(1), fft%n(2), fft%n(3))
    FLOAT, intent(out) :: r(fft%n(1), fft%n(2), fft%n(3))

    call push_sub('fftw3.dfft_backward')
    
    call DFFTW(execute_dft_c2r) (fft%planb, c, r)

    ! multiply by 1/(N1*N2*N2)
    call lalg_scal(fft%n(1), fft%n(2), fft%n(3), &
      M_ONE/real(fft%n(1)*fft%n(2)*fft%n(3), REAL_PRECISION), r)

    call pop_sub()

  end subroutine dfft_backward


  ! ---------------------------------------------------------
  subroutine dfft_backward1(fft, c, r)
    type(fft_t), intent(in) :: fft
    CMPLX, intent(in)  :: c(fft%n(1))
    FLOAT, intent(out) :: r(fft%n(1))

    call push_sub('fftw3.dfft_backward1')
    
    call DFFTW(execute_dft_c2r) (fft%planb, c, r)

    ! multiply by 1/(N1*N2*N2)
    call lalg_scal(fft%n(1), &
      M_ONE/real(fft%n(1), REAL_PRECISION), r)

    call pop_sub()

  end subroutine dfft_backward1


  ! ---------------------------------------------------------
  ! first the complex versions
  subroutine zfft_forward(fft, in, out)
    type(fft_t), intent(in) :: fft
    CMPLX, intent(in)  :: in(:,:,:)
    CMPLX, intent(out) :: out(:,:,:)

    call push_sub('fftw3.zfft_forward')
    
    call DFFTW(execute_dft) (fft%planf, in, out)
    
    call pop_sub()
    
  end subroutine zfft_forward


  ! ---------------------------------------------------------
  subroutine zfft_backward(fft, in, out)
    type(fft_t), intent(in) :: fft
    CMPLX, intent(in)  ::  in(fft%n(1), fft%n(2), fft%n(3))
    CMPLX, intent(out) :: out(fft%n(1), fft%n(2), fft%n(3))

    call push_sub('fftw3.zfft_backward')

    call DFFTW(execute_dft) (fft%planb, in, out)

    ! multiply by 1/(N1*N2*N2)
    call lalg_scal(fft%n(1), fft%n(2), fft%n(3), &
      M_z1/real(fft%n(1)*fft%n(2)*fft%n(3), REAL_PRECISION), out)

    call pop_sub()

  end subroutine zfft_backward


  ! ---------------------------------------------------------
  ! convert between array index and G vector
  function pad_feq(i, n, mode)
    integer, intent(in) :: i,n
    logical, intent(in) :: mode
    integer :: pad_feq

    ! no push_sub: called too frequently

    if(mode) then      ! index to frequency number
      if( i <= n/2 + 1 ) then
        pad_feq = i - 1
      else
        pad_feq = i - n - 1
      end if
    else
      if( i >= 0 ) then
        pad_feq = i + 1
      else
        pad_feq = i + n + 1
      end if
    end if

    return
  end function pad_feq

end module fft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
