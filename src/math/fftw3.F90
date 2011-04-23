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

#include "global.h"

#define FFT_MAX 10
#define FFT_NULL -1

#if defined(SINGLE_PRECISION)
#  define DFFTW(x) sfftw_ ## x
#else
#  define DFFTW(x) dfftw_ ## x
#endif

module fft_m
  use c_pointer_m
  use datasets_m
  use global_m
  use lalg_basic_m
  use loct_math_m
  use messages_m
  use parser_m
  use profiling_m
  use varinfo_m

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
    integer     :: slot       !< in which slot do we have this fft

    integer     :: n(3)       !< size of the fft
    integer     :: is_real    !< is the fft real or complex
    type(c_ptr) :: planf      !< the plan for forward transforms
    type(c_ptr) :: planb      !< the plan for backward transforms
  end type fft_t

  integer     :: fft_refs(FFT_MAX)
  type(fft_t) :: fft_array(FFT_MAX)
  logical     :: fft_optimize
  integer     :: fft_prepare_plan

contains

  ! ---------------------------------------------------------
  !> initialize the table
  subroutine fft_all_init()
    integer :: ii

    PUSH_SUB(fft_all_init)

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
    do ii = 1, FFT_MAX
      fft_refs(ii) = FFT_NULL
    end do

    !%Variable FFTPreparePlan
    !%Type integer
    !%Default fftw_measure
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
    if(.not. varinfo_valid_option('FFTPreparePlan', fft_prepare_plan)) call input_error('FFTPreparePlan')

    POP_SUB(fft_all_init)
  end subroutine fft_all_init


  ! ---------------------------------------------------------
  !> delete all plans
  subroutine fft_all_end()
    integer :: ii

    PUSH_SUB(fft_all_end)

    do ii = 1, FFT_MAX
      if(fft_refs(ii) /= FFT_NULL) then
        call DFFTW(destroy_plan) (fft_array(ii)%planf)
        call DFFTW(destroy_plan) (fft_array(ii)%planb)
        fft_refs(ii) = FFT_NULL
      end if
    end do

    call DFFTW(cleanup)
    POP_SUB(fft_all_end)
  end subroutine fft_all_end


  ! ---------------------------------------------------------
  subroutine fft_init(nn, dim, is_real, fft, optimize)
    integer,           intent(inout) :: nn(1:3)
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: is_real
    type(fft_t),       intent(out)   :: fft
    logical, optional, intent(in)    :: optimize

    FLOAT, allocatable :: rin(:, :, :)
    CMPLX, allocatable :: cin(:, :, :), cout(:, :, :)

    integer :: ii, jj, fft_dim, idir
    logical :: optimize_
    character(len=100) :: str_tmp

    PUSH_SUB(fft_init)

    ! First, figure out the dimensionality of the FFT.
    fft_dim = 0
    do ii = 1, dim
      if(nn(ii) <= 1) exit
      fft_dim = fft_dim + 1
    end do

    if(fft_dim .eq. 0) then
      message(1) = "Internal error in fft_init: apparently, a 1x1x1 FFT is required."
      call messages_fatal(1)
    end if

    if(fft_dim > 3) call messages_not_implemented('FFT for dimension > 3')

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
    do ii = FFT_MAX, 1, -1
      if(fft_refs(ii) /= FFT_NULL) then
        if(all(nn(1:dim) == fft_array(ii)%n(1:dim)) .and. is_real == fft_array(ii)%is_real) then
          fft = fft_array(ii)              ! return a copy
          fft_refs(ii) = fft_refs(ii) + 1  ! increment the ref count
          POP_SUB(fft_init)
          return
        end if
      else
        jj = ii
      end if
    end do

    if(jj == 0) then
      message(1) = "Not enough slots for FFTs."
      message(2) = "Please increase FFT_MAX in fftw3.F90 and recompile."
      call messages_fatal(2)
    end if

    ! jj now contains an empty slot
    fft_refs(jj)          = 1
    fft_array(jj)%slot    = jj
    fft_array(jj)%n(1:3)  = nn(1:3)
    fft_array(jj)%is_real = is_real
    if(is_real == fft_real) then
      SAFE_ALLOCATE(rin(1:nn(1), 1:nn(2), 1:nn(3)))
      SAFE_ALLOCATE(cout(1:nn(1)/2+1, 1:nn(2), 1:nn(3)))

      select case(fft_dim)
      case(3)
        call DFFTW(plan_dft_r2c_3d) (fft_array(jj)%planf, nn(1), nn(2), nn(3), rin, cout, fft_prepare_plan+fftw_unaligned)
        call DFFTW(plan_dft_c2r_3d) (fft_array(jj)%planb, nn(1), nn(2), nn(3), cout, rin, fft_prepare_plan+fftw_unaligned)
      case(2)
        call DFFTW(plan_dft_r2c_2d) (fft_array(jj)%planf, nn(1), nn(2), rin, cout, fft_prepare_plan+fftw_unaligned)
        call DFFTW(plan_dft_c2r_2d) (fft_array(jj)%planb, nn(1), nn(2), cout, rin, fft_prepare_plan+fftw_unaligned)
      case(1)
        call DFFTW(plan_dft_r2c_1d) (fft_array(jj)%planf, nn(1), rin, cout, fft_prepare_plan+fftw_unaligned)
        call DFFTW(plan_dft_c2r_1d) (fft_array(jj)%planb, nn(1), cout, rin, fft_prepare_plan+fftw_unaligned)
      end select
      SAFE_DEALLOCATE_A(rin)
      SAFE_DEALLOCATE_A(cout)
    else
      SAFE_ALLOCATE( cin(1:nn(1), 1:nn(2), 1:nn(3)))
      SAFE_ALLOCATE(cout(1:nn(1), 1:nn(2), 1:nn(3)))
      select case(fft_dim)
      case(3)
        call DFFTW(plan_dft_3d) (fft_array(jj)%planf, nn(1), nn(2), nn(3), cin, cout, fftw_forward,  fft_prepare_plan)
        call DFFTW(plan_dft_3d) (fft_array(jj)%planb, nn(1), nn(2), nn(3), cin, cout, fftw_backward, fft_prepare_plan)
      case(2)
        call DFFTW(plan_dft_2d) (fft_array(jj)%planf, nn(1), nn(2), cin, cout, fftw_forward,  fft_prepare_plan)
        call DFFTW(plan_dft_2d) (fft_array(jj)%planb, nn(1), nn(2), cin, cout, fftw_backward, fft_prepare_plan)
      case(1)
        call DFFTW(plan_dft_1d) (fft_array(jj)%planf, nn(1), cin, cout, fftw_forward,  fft_prepare_plan)
        call DFFTW(plan_dft_1d) (fft_array(jj)%planb, nn(1), cin, cout, fftw_backward, fft_prepare_plan)
      end select
      SAFE_DEALLOCATE_A(cin)
      SAFE_DEALLOCATE_A(cout)
    end if
    fft = fft_array(jj)

    write(message(1), '(a)') "Info: FFT allocated with size ("
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
    call messages_info(1)

    POP_SUB(fft_init)
  end subroutine fft_init


  ! ---------------------------------------------------------
  subroutine fft_copy(fft_i, fft_o)
    type(fft_t), intent( in) :: fft_i
    type(fft_t), intent(out) :: fft_o

    PUSH_SUB(fft_copy)

    ASSERT(fft_i%slot>=1.and.fft_i%slot<=FFT_MAX)
    ASSERT(fft_refs(fft_i%slot) > 0)

    fft_o = fft_i
    fft_refs(fft_i%slot) = fft_refs(fft_i%slot) + 1

    POP_SUB(fft_copy)
  end subroutine fft_copy


  ! ---------------------------------------------------------
  subroutine fft_end(fft)
    type(fft_t), intent(inout) :: fft

    integer :: ii

    PUSH_SUB(fft_end)

    ii = fft%slot
    if(fft_refs(ii) == FFT_NULL) then
      message(1) = "Trying to deallocate FFT that has not been allocated."
      call messages_warning(1)
    else
      if(fft_refs(ii) > 1) then
        fft_refs(ii) = fft_refs(ii) - 1
      else
        fft_refs(ii) = FFT_NULL
        call DFFTW(destroy_plan) (fft_array(ii)%planf)
        call DFFTW(destroy_plan) (fft_array(ii)%planb)
        write(message(1), '(a,i4)') "Info: FFT deallocated from slot ", ii
        call messages_info(1)
      end if
    end if

    POP_SUB(fft_end)
  end subroutine fft_end


  ! ---------------------------------------------------------
  subroutine fft_getdim_real(fft, dd)
    type(fft_t), intent(in)  :: fft
    integer,     intent(out) :: dd(1:3)

    PUSH_SUB(fft_getdim_real)
    dd = fft%n

    POP_SUB(fft_getdim_real)
  end subroutine fft_getdim_real


  ! ---------------------------------------------------------
  subroutine fft_getdim_complex(fft, dd)
    type(fft_t), intent(in)  :: fft
    integer,     intent(out) :: dd(1:3)

    PUSH_SUB(fft_getdim_complex)

    dd = fft%n
    if(fft%is_real == fft_real)  dd(1) = dd(1)/2 + 1

    POP_SUB(fft_getdim_complex)
  end subroutine fft_getdim_complex


  ! ---------------------------------------------------------
  ! these routines simply call fftw
  ! first the real to complex versions
  subroutine dfft_forward(fft, rr, cc)
    type(fft_t), intent(in)  :: fft
    FLOAT,       intent(in)  :: rr(:,:,:)
    CMPLX,       intent(out) :: cc(:,:,:)

    PUSH_SUB(dfft_forward)

    call DFFTW(execute_dft_r2c) (fft%planf, rr, cc)

    POP_SUB(dfft_forward)

  end subroutine dfft_forward


  ! ---------------------------------------------------------
  ! these routines simply call fftw
  ! first the real to complex versions
  subroutine dfft_forward1(fft, rr, cc)
    type(fft_t), intent(in)  :: fft
    FLOAT,       intent(in)  :: rr(:)
    CMPLX,       intent(out) :: cc(:)

    PUSH_SUB(dfft_forward1)

    call DFFTW(execute_dft_r2c) (fft%planf, rr, cc)

    POP_SUB(dfft_forward1)

  end subroutine dfft_forward1


  ! ---------------------------------------------------------
  subroutine dfft_backward(fft, cc, rr)
    type(fft_t), intent(in) :: fft
    CMPLX, intent(in)  :: cc(fft%n(1), fft%n(2), fft%n(3))
    FLOAT, intent(out) :: rr(fft%n(1), fft%n(2), fft%n(3))

    PUSH_SUB(dfft_backward)
    
    call DFFTW(execute_dft_c2r) (fft%planb, cc, rr)

    ! multiply by 1/(N1*N2*N2)
    call lalg_scal(fft%n(1), fft%n(2), fft%n(3), &
      M_ONE / (fft%n(1)*fft%n(2)*fft%n(3)), rr)

    POP_SUB(dfft_backward)

  end subroutine dfft_backward


  ! ---------------------------------------------------------
  subroutine dfft_backward1(fft, cc, rr)
    type(fft_t), intent(in) :: fft
    CMPLX, intent(in)  :: cc(fft%n(1))
    FLOAT, intent(out) :: rr(fft%n(1))

    PUSH_SUB(dfft_backward1)
    
    call DFFTW(execute_dft_c2r) (fft%planb, cc, rr)

    ! multiply by 1/(N1*N2*N2)
    call lalg_scal(fft%n(1), &
      M_ONE / fft%n(1), rr)

    POP_SUB(dfft_backward1)

  end subroutine dfft_backward1


  ! ---------------------------------------------------------
  ! first the complex versions
  subroutine zfft_forward(fft, in, out)
    type(fft_t), intent(in)  :: fft
    CMPLX,       intent(in)  :: in(:,:,:)
    CMPLX,       intent(out) :: out(:,:,:)

    PUSH_SUB(zfft_forward)
    
    call DFFTW(execute_dft) (fft%planf, in, out)
    
    POP_SUB(zfft_forward)
    
  end subroutine zfft_forward


  ! ---------------------------------------------------------
  subroutine zfft_backward(fft, in, out)
    type(fft_t), intent(in)  :: fft
    CMPLX,       intent(in)  :: in (fft%n(1), fft%n(2), fft%n(3))
    CMPLX,       intent(out) :: out(fft%n(1), fft%n(2), fft%n(3))

    PUSH_SUB(zfft_backward)

    call DFFTW(execute_dft) (fft%planb, in, out)

    ! multiply by 1/(N1*N2*N2)
    call lalg_scal(fft%n(1), fft%n(2), fft%n(3), &
      M_z1 / (fft%n(1)*fft%n(2)*fft%n(3)), out)

    POP_SUB(zfft_backward)

  end subroutine zfft_backward


  ! ---------------------------------------------------------
  ! convert between array index and G vector
  function pad_feq(ii, nn, mode)
    integer, intent(in) :: ii,nn
    logical, intent(in) :: mode
    integer :: pad_feq

    ! no push_sub: called too frequently

    if(mode) then      ! index to frequency number
      if( ii <= nn/2 + 1 ) then
        pad_feq = ii - 1
      else
        pad_feq = ii - nn - 1
      end if
    else
      if( ii >= 0 ) then
        pad_feq = ii + 1
      else
        pad_feq = ii + nn + 1
      end if
    end if

    return
  end function pad_feq

end module fft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
