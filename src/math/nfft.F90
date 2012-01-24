!! Copyright (C) 2011 U. De Giovannini
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
!! $Id: fftw3.F90 7027 2010-10-11 22:47:20Z dstrubbe $


#include "global.h"

#define NFFT_MAX 10
#define NFFT_NULL -1


module nfft_m

#if !defined(HAVE_NFFT) 
  integer, public :: nfft_dummy ! this avoids compilers complaining about empty module 
#else

  use c_pointer_m
  use datasets_m
  use fft_m   
  use fftw_m
  use global_m
  use iso_c_binding
  use loct_math_m
  use messages_m
  use parser_m
  use varinfo_m

  implicit none

  private
  public ::        &
    nfft_t,         &
    nfft_init,      &
    nfft_end,       &
    nfft_precompute,&
    znfft_forward,  &
    znfft_backward

  ! global constants
  integer, public, parameter ::         &
    nfft_real    = 0,                   &
    nfft_complex = 1

  !NFFT flags
  integer, public, parameter ::        &
    nfft_PRE_PHI_HUT       =        0, &
    nfft_FG_PSI            =        2, &
    nfft_PRE_LIN_PSI       =        4, &
    nfft_PRE_FG_PSI        =        8, &
    nfft_PRE_PSI           =       16, &
    nfft_PRE_FULL_PSI      =       32, &
    nfft_MALLOC_X          =       64, &
    nfft_MALLOC_F_HAT      =      128, &
    nfft_MALLOC_F          =      256, &
    nfft_FFT_OUT_OF_PLACE  =      512, &
    nfft_FFTW_INIT         =     1024


  type nfft_t
    integer           :: slot       ! in which slot do we have this fft

    integer           :: N(MAX_DIM) ! size of the nfft bandwidths 
    integer           :: M          ! Number of the nfft nodes 
    integer           :: is_real    ! is the fft real or complex
    integer           :: dim        ! the dimension 
    ! Guru options
    logical           :: guru       ! use guru options
    integer           :: precompute ! precompute strategy
    integer           :: mm         ! Window function cut-off parameter 
    FLOAT             :: sigma      ! Oversampling factor 

    type(c_ptr)       :: plan       ! the plan    

  end type nfft_t



  integer      :: nfft_refs(NFFT_MAX)
  type(nfft_t) :: nfft_array(NFFT_MAX)
  logical      :: nfft_optimize
  integer      :: nfft_prepare_plan


interface
! Direct access to NFFT C-functions
  subroutine nfft_init_1d(plan, N, M) &
     bind(C,name='nfft_init_1d')
      import 
      type(c_ptr),          intent(inout)   :: plan
      integer(c_int), value, intent (in)    :: N
      integer(c_int), value, intent (in)    :: M
!      integer(c_int),  intent (in)    :: N
!      integer(c_int),  intent (in)    :: M
  end subroutine nfft_init_1d

  subroutine nfft_init_2d(plan, N1, N2, M) &
     bind(C,name='nfft_init_2d')
      import 
      type(c_ptr),          intent(inout)   :: plan
      integer(c_int), value, intent (in)    :: N1
      integer(c_int), value, intent (in)    :: N2
      integer(c_int), value, intent (in)    :: M
  end subroutine nfft_init_2d

  subroutine nfft_init_3d(plan, N1, N2,N3, M) &
     bind(C,name='nfft_init_3d')
      import 
      type(c_ptr),          intent(inout)   :: plan
      integer(c_int), value, intent (in)    :: N1
      integer(c_int), value, intent (in)    :: N2
      integer(c_int), value, intent (in)    :: N3
      integer(c_int), value, intent (in)    :: M
  end subroutine nfft_init_3d

!void nfft_init_guru(nfft_plan *ths, int d, int *N, int M_total, int *n,
!                        int m, unsigned nfft_flags, unsigned fftw_flags)

  subroutine nfft_init_guru(plan, d, N, M_total, nn, mm, nfft_flags, fftw_flags) &
     bind(C,name='nfft_init_guru')
      import 
      type(c_ptr),          intent(inout)   :: plan
      integer(c_int), value, intent (in)    :: d
      integer(c_int),        intent (in)    :: N(*)
      integer(c_int), value, intent (in)    :: M_total
      integer(c_int),        intent (in)    :: nn(*)
      integer(c_int), value, intent (in)    :: mm
      integer(c_int), value, intent (in)    :: nfft_flags
      integer(c_int), value, intent (in)    :: fftw_flags
  end subroutine nfft_init_guru

  subroutine nfft_check(plan) &
     bind(C,name='nfft_check')
      import 
      type(c_ptr),          intent(inout)   :: plan
  end subroutine nfft_check

  subroutine nfft_finalize(plan) &
    bind(c,name='nfft_finalize')
      import
      type(c_ptr),   intent(inout) :: plan
  end subroutine nfft_finalize

  subroutine nfft_adjoint(plan) &
    bind(c,name='nfft_adjoint')
      import
      type(c_ptr),   intent(in) :: plan
  end subroutine nfft_adjoint

  subroutine nfft_trafo(plan) &
    bind(c,name='nfft_trafo')
      import
      type(c_ptr),   intent(in) :: plan
  end subroutine nfft_trafo

! OCTOPUS NFFT C-wrappers and helper functions

  subroutine oct_nfft_precompute_one_psi_1d(plan, M, X1) &
    bind(c,name='oct_nfft_precompute_one_psi_1d')
      import
      type(c_ptr),            intent(inout) :: plan
      integer(c_int), value,  intent(in)    :: M
      real(c_double),         intent(in)    :: X1(*)
  end subroutine oct_nfft_precompute_one_psi_1d

  subroutine oct_nfft_precompute_one_psi_2d(plan, M, X1, X2) &
    bind(c,name='oct_nfft_precompute_one_psi_2d')
      import
      type(c_ptr),            intent(inout) :: plan
      integer(c_int), value,  intent(in)    :: M
      real(c_double),         intent(in)    :: X1(*)
      real(c_double),         intent(in)    :: X2(*)
  end subroutine oct_nfft_precompute_one_psi_2d

  subroutine oct_nfft_precompute_one_psi_3d(plan, M, X1, X2, X3) &
    bind(c,name='oct_nfft_precompute_one_psi_3d')
      import
      type(c_ptr),            intent(inout) :: plan
      integer(c_int), value,  intent(in)    :: M
      real(c_double),         intent(in)    :: X1(*)
      real(c_double),         intent(in)    :: X2(*)
      real(c_double),         intent(in)    :: X3(*)
  end subroutine oct_nfft_precompute_one_psi_3d

  subroutine oct_znfft_forward(plan, M, in, out) &
    bind(c,name='oct_znfft_forward')
      import
      type(c_ptr),                intent(in)  :: plan
      integer(c_int), value,      intent(in)  :: M
      complex(c_double_complex),  intent(in)  :: in(*)
      complex(c_double_complex),  intent(out) :: out(*)
  end subroutine oct_znfft_forward

  subroutine oct_set_f(plan, M, dim, val, ix, iy, iz) &
    bind(c,name='oct_set_f')
      import
      type(c_ptr),                      intent(in) :: plan
      integer(c_int), value,            intent(in) :: M
      integer(c_int), value,            intent(in) :: dim
      complex(c_double_complex), value, intent(in) :: val
      integer(c_int), value,            intent(in) :: ix
      integer(c_int), value,            intent(in) :: iy
      integer(c_int), value,            intent(in) :: iz
  end subroutine oct_set_f

  subroutine oct_get_f(plan, M, dim, val, ix, iy, iz) &
    bind(c,name='oct_get_f')
      import
      type(c_ptr),               intent(in) :: plan
      integer(c_int), value,     intent(in) :: M
      integer(c_int), value,     intent(in) :: dim
      complex(c_double_complex), intent(in) :: val
      integer(c_int), value,     intent(in) :: ix
      integer(c_int), value,     intent(in) :: iy
      integer(c_int), value,     intent(in) :: iz
  end subroutine oct_get_f

  subroutine oct_set_f_hat(plan, M, dim, val, ix, iy, iz) &
    bind(c,name='oct_set_f_hat')
      import
      type(c_ptr),                      intent(in) :: plan
      integer(c_int), value,            intent(in) :: M
      integer(c_int), value,            intent(in) :: dim
      complex(c_double_complex), value, intent(in) :: val
      integer(c_int), value,            intent(in) :: ix
      integer(c_int), value,            intent(in) :: iy
      integer(c_int), value,            intent(in) :: iz
  end subroutine oct_set_f_hat

  subroutine oct_get_f_hat(plan, M, dim, val, ix, iy, iz) &
    bind(c,name='oct_get_f_hat')
      import
      type(c_ptr),               intent(in) :: plan
      integer(c_int), value,     intent(in) :: M
      integer(c_int), value,     intent(in) :: dim
      complex(c_double_complex), intent(in) :: val
      integer(c_int), value,     intent(in) :: ix
      integer(c_int), value,     intent(in) :: iy
      integer(c_int), value,     intent(in) :: iz
  end subroutine oct_get_f_hat

end interface

contains

  ! ---------------------------------------------------------
  ! GURU options
  subroutine nfft_guru_options(nfft)
    type(nfft_t), intent(inout) :: nfft

    PUSH_SUB(nfft_guru_options)

    !%Variable NFFTGuruInterface
    !%Type logical
    !%Default false
    !%Section Mesh::FFTs
    !%Description
    !% Perform NFFT with guru interface. This permits the fine tuning of several critical parameters.
    !%End
    call parse_logical(datasets_check('NFFTGuruInterface'), .false., nfft%guru)
 

    !%Variable NFFTCutoff
    !%Type integer
    !%Default 6
    !%Section Mesh::FFTs
    !%Description
    !% Cut-off parameter of the window function. 
    !% See NFFT manual for details.
    !%End
    call parse_integer(datasets_check('NFFTCutoff'), 6, nfft%mm)


    !%Variable NFFTOversampling
    !%Type float
    !%Default 2
    !%Section Mesh::FFTs
    !%Description
    !% NFFT oversampling factor (sigma). This will rule the size of the FFT under the hood.
    !%End
    call parse_float(datasets_check('NFFTOversampling'), M_TWO, nfft%sigma)

    !%Variable NFFTPrecompute
    !%Type integer
    !%Default NFFT_PRE_PSI
    !%Section Mesh::FFTs
    !%Description
    !% NFFT precomputation strategy.
    !%Option NFFT_PRE_LIN_PSI 4
    !% This method implements a linear interpolation from a lookup table.
    !%Option NFFT_PRE_PSI 16
    !% This method uses a medium amount of memory to store d*(2*m+1)*M real numbers and requires at most 
    !% 2(2m + 1)d extra multiplications for each node.
    !% This is the default option.
    !%Option NFFT_PRE_FULL_PSI 32
    !% Is the fastest method but requires a large amount of memory as it requires to store (2*m+1)^d*M  
    !% real numbers. No extra operations are needed during matrix vector multiplication.
    !%End
    call parse_integer(datasets_check('NFFTPrecompute'), 16, nfft%precompute)
     if(.not.varinfo_valid_option('NFFTPrecompute', nfft%precompute)) call input_error('NFFTPrecompute')
!    call messages_print_var_option(stdout, "NFFTPrecompute", nfft%precompute)

!    if(.not.varinfo_valid_option('NFFTPrecompute', nfft%precompute, is_flag=.true.)) then
!      call input_error('NFFTPrecompute')
 !   end if


    POP_SUB(nfft_guru_options)
  end subroutine nfft_guru_options



  ! ---------------------------------------------------------
  subroutine nfft_init(nn, dim, mm, is_real, nfft, optimize)
    integer,           intent(inout) :: nn(1:MAX_DIM)
    integer,           intent(inout) :: mm
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: is_real
    type(nfft_t),      intent(inout) :: nfft
    logical, optional, intent(in)    :: optimize

    integer :: ii, jj, nfft_dim, idir,my_nn(1:MAX_DIM)
    logical :: optimize_
    character(len=100) :: str_tmp
    integer :: nfft_flags


    PUSH_SUB(nfft_init)

    nfft_dim = 0
    do ii = 1, dim
      if(nn(ii) <= 1) exit
      nfft_dim = nfft_dim + 1
!      nn(ii)=2**ceiling(log(real(nn(ii)))/log(M_TWO))
      if(int(nn(ii)/2)*2 .ne. nn(ii)) nn(ii)=nn(ii)+1
    end do

    if(nfft_dim .eq. 0) then
      message(1) = "Internal error in fft_init: apparently, a 1x1x1 NFFT is required."
      call messages_fatal(1)
    end if

    if(nfft_dim .gt. 3) then
      message(1) = "NFFT for dimension greater than 3 not implemented."
      call messages_fatal(1)
    end if


    nfft%dim=nfft_dim
    nfft%M = mm
    nfft%N = nn   


    call nfft_guru_options(nfft)


    ! Why do we not use the value of the flag optimize here?? -DAS
    if(nfft%guru) then 

      do ii = 1, nfft_dim
        my_nn(ii) = nn(ii)*nfft%sigma
        call loct_fft_optimize(my_nn(ii), 1) ! ask for an odd number
      end do
      nfft_flags =  nfft_PRE_PHI_HUT  + nfft_MALLOC_X +nfft_MALLOC_F_HAT +&
                    nfft_MALLOC_F + nfft_FFTW_INIT + nfft_FFT_OUT_OF_PLACE

      nfft_flags = nfft_flags + nfft%precompute

      call  nfft_init_guru(nfft%plan, nfft_dim, nn, mm**nfft_dim, my_nn, nfft%mm, &
                    nfft_flags, FFTW_MEASURE + FFTW_DESTROY_INPUT)

    else

      select case(nfft_dim)
      case(3)
        call nfft_init_3d(nfft%plan, nn(1), nn(2),nn(3), mm*mm*mm)
      case(2)
        call nfft_init_2d(nfft%plan, nn(1), nn(2), mm*mm)
      case(1)
        call  nfft_init_1d(nfft%plan,nn(1),mm)
      end select

    end if


    write(message(1), '(a)') "Info: NFFT allocated with size ("
    do idir = 1, nfft_dim
      write(str_tmp, '(i7,a)') nn(idir)
      if(idir == nfft_dim) then
        message(1) = trim(message(1)) // trim(str_tmp) // ";; "
      else
        message(1) = trim(message(1)) // trim(str_tmp) // ","
      end if
    end do
    write(str_tmp, '(i7,a)') mm
    message(1) = trim(message(1)) // trim(str_tmp) // ") "
    call messages_info(1)


    if(nfft%guru) then

      write(message(1), '(a)') "Info: NFFT use FFT with size ("
      do idir = 1, nfft_dim
        write(str_tmp, '(i7,a)') my_nn(idir)
        if(idir == nfft_dim) then
          message(1) = trim(message(1)) // trim(str_tmp) // ") "
        else
          message(1) = trim(message(1)) // trim(str_tmp) // ","
        end if
      end do
      call messages_info(1)

      write(message(1), '(a,i3)') "Info: NFFT window function cut-off parameter m = ",nfft%mm
      call messages_info(1)

      write(message(1), '(a)') "Info: NFFT precomputation strategy "
      select case(nfft%precompute)
      case(nfft_PRE_LIN_PSI)
        write(str_tmp, '(a)') " NFFT_PRE_LIN_PSI" 
      case(nfft_PRE_PSI)
        write(str_tmp, '(a)') " NFFT_PRE_PSI" 
      case(nfft_PRE_FULL_PSI)
        write(str_tmp, '(a)') " NFFT_PRE_FULL_PSI" 
      end select
      message(1) = trim(message(1)) // trim(str_tmp)     
      call messages_info(1)


    end if


    POP_SUB(nfft_init)
  end subroutine nfft_init

  ! ---------------------------------------------------------
  subroutine nfft_end(nfft)
    type(nfft_t), intent(inout) :: nfft

    integer :: ii

    PUSH_SUB(nfft_end)

    call nfft_finalize(nfft%plan);

    message(1) = "Info: NFFT deallocated."
    call messages_info(1)

   
    POP_SUB(nfft_end)
  end subroutine nfft_end


  !----------------------------------------------------------
  ! We allow different spacing for each direction x,y,z
  ! ---------------------------------------------------------
  subroutine nfft_precompute(nfft, X1,X2,X3)
    FLOAT,           intent(in)    :: X1(:)
    FLOAT, optional, intent(in)    :: X2(:)
    FLOAT, optional, intent(in)    :: X3(:)
    type(nfft_t),    intent(inout) :: nfft


    real(c_double) :: x1_(1:nfft%M), x2_(1:nfft%M), x3_(1:nfft%M)
    FLOAT :: length, cc, eps

    PUSH_SUB(nfft_precompute)
 
    eps = 1.000001 ! the sample nodes must be in [0.5,0.5)
 
    select case(nfft%dim)
      case(3)
        length = (maxval(X1)-minval(X1))*eps
        cc = (minval(X1)+maxval(X1))/M_TWO
        x1_ =(X1-cc)/length
        length = (maxval(X2)-minval(X2))*eps
        cc = (minval(X2)+maxval(X2))/M_TWO
        x2_ =(X2-cc)/length
        length = (maxval(X3)-minval(X3))*eps
        cc = (minval(X3)+maxval(X3))/M_TWO
        x3_ =(X3-cc)/length
        call oct_nfft_precompute_one_psi_3d(nfft%plan, nfft%M, x1_, x2_, x3_)

      case(2)
        length = (maxval(X1)-minval(X1))*eps
        cc = (minval(X1)+maxval(X1))/M_TWO
        x1_ =(X1-cc)/length
        length = (maxval(X2)-minval(X2))*eps
        cc = (minval(X2)+maxval(X2))/M_TWO
        x2_ =(X2-cc)/length
        call oct_nfft_precompute_one_psi_2d(nfft%plan, nfft%M, x1_, x2_)

      case(1)
        length = (maxval(X1)-minval(X1))*eps
        cc = (minval(X1)+maxval(X1))/M_TWO
        x1_ =(X1-cc)/length
        call oct_nfft_precompute_one_psi_1d(nfft%plan,nfft%M,x1_)
 
    end select

    ! check the plan
    call nfft_check(nfft%plan)


    write(message(1), '(a)') "Info: NFFT plan precomputed."
    call messages_info(1)


    PUSH_SUB(nfft_precompute)
  end subroutine nfft_precompute

  subroutine znfft_forward(nfft, in, out)
    type(nfft_t), intent(in)  :: nfft
    CMPLX,        intent(in)  :: in(:,:,:)
    CMPLX,        intent(out) :: out(:,:,:)

    integer :: ix, iy, iz, MM(MAX_DIM)

    PUSH_SUB(znfft_forward)
    
    select case(nfft%dim)
      case (1)
        MM(1) = nfft%M
        MM(2) = 1
        MM(3) = 1

      case (2)
        MM(1) = nfft%M
        MM(2) = nfft%M
        MM(3) = 1

      case (3)
        MM(1) = nfft%M
        MM(2) = nfft%M
        MM(3) = nfft%M

    end select 

    do ix = 1, nfft%N(1)
      do iy = 1, nfft%N(2)
        do iz = 1, nfft%N(3)
          call oct_set_f_hat(nfft%plan, nfft%N(1), nfft%dim, in(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do

    call nfft_trafo(nfft%plan)

    do ix = 1, MM(1)
      do iy = 1, MM(2)
        do iz = 1, MM(3)
          call oct_get_f(nfft%plan, nfft%M, nfft%dim, out(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do

    POP_SUB(znfft_forward)
    
  end subroutine znfft_forward


  ! ---------------------------------------------------------
  subroutine znfft_backward(nfft, in, out)
    type(nfft_t), intent(in)  :: nfft
    CMPLX,        intent(in)  :: in (:,:,:)
    CMPLX,        intent(out) :: out(:,:,:)

    integer :: ix, iy, iz, MM(MAX_DIM)

    PUSH_SUB(znfft_backward)

    select case(nfft%dim)
      case (1)
        MM(1) = nfft%M
        MM(2) = 1
        MM(3) = 1

      case (2)
        MM(1) = nfft%M
        MM(2) = nfft%M
        MM(3) = 1

      case (3)
        MM(1) = nfft%M
        MM(2) = nfft%M
        MM(3) = nfft%M

    end select

    do ix = 1, MM(1)
      do iy = 1, MM(2)
        do iz = 1, MM(3)
          call oct_set_f(nfft%plan, nfft%M, nfft%dim, in(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do

    call nfft_adjoint(nfft%plan)

    do ix = 1,nfft%N(1)
      do iy = 1, nfft%N(2)
        do iz = 1, nfft%N(3)
          call oct_get_f_hat(nfft%plan, nfft%M, nfft%dim, out(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do


!    call lalg_scal(nfft%N(1), nfft%N(2), nfft%N(3), &
!      M_z1 / (nfft%N(1)*nfft%N(2)*nfft%N(3)), out)

!   out =out /nfft%N(1)

    POP_SUB(znfft_backward)

  end subroutine znfft_backward

#endif

end module nfft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
