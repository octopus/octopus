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

module nfft_m

#if !defined(HAVE_NFFT) 
  integer, public :: nfft_dummy ! this avoids compilers complaining about empty module 
#else


  use c_pointer_m
  use datasets_m
  use fftw_m
  use global_m
  use loct_math_m
  use messages_m
  use parser_m
  use varinfo_m

  implicit none

  private
  public ::          &
    nfft_t,          &
    nfft_copy_info,  &
    nfft_init,       &
    nfft_end,        &
    nfft_precompute, &
    nfft_write_info, &
    znfft_forward,   &
    znfft_backward,  &
    dnfft_forward,   &
    dnfft_backward

  ! global constants
  integer, public, parameter ::         &
    nfft_real    = 0,                   &
    nfft_complex = 1

  !NFFT flags
  integer, public, parameter ::        &
    NFFT_PRE_PHI_HUT       =        0, &
    NFFT_FG_PSI            =        2, &
    NFFT_PRE_LIN_PSI       =        4, &
    NFFT_PRE_FG_PSI        =        8, &
    NFFT_PRE_PSI           =       16, &
    NFFT_PRE_FULL_PSI      =       32, &
    NFFT_MALLOC_X          =       64, &
    NFFT_MALLOC_F_HAT      =      128, &
    NFFT_MALLOC_F          =      256, &
    NFFT_FFT_OUT_OF_PLACE  =      512, &
    NFFT_FFTW_INIT         =     1024


  type nfft_t

    integer           :: N(3)       !> size of the nfft bandwidths 
    integer           :: M          !> Number of the nfft nodes 
    integer           :: is_real    !> is the fft real or complex
    integer           :: dim        !> the dimension 
    integer           :: fftN(3)    !> size of the fft used   
    FLOAT             :: norm       !> Normalization  

    ! Guru options
    logical           :: set_defaults = .false. !> set default values from the code
    logical           :: guru                   !> use guru options?
    integer           :: precompute             !> precompute strategy
    integer           :: mm                     !> Window function cut-off parameter 
    FLOAT             :: sigma                  !> Oversampling factor 
         
    type(c_ptr)       :: plan                   !> the NFFT plan    

  end type nfft_t


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
    call parse_logical(datasets_check('NFFTGuruInterface'),  nfft%guru, nfft%guru)
 

    !%Variable NFFTCutoff
    !%Type integer
    !%Default 6
    !%Section Mesh::FFTs
    !%Description
    !% Cut-off parameter of the window function. 
    !% See NFFT manual for details.
    !%End
    call parse_integer(datasets_check('NFFTCutoff'), nfft%mm, nfft%mm)


    !%Variable NFFTOversampling
    !%Type float
    !%Default 2
    !%Section Mesh::FFTs
    !%Description
    !% NFFT oversampling factor (sigma). This will rule the size of the FFT under the hood.
    !%End
    call parse_float(datasets_check('NFFTOversampling'), nfft%sigma, nfft%sigma)

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
    call parse_integer(datasets_check('NFFTPrecompute'), nfft%precompute, nfft%precompute)
     if(.not.varinfo_valid_option('NFFTPrecompute', nfft%precompute)) call input_error('NFFTPrecompute')
!    call messages_print_var_option(stdout, "NFFTPrecompute", nfft%precompute)

!     if(.not.varinfo_valid_option('NFFTPrecompute', nfft%precompute, is_flag=.true.)) then
!       call input_error('NFFTPrecompute')
!     end if


    POP_SUB(nfft_guru_options)
  end subroutine nfft_guru_options



  ! ---------------------------------------------------------
  subroutine nfft_init(nfft, nn, dim, mm, is_real, optimize)
    type(nfft_t),      intent(inout) :: nfft
    integer,           intent(inout) :: nn(3) !> nfft bandwidths 
    integer,           intent(inout) :: mm    !> nfft nodes 
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: is_real
    logical, optional, intent(in)    :: optimize

    integer :: ii, jj, idir, my_nn(3)
    logical :: optimize_
    integer :: nfft_flags


    PUSH_SUB(nfft_init)

    optimize_ = optional_default(optimize, .true.)

    nfft%dim = dim
    nfft%M = mm
    nfft%N = nn   

    if(.not. nfft%set_defaults) then
      !Set defaults
      nfft%guru = .false.
      nfft%mm = 6 
      nfft%sigma = M_TWO
      nfft%precompute = NFFT_PRE_PSI
    end if
    
    call nfft_guru_options(nfft)

    my_nn = 0
    do ii = 1, dim
      my_nn(ii) = nn(ii)*nfft%sigma
      if(optimize_ .or. (.not. nfft%guru)) call loct_fft_optimize(my_nn(ii), 1) ! ask for an odd number
    end do
    
    nfft%fftN(1:dim) = my_nn(1:dim)

    if(nfft%guru) then ! Guru interface 
      nfft_flags =  nfft_PRE_PHI_HUT  + nfft_MALLOC_X +nfft_MALLOC_F_HAT +&
                    nfft_MALLOC_F + nfft_FFTW_INIT + nfft_FFT_OUT_OF_PLACE

      nfft_flags = nfft_flags + nfft%precompute
      call  oct_nfft_init_guru(nfft%plan, dim, nn, mm**dim, my_nn, nfft%mm, &
                    nfft_flags, FFTW_MEASURE + FFTW_DESTROY_INPUT)

    else ! Default interfaces

      select case(dim) 
      case(3)
        call oct_nfft_init_3d(nfft%plan, nn(1), nn(2),nn(3), mm*mm*mm)
      case(2)
        call oct_nfft_init_2d(nfft%plan, nn(1), nn(2), mm*mm)
      case(1)
        call oct_nfft_init_1d(nfft%plan,nn(1),mm)
      end select

    end if


    POP_SUB(nfft_init)
  end subroutine nfft_init

  ! ---------------------------------------------------------
  subroutine nfft_write_info(nfft)
    type(nfft_t), intent(inout) :: nfft

    integer :: idir
    character(len=100) :: str_tmp

    PUSH_SUB(nfft_write_info)

    write(message(1), '(a)') "Info: NFFT dimensions ("
    do idir = 1, nfft%dim
      write(str_tmp, '(i7,a)') nfft%N(idir)
      if(idir == nfft%dim) then
        message(1) = trim(message(1)) // trim(str_tmp) // ";; "
      else
        message(1) = trim(message(1)) // trim(str_tmp) // ","
      end if
    end do
    write(str_tmp, '(i7,a)') nfft%M 
    message(1) = trim(message(1)) // trim(str_tmp) // ") "
    call messages_info(1)

    write(message(1), '(a)') "Info: NFFT uses FFTW3 with size ("
    do idir = 1, nfft%dim
      write(str_tmp, '(i7,a)') nfft%fftN(idir)
      if(idir == nfft%dim) then
        message(1) = trim(message(1)) // trim(str_tmp) // ") "
      else
        message(1) = trim(message(1)) // trim(str_tmp) // ","
      end if
    end do
    call messages_info(1)

    write(message(1), '(a,f8.4)') "Info: NFFT oversampling factor sigma = ",nfft%sigma
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

   
    POP_SUB(nfft_write_info)
  end subroutine nfft_write_info


  ! ---------------------------------------------------------
  subroutine nfft_end(nfft)
    type(nfft_t), intent(inout) :: nfft

    integer :: ii

    PUSH_SUB(nfft_end)

    call oct_nfft_finalize(nfft%plan);

   
    POP_SUB(nfft_end)
  end subroutine nfft_end

  ! ---------------------------------------------------------
  ! This routine is intend to copy the configuration parameters 
  ! rather the whole structure.
  ! ---------------------------------------------------------  
  subroutine nfft_copy_info(in, out)
    type(nfft_t), intent(in)  :: in
    type(nfft_t), intent(out) :: out


    PUSH_SUB(nfft_copy_info)

    out%N = in%N       
    out%M = in%M       
    out%is_real= in%is_real
    out%dim = in%dim
    out%fftN = in%fftN
    out%norm = in%norm

    out%set_defaults = in%set_defaults
    out%guru = in%guru            
    out%precompute = in%precompute
    out%mm = in%mm               
    out%sigma = in%sigma         

   
    POP_SUB(nfft_copy_info)
  end subroutine nfft_copy_info


  !----------------------------------------------------------
  ! Precompute the plan according to the position the grid nodes in real space 
  ! x axis is X1, y axis is X2, z axis is X3
  ! NOTE: We only allow different spacing for each direction x,y,z
  ! the NFFT interface however is more general 
  ! ---------------------------------------------------------
  subroutine nfft_precompute(nfft, X1, X2, X3)
    type(nfft_t),    intent(inout) :: nfft
    FLOAT,           intent(in)    :: X1(:)
    FLOAT, optional, intent(in)    :: X2(:)
    FLOAT, optional, intent(in)    :: X3(:)
    


    FLOAT   :: x1_(1:nfft%M), x2_(1:nfft%M), x3_(1:nfft%M)
    FLOAT   :: length, cc, eps, dX(1:nfft%M-1,1:3)
    integer :: ii

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
         
        ! Set the normalization factor  
        do ii = 1, nfft%M-1
          dX(ii,1)= abs(x1_(ii+1)-x1_(ii))
          dX(ii,2)= abs(x2_(ii+1)-x2_(ii))
          dX(ii,3)= abs(x3_(ii+1)-x3_(ii))
        end do
        nfft%norm = M_ONE/(minval(dX(:,1)) * minval(dX(:,2)) * minval(dX(:,3)))

      case(2)
        length = (maxval(X1)-minval(X1))*eps
        cc = (minval(X1)+maxval(X1))/M_TWO
        x1_ =(X1-cc)/length
        length = (maxval(X2)-minval(X2))*eps
        cc = (minval(X2)+maxval(X2))/M_TWO
        x2_ =(X2-cc)/length
        call oct_nfft_precompute_one_psi_2d(nfft%plan, nfft%M, x1_, x2_)

        ! Set the normalization factor  
        do ii = 1, nfft%M-1
          dX(ii,1)= abs(x1_(ii+1)-x1_(ii))
          dX(ii,2)= abs(x2_(ii+1)-x2_(ii))
        end do
        nfft%norm = M_ONE/(minval(dX(:,1)) * minval(dX(:,2)))


      case(1)
        length = (maxval(X1)-minval(X1))*eps
        cc = (minval(X1)+maxval(X1))/M_TWO
        x1_ =(X1-cc)/length
        call oct_nfft_precompute_one_psi_1d(nfft%plan,nfft%M,x1_)

        ! Set the normalization factor  
        do ii = 1, nfft%M-1
          dX(ii,1)= abs(x1_(ii+1)-x1_(ii))
        end do
        nfft%norm = M_ONE/(minval(dX(:,1)))
 
    end select

    ! check the plan
    call oct_nfft_check(nfft%plan)


    write(message(1), '(a)') "Info: NFFT plan precomputed."
    call messages_info(1)


    PUSH_SUB(nfft_precompute)
  end subroutine nfft_precompute


#include "undef.F90"
#include "real.F90"
#include "nfft_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "nfft_inc.F90"



#endif




end module nfft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
