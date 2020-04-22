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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!


#include "global.h"

module nfft_oct_m
  use fftw_params_oct_m
  use global_oct_m
  use, intrinsic :: iso_c_binding
  use loct_math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use varinfo_oct_m
  implicit none
  
  private
  
  public ::          &
    nfft_t,          &
    nfft_copy_info,  &
    nfft_init,       &
    nfft_end,        &
    nfft_precompute, &
    nfft_write_info, &
    nfft_guru_options, &
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
    private

    integer           :: N(3)       !> size of the nfft bandwidths
    integer           :: M(3)          !> Number of the nfft nodes
    integer           :: is_real    !> is the fft real or complex
    integer           :: dim        !> the dimension
    integer           :: fftN(3)    !> size of the fft used
    FLOAT, public     :: norm       !> Normalization

    ! Guru options
    logical, public   :: set_defaults = .false. !> the defaults can be overriden
    logical, public   :: guru                   !> use guru options?
    integer, public   :: precompute             !> precompute strategy
    integer, public   :: mm                     !> Window function cut-off parameter
    FLOAT,   public   :: sigma                  !> Oversampling factor

    type(c_ptr)       :: plan                   !> the NFFT plan

  end type nfft_t


contains



  ! ---------------------------------------------------------
  ! GURU options
  subroutine nfft_guru_options(nfft, namespace)
    type(nfft_t),      intent(inout) :: nfft
    type(namespace_t), intent(in)    :: namespace

    PUSH_SUB(nfft_guru_options)

    !%Variable NFFTGuruInterface
    !%Type logical
    !%Default false
    !%Section Mesh::FFTs
    !%Description
    !% Perform NFFT with guru interface. This permits the fine tuning of several critical parameters.
    !%End
    call parse_variable(namespace, 'NFFTGuruInterface',  .false., nfft%guru)


    !%Variable NFFTCutoff
    !%Type integer
    !%Default 6
    !%Section Mesh::FFTs
    !%Description
    !% Cut-off parameter of the window function.
    !% See NFFT manual for details.
    !%End
    call parse_variable(namespace, 'NFFTCutoff', 6, nfft%mm)


    !%Variable NFFTOversampling
    !%Type float
    !%Default 2
    !%Section Mesh::FFTs
    !%Description
    !% NFFT oversampling factor (sigma). This will rule the size of the FFT under the hood.
    !%End
    call parse_variable(namespace, 'NFFTOversampling', M_TWO, nfft%sigma)

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
    call parse_variable(namespace, 'NFFTPrecompute', NFFT_PRE_PSI, nfft%precompute)
     if(.not.varinfo_valid_option('NFFTPrecompute', nfft%precompute)) call messages_input_error(namespace, 'NFFTPrecompute')
!    call messages_print_var_option(stdout, "NFFTPrecompute", nfft%precompute)

!     if(.not.varinfo_valid_option('NFFTPrecompute', nfft%precompute, is_flag=.true.)) then
!       call messages_input_error('NFFTPrecompute')
!     end if


    POP_SUB(nfft_guru_options)
  end subroutine nfft_guru_options



  ! ---------------------------------------------------------
  subroutine nfft_init(nfft, nfft_options, N, dim, M, is_real, optimize)
    type(nfft_t),      intent(inout) :: nfft
    type(nfft_t),      intent(in)    :: nfft_options
    integer,           intent(inout) :: N(3) !> nfft bandwidths
    integer,           intent(inout) :: M(3) !> nfft nodes
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: is_real
    logical, optional, intent(in)    :: optimize

    integer :: ii, my_N(3)
    logical :: optimize_
    integer :: nfft_flags


    PUSH_SUB(nfft_init)

    optimize_ = optional_default(optimize, .true.)

    nfft%dim = dim
    nfft%M(:) = M(:)
    nfft%N(:) = N(:)

    if(.not. nfft%set_defaults) then
      nfft%guru = nfft_options%guru
      nfft%mm = nfft_options%mm
      nfft%sigma = nfft_options%sigma
      nfft%precompute = nfft_options%precompute
    end if
    
    ! set unused dimensions to 1
    nfft%M(dim+1:3) = 1

    my_N = 0
    do ii = 1, dim
      my_N(ii) = N(ii)*nfft%sigma
      if(optimize_ .or. (.not. nfft%guru)) call loct_fft_optimize(my_N(ii), 1) ! ask for an odd number
    end do

    nfft%fftN(1:dim) = my_N(1:dim)

    if(nfft%guru) then ! Guru interface
      nfft_flags =  nfft_PRE_PHI_HUT  + nfft_MALLOC_X +nfft_MALLOC_F_HAT +&
                    nfft_MALLOC_F + nfft_FFTW_INIT + nfft_FFT_OUT_OF_PLACE

      nfft_flags = nfft_flags + nfft%precompute

      call  oct_nfft_init_guru(nfft%plan, dim, N, M(1)*M(2)*M(3), my_N, nfft%mm, &
                    nfft_flags, FFTW_MEASURE + FFTW_DESTROY_INPUT)

    else ! Default interfaces

      select case(dim)
      case(3)
        call oct_nfft_init_3d(nfft%plan, N(1), N(2),N(3), M(1)*M(2)*M(3))
      case(2)
        call oct_nfft_init_2d(nfft%plan, N(1), N(2), M(1)*M(2))
      case(1)
        call oct_nfft_init_1d(nfft%plan,N(1),M(1))
      end select

    end if


    POP_SUB(nfft_init)
  end subroutine nfft_init

  ! ---------------------------------------------------------
  subroutine nfft_write_info(nfft)
    type(nfft_t), intent(inout) :: nfft

    integer :: idir
!    integer :: mm

    PUSH_SUB(nfft_write_info)

    call messages_write("Info: NFFT parameters")
    call messages_new_line()
    call messages_write("      Fourier coefficients      N = ")
    do idir = 1,  nfft%dim
      call messages_write(nfft%N(idir))
      if(idir < nfft%dim) call messages_write(" x ")
    end do
    call messages_new_line()

    call messages_write("      Spatial nodes             M = ")

!     mm = nfft%M(1)*nfft%M(2)*nfft%M(3)
!
!     call messages_write(mm)
!     call messages_new_line()
    do idir = 1,  nfft%dim
      call messages_write(nfft%M(idir))
      if(idir < nfft%dim) call messages_write(" x ")
    end do
    call messages_new_line()


    call messages_write("      Oversampling factor   sigma = ")
    call messages_write(nfft%sigma)
    call messages_new_line()

    call messages_write("      FFT grid size             n = ")
    do idir = 1,  nfft%dim
      call messages_write(nfft%fftN(idir))
      if(idir < nfft%dim) call messages_write(" x ")
    end do
    call messages_new_line()

    call messages_write("      Real space cutoff         m = ")
    call messages_write(nfft%mm)
    call messages_new_line()

    call messages_write("      Pre-computation strategy    = ")
    select case(nfft%precompute)
    case(nfft_PRE_LIN_PSI)
      call messages_write(" NFFT_PRE_LIN_PSI")
    case(nfft_PRE_PSI)
      call messages_write(" NFFT_PRE_PSI")
    case(nfft_PRE_FULL_PSI)
      call messages_write(" NFFT_PRE_FULL_PSI")
    end select

    call messages_info()


    POP_SUB(nfft_write_info)
  end subroutine nfft_write_info


  ! ---------------------------------------------------------
  subroutine nfft_end(nfft)
    type(nfft_t), intent(inout) :: nfft

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



    FLOAT   :: x1_(1:nfft%M(1)), x2_(1:nfft%M(2)), x3_(1:nfft%M(3))
    FLOAT   :: length, cc, eps, dX1(1:nfft%M(1)-1),  dX2(1:nfft%M(2)-1), dX3(1:nfft%M(3)-1)

    integer :: ii

    PUSH_SUB(nfft_precompute)

!     eps = 1.000001 ! the sample nodes must be in [0.5,0.5)
    eps = M_ONE+M_EPSILON ! the sample nodes must be in [0.5,0.5)

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
        do ii = 1, nfft%M(1)-1
          dX1(ii)= abs(x1_(ii+1)-x1_(ii))
        end  do
        do ii = 1, nfft%M(2)-1
          dX2(ii)= abs(x2_(ii+1)-x2_(ii))
        end do
        do ii = 1, nfft%M(3)-1
          dX3(ii)= abs(x3_(ii+1)-x3_(ii))
        end do
        nfft%norm = M_ONE/(minval(dX1(:)) * minval(dX2(:)) * minval(dX3(:)))

      case(2)
        length = (maxval(X1)-minval(X1))*eps
        cc = (minval(X1)+maxval(X1))/M_TWO
        x1_ =(X1-cc)/length
        length = (maxval(X2)-minval(X2))*eps
        cc = (minval(X2)+maxval(X2))/M_TWO
        x2_ =(X2-cc)/length
        call oct_nfft_precompute_one_psi_2d(nfft%plan, nfft%M, x1_, x2_)

        ! Set the normalization factor
        do ii = 1, nfft%M(1)-1
          dX1(ii)= abs(x1_(ii+1)-x1_(ii))
        end do
        do ii = 1, nfft%M(2)-1
          dX2(ii)= abs(x2_(ii+1)-x2_(ii))
        end do
        nfft%norm = M_ONE/(minval(dX1(:)) * minval(dX2(:)))


      case(1)
        length = (maxval(X1)-minval(X1))*eps
        cc = (minval(X1)+maxval(X1))/M_TWO
        x1_ =(X1-cc)/length
        call oct_nfft_precompute_one_psi_1d(nfft%plan,nfft%M(1),x1_)

        ! Set the normalization factor
        do ii = 1, nfft%M(1)-1
          dX1(ii)= abs(x1_(ii+1)-x1_(ii))
        end do
        nfft%norm = M_ONE/(minval(dX1(:)))

    end select

    ! check the plan
    call oct_nfft_check(nfft%plan)


    write(message(1), '(a)') "Info: NFFT plan precomputed."
    call messages_info(1)


    POP_SUB(nfft_precompute)
  end subroutine nfft_precompute

  !--------------------------------------------
  subroutine dnfft_forward(nfft, in, out)
    type(nfft_t), intent(in)  :: nfft
    FLOAT,        intent(in)  :: in(:,:,:)
    CMPLX,        intent(out) :: out(:,:,:)

    CMPLX, allocatable :: zin(:,:,:)
    integer:: b(6)

    PUSH_SUB(dnfft_forward)

    b(1) = lbound(in, dim=1)
    b(2) = ubound(in, dim=1)
    b(3) = lbound(in, dim=2)
    b(4) = ubound(in, dim=3)
    b(5) = lbound(in, dim=3)
    b(6) = ubound(in, dim=3)

!    SAxFE_ALLOCATE(zin(b(1):b(2),b(3):b(4),b(5):b(6)))
    allocate(zin(b(1):b(2),b(3):b(4),b(5):b(6)))
    zin = in
    call znfft_forward(nfft, zin, out)

    deallocate(zin)
!     SAxFE_DEALLOCATE_A(zin)

    POP_SUB(dnfft_forward)
  end subroutine dnfft_forward

  !--------------------------------------------
  subroutine dnfft_backward(nfft, in, out)
    type(nfft_t), intent(in)  :: nfft
    CMPLX,        intent(in)  :: in (:,:,:)
    FLOAT,        intent(out) :: out(:,:,:)

    CMPLX, allocatable :: zout(:,:,:)
    integer:: b(6)

    PUSH_SUB(dnfft_backward)

    b(1) = lbound(out, dim=1)
    b(2) = ubound(out, dim=1)
    b(3) = lbound(out, dim=2)
    b(4) = ubound(out, dim=3)
    b(5) = lbound(out, dim=3)
    b(6) = ubound(out, dim=3)

    allocate(zout(b(1):b(2),b(3):b(4),b(5):b(6)))

    call znfft_backward(nfft, in, zout)
    out = zout
    deallocate(zout)

    POP_SUB(dnfft_backward)
  end subroutine dnfft_backward

  !--------------------------------------------
  subroutine znfft_forward(nfft, in, out)
    type(nfft_t), intent(in)  :: nfft
    CMPLX,        intent(in)  :: in(:,:,:)
    CMPLX,        intent(out) :: out(:,:,:)

    integer :: ix, iy, iz

    PUSH_SUB(znfft_forward)

    do ix = 1, nfft%N(1)
      do iy = 1, nfft%N(2)
        do iz = 1, nfft%N(3)
          call zoct_set_f_hat(nfft%plan, nfft%dim, in(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do

    call oct_nfft_trafo(nfft%plan)

    do ix = 1, nfft%M(1)
      do iy = 1, nfft%M(2)
        do iz = 1, nfft%M(3)
          call zoct_get_f(nfft%plan, nfft%M, nfft%dim, out(ix,iy,iz), ix, iy, iz)
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

    integer :: ix, iy, iz

    PUSH_SUB(znfft_backward)

    do ix = 1, nfft%M(1)
      do iy = 1, nfft%M(2)
        do iz = 1, nfft%M(3)
          call zoct_set_f(nfft%plan, nfft%M, nfft%dim, in(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do

    call oct_nfft_adjoint(nfft%plan)

    do ix = 1,nfft%N(1)
      do iy = 1, nfft%N(2)
        do iz = 1, nfft%N(3)
          call zoct_get_f_hat(nfft%plan, nfft%dim, out(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do

    out = out/nfft%norm

    POP_SUB(znfft_backward)
  end subroutine znfft_backward

end module nfft_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
