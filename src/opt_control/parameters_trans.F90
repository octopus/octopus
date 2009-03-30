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
!! $Id: opt_control.F90 2870 2007-04-28 06:26:47Z acastro $


  ! ---------------------------------------------------------
  subroutine parameters_par_to_x(par, xx)
    type(oct_control_parameters_t), intent(inout) :: par
    REAL_DOUBLE,                    intent(inout) :: xx(:)
    integer :: j, n
    FLOAT, allocatable :: ep(:), e(:), y(:), a(:), x(:)

    call push_sub('parameters.parameters_par_to_x')

    n = par%dim
    ALLOCATE(e(n), n)
    ALLOCATE(ep(n), n)
    ALLOCATE(x(size(xx)), size(xx))
    x = xx

    ASSERT(par%current_representation .ne. ctr_real_time)

    do j =  1, n
      ep(j) = tdf(par%f(1), j)
    end do
    e = matmul(par%utransf, ep)

    if(par_common%representation .eq. ctr_zero_fourier_series) then
      ALLOCATE(a(n-1), n-1)
      ALLOCATE(y(n-1), n-1)
      a = M_ZERO
      a(1:n/2-1) = M_ONE

      call hypersphere_cut(ep(2:n), a, y)
      call cartesian2hyperspherical(y, x(1:n-2))
      deallocate(a, y)
     else
      call cartesian2hyperspherical(e, x(1:n-1))
    end if
    
    xx = x
    deallocate(e, ep, x) 
    call pop_sub()
  end subroutine parameters_par_to_x
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_x_to_par(par, xx)
    type(oct_control_parameters_t), intent(inout) :: par
    REAL_DOUBLE,                          intent(in)    :: xx(:)

    FLOAT, allocatable :: y(:), a(:), e(:), ep(:), x(:)
    integer :: n

    call push_sub('parameters.parameters_x_to_par')

    ASSERT(par%current_representation .ne. ctr_real_time)

    n = par%dim
    ALLOCATE(e(n), n)
    ALLOCATE(ep(n), n)
    ALLOCATE(x(size(xx)), size(xx))
    x = xx

    if(par_common%representation .eq. ctr_zero_fourier_series) then
      call hyperspherical2cartesian(x, ep(2:n))
      ALLOCATE(a(n-1), n-1)
      ALLOCATE(y(n-1), n-1)
      a = M_ZERO
      a(1:n/2-1) = M_ONE
      call hypersphere_cut_back(ep(2:n), a, e(2:n))
      e(1) = -sum(e(2:n/2))
      e = sqrt(par_common%targetfluence) * e
      ep = matmul(par%utransfi, e)
      call tdf_set_numerical(par%f(1), ep)
      deallocate(a, y)
    else
      call hyperspherical2cartesian(x(1:n-1), e)
      e = sqrt(par_common%targetfluence) * e
      ep = matmul(par%utransfi, e)
      call tdf_set_numerical(par%f(1), ep)
    end if

    deallocate(e, ep, x)
    call pop_sub()
  end subroutine parameters_x_to_par
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_trans_matrix(par)
    type(oct_control_parameters_t), intent(inout) :: par

    integer :: i, mm, nn
    FLOAT :: t, det
    type(tdf_t) :: fn, fm
    FLOAT, allocatable :: eigenvec(:, :), eigenval(:)

    if(par_common%representation .eq. ctr_real_time) return

    call push_sub('parameters.parameters_trans_matrix')

    ALLOCATE(par%utransf (par%dim, par%dim), par%dim*par%dim)
    ALLOCATE(par%utransfi(par%dim, par%dim), par%dim*par%dim)
    par%utransf  = M_ZERO
    par%utransfi = M_ZERO
    forall(mm = 1:par%dim) par%utransf(mm, mm) = M_ONE
    forall(mm = 1:par%dim) par%utransfi(mm, mm) = M_ONE

    if( par_common%mode .eq. parameter_mode_f ) then
      ! If the object to optimize is the envelope of the
      ! the laser pulse. Being e(t) the laser pulse, it is assumed that it
      ! has the form:
      !   e(t) = f(t) cos(w0*t),
      ! where f(t) is the envelope. This is then expanded in a basis set:
      !   f(t) = sum_{n=1}^N f_n g_n(t).
      ! The fluence F[e] is then given by:
      !   F[e] = sum_{m=1}^N sum_{n=1}^N f_m f_n S_{nm},
      !   S_{nm} = \int_{0}^{T} dt g_n(t) g_m(t) cos^2(w0*t).
      ! The following lines of code calculate a matrix U, placed in par%utransf,
      ! that performs the transformation of variables that takes {f_n} to
      ! {h_n}, (\vec{h} = U\vec{f}), such that
      !   F[e] = sum-{m=1}^N h_n^2.
      ! The inverse matrix U^{-1} is placed in par%utransfi.
      !
      ! This scan probably be optimized in some way?
      ALLOCATE(eigenvec(par%dim, par%dim), par%dim*par%dim)
      ALLOCATE(eigenval(par%dim), par%dim)
      par%utransf  = M_ZERO
      par%utransfi = M_ZERO

      do mm = 1, par%dim
        call tdf_init_numerical(fm, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
          par_common%omegamax, rep = par_common%representation)
        call tdf_set_numerical(fm, mm, M_ONE)
        select case(par_common%representation)
          case(ctr_sine_fourier_series); call tdf_sineseries_to_numerical(fm)
          case(ctr_fourier_series);      call tdf_fourier_to_numerical(fm)
          case(ctr_zero_fourier_series); call tdf_zerofourier_to_numerical(fm)
        end select
        do i = 1, tdf_niter(fm) + 1
          t = (i-1)*tdf_dt(fm)
          call tdf_set_numerical(fm, i, tdf(fm, i)*cos(par%w0*t))
        end do

        do nn = mm, par%dim
          call tdf_init_numerical(fn, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
            par_common%omegamax, rep = par_common%representation)
          call tdf_set_numerical(fn, nn, M_ONE)
          select case(par_common%representation)
            case(ctr_sine_fourier_series); call tdf_sineseries_to_numerical(fn)
            case(ctr_fourier_series);      call tdf_fourier_to_numerical(fn)
            case(ctr_zero_fourier_series); call tdf_zerofourier_to_numerical(fn)
          end select
          do i = 1, tdf_niter(fn) + 1
            t = (i-1)*tdf_dt(fn)
            call tdf_set_numerical(fn, i, tdf(fn, i)*cos(par%w0*t))
          end do
          par%utransf(mm, nn) = tdf_dot_product(fm, fn)
          call tdf_end(fn)
        end do
        call tdf_end(fm)
      end do

      do mm = 1, par%dim
        do nn = 1, mm - 1
          par%utransf(mm, nn) = par%utransf(nn, mm)
        end do
      end do

      call lalg_eigensolve(par%dim, par%utransf, eigenvec, eigenval)
      do mm = 1, par%dim
        do nn = 1, par%dim
          eigenvec(mm, nn) = eigenvec(mm, nn) * sqrt(eigenval(nn))
        end do
      end do

      par%utransf = transpose(eigenvec)
      par%utransfi = par%utransf
      det =  lalg_inverter(par%dim, par%utransfi)

    end if

    call pop_sub()
  end subroutine parameters_trans_matrix
  ! ---------------------------------------------------------
