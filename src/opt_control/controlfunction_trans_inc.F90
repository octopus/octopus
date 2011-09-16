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
  subroutine controlfunction_basis_to_theta(par)
    type(controlfunction_t), intent(inout) :: par
    integer :: j, n, dof
    FLOAT, allocatable :: ep(:), e(:), y(:), a(:), x(:)

    PUSH_SUB(controlfunction_basis_to_theta)

    ASSERT(par%current_representation .ne. ctr_real_time)

    select case(par%current_representation)
    case(ctr_fourier_series_h, ctr_zero_fourier_series_h)
      n = par%dim
      dof = par%dof
      SAFE_ALLOCATE( e(1:n))
      SAFE_ALLOCATE(ep(1:n))
      SAFE_ALLOCATE(x(1:dof))

      forall(j = 1: n) ep(j) = tdf(par%f(1), j)
      e = matmul(par%utransf, ep)

      if(cf_common%representation .eq. ctr_zero_fourier_series_h) then
        SAFE_ALLOCATE(y(1:n-1))
        y = matmul(par%hypersphere_transform, ep(2:n))
        call cartesian2hyperspherical(y, x(1:n-2))
        SAFE_DEALLOCATE_A(y)
       else
        call cartesian2hyperspherical(e, x(1:n-1))
      end if

      par%theta = x
      SAFE_DEALLOCATE_A(e)
      SAFE_DEALLOCATE_A(ep)
      SAFE_DEALLOCATE_A(x)

    case(ctr_fourier_series)
      forall(j = 1: par%dim) par%theta(j) = tdf(par%f(1), j)

    case(ctr_zero_fourier_series)
      ! In this case, the transformation is (n = par%dim):
      ! theta(1) = a(2)
      ! theta(2) = a(3)
      ! ...      = ...
      ! theta(n/2-1)   = a(n/2)
      ! theta(n/2) = b(1)
      ! ...      = ...
      ! theta(n-1) = b(n/2)
      ! where a are the coefficients of the cosines in the Fourier series, and b are
      ! the coefficients of the sines

      SAFE_ALLOCATE(e(1:par%dim))
      forall(j = 1: par%dim) e(j) = tdf(par%f(1), j)

      do j = 2, par%dim
        par%theta(j-1) = e(j)
      end do

      SAFE_DEALLOCATE_A(e)
      
    end select

    POP_SUB(controlfunction_basis_to_theta)
  end subroutine controlfunction_basis_to_theta
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_theta_to_basis(par)
    type(controlfunction_t), intent(inout) :: par

    FLOAT, allocatable :: y(:), a(:), e(:), ep(:), x(:)
    integer :: n, dof, j

    PUSH_SUB(controlfunction_theta_to_basis)

    ASSERT(par%current_representation .ne. ctr_real_time)


    select case(par%current_representation)
    case(ctr_fourier_series_h, ctr_zero_fourier_series_h)

      n = par%dim
      dof = par%dof
      SAFE_ALLOCATE( e(1:n))
      SAFE_ALLOCATE(ep(1:n))
      SAFE_ALLOCATE(x(1:dof))
      x = par%theta

      if(cf_common%representation .eq. ctr_zero_fourier_series_h) then
        call hyperspherical2cartesian(x, ep(2:n))
        e(2:n) = matmul(par%ihypersphere_transform, ep(2:n))
        e(1) = -sum(e(2:n/2))
        e = sqrt(cf_common%targetfluence) * e
        ep = matmul(par%utransfi, e)
        call tdf_set_numerical(par%f(1), ep)
      else
        call hyperspherical2cartesian(x(1:n-1), e)
        e = sqrt(cf_common%targetfluence) * e
        ep = matmul(par%utransfi, e)
        call tdf_set_numerical(par%f(1), ep)
      end if

      SAFE_DEALLOCATE_A(e)
      SAFE_DEALLOCATE_A(ep)
      SAFE_DEALLOCATE_A(x)

    case(ctr_fourier_series)

      call tdf_set_numerical(par%f(1), par%theta)

    case(ctr_zero_fourier_series)

      ! In this case, the transformation is (n = par%dim):
      ! a(1) = -sum(a(2)...a(n/2))  represents the constraint
      ! theta(1) = a(2)
      ! theta(2) = a(3)
      ! ...      = ...
      ! theta(n/2-1)   = a(n/2)
      ! theta(n/2) = b(1)
      ! ...      = ...
      ! theta(n-1) = b(n/2)
      ! where a are the coefficients of the cosines in the Fourier series, and b are
      ! the coefficients of the sines (Note that in the next lines "a" are the coefficients
      ! of both sines and cosines.

      n = par%dim
      SAFE_ALLOCATE(a(1:n))

      do j = 2, n
        a(j) = par%theta(j-1)
      end do
       a(1) = -sum(a(2:n/2))
      call tdf_set_numerical(par%f(1), a)

      SAFE_DEALLOCATE_A(a)

    end select

    POP_SUB(controlfunction_theta_to_basis)
  end subroutine controlfunction_theta_to_basis
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_get_theta(par, theta)
    type(controlfunction_t), intent(in) :: par
    FLOAT, intent(inout) :: theta(:)

    PUSH_SUB(controlfunction_get_theta)
    theta = par%theta

    POP_SUB(controlfunction_get_theta)
  end subroutine controlfunction_get_theta
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_set_theta(par, theta)
    type(controlfunction_t), intent(inout) :: par
    FLOAT, intent(in) :: theta(:)

    PUSH_SUB(controlfunction_set_theta)
    par%theta = theta

    POP_SUB(controlfunction_set_theta)
  end subroutine controlfunction_set_theta
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_trans_matrix(par)
    type(controlfunction_t), intent(inout) :: par

    integer :: i, mm, nn, n, j, k
    FLOAT :: t, det, w1
    type(tdf_t) :: fn, fm
    FLOAT, allocatable :: neigenvec(:, :), eigenvec(:, :), eigenval(:), a(:), alpha(:, :)

    PUSH_SUB(controlfunction_trans_matrix)

    select case(cf_common%representation)

    case(ctr_real_time) 

      POP_SUB(controlfunction_trans_matrix)
      return

    case(ctr_fourier_series)

      SAFE_ALLOCATE(par%utransf (1:par%dim, 1:par%dim))
      SAFE_ALLOCATE(par%utransfi(1:par%dim, 1:par%dim))
      par%utransf  = M_ZERO
      par%utransfi = M_ZERO
      forall(mm = 1:par%dim) par%utransf(mm, mm) = M_ONE
      forall(mm = 1:par%dim) par%utransfi(mm, mm) = M_ONE

      if( cf_common%mode .eq. controlfunction_mode_f ) then
        ! If the object to optimize is the envelope of the
        ! the laser pulse. Being e(t) the laser pulse, it is assumed that it
        ! has the form:
        !   e(t) = f(t) cos(w0*t),
        ! where f(t) is the envelope. This is then expanded in a basis set:
        !   f(t) = sum_{n=1}^N f_n g_n(t).
        ! The fluence F[e] is then given by:
        !   F[e] = sum_{m=1}^N sum_{n=1}^N f_m f_n S_{nm},
        !   S_{nm} = \int_{0}^{T} dt g_n(t) g_m(t) cos^2(w0*t).
        ! The following lines of code calculate this matrix S and place it in par%utransf,
        ! This can probably be optimized in some way?
        par%utransf  = M_ZERO
        par%utransfi = M_ZERO

        do mm = 1, par%dim
          call tdf_init_numerical(fm, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
            cf_common%omegamax, rep = TDF_FOURIER_SERIES)
          call tdf_set_numerical(fm, mm, M_ONE)
          call tdf_fourier_to_numerical(fm)
          do i = 1, tdf_niter(fm) + 1
            t = (i-1)*tdf_dt(fm)
            call tdf_set_numerical(fm, i, tdf(fm, i)*cos(par%w0*t))
          end do

          do nn = mm, par%dim
            call tdf_init_numerical(fn, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
              cf_common%omegamax, rep = TDF_FOURIER_SERIES)
            call tdf_set_numerical(fn, nn, M_ONE)
            call tdf_fourier_to_numerical(fn)
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

      end if

    case(ctr_zero_fourier_series)

      SAFE_ALLOCATE(par%utransf (1:par%dof, 1:par%dof))
      SAFE_ALLOCATE(par%utransfi(1:par%dof, 1:par%dof))

      par%utransf  = M_ZERO
      par%utransfi = M_ZERO

      w1 = (M_TWO*M_PI/(tdf_dt(par%f(1))*tdf_niter(par%f(1))))

      do mm = 1, par%dof
        call tdf_init_numerical(fm, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
          cf_common%omegamax, rep = TDF_ZERO_FOURIER)
        call tdf_set_numerical(fm, mm+1, M_ONE)
        call tdf_fourier_to_numerical(fm)
        if(mm <= par%dof/2) then
          do i = 1, tdf_niter(fm) + 1
            t = (i-1)*tdf_dt(fm)
            call tdf_set_numerical(fm, i, &
              (tdf(fm, i)-sqrt(M_TWO/(tdf_dt(fm)*tdf_niter(fm)))*cos(w1*t))*cos(par%w0*t))
          end do
        else
          do i = 1, tdf_niter(fm) + 1
            t = (i-1)*tdf_dt(fm)
            call tdf_set_numerical(fm, i, tdf(fm, i)*cos(par%w0*t))
          end do
        end if

        do nn = mm, par%dof
          call tdf_init_numerical(fn, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
            cf_common%omegamax, rep = TDF_ZERO_FOURIER)
          call tdf_set_numerical(fn, nn+1, M_ONE)
          call tdf_fourier_to_numerical(fn)
          if(nn <= par%dof/2) then
            do i = 1, tdf_niter(fn) + 1
              t = (i-1)*tdf_dt(fn)
              call tdf_set_numerical(fn, i, &
                (tdf(fn, i)-sqrt(M_TWO/(tdf_dt(fm)*tdf_niter(fm)))*cos(w1*t))*cos(par%w0*t))
            end do
          else
            do i = 1, tdf_niter(fn) + 1
              t = (i-1)*tdf_dt(fn)         
              call tdf_set_numerical(fn, i, tdf(fn, i)*cos(par%w0*t))
            end do
          end if
          par%utransf(mm, nn) = tdf_dot_product(fm, fn)
          call tdf_end(fn)
        end do
        call tdf_end(fm)
      end do

      do mm = 1, par%dof
        do nn = 1, mm - 1
          par%utransf(mm, nn) = par%utransf(nn, mm)
        end do
      end do

    case(ctr_zero_fourier_series_h)

      n = par%dim
      SAFE_ALLOCATE(par%hypersphere_transform (n-1, 1:n-1))
      SAFE_ALLOCATE(par%ihypersphere_transform (n-1, 1:n-1))
      SAFE_ALLOCATE(alpha(1:n-1, 1:n-1))
      SAFE_ALLOCATE(eigenvec(1:n-1, 1:n-1))
      SAFE_ALLOCATE(neigenvec(1:n-1, 1:n-1))
      SAFE_ALLOCATE(eigenval(1:n-1))
      ! k = 1
      eigenvec = M_ZERO
      eigenvec(1:n/2-1, 1) = M_ONE
      eigenvec(n/2:n-1, 1) = M_ZERO
      eigenval(1) = n/2
      ! k = 2, ...., n/2-1
      do k = 2, n/2-1
        eigenval(k) = M_ONE
        eigenvec(:, k) = M_ZERO
        eigenvec(1, k) = -M_ONE
        eigenvec(k, k) = M_ONE
      end do
      ! k = n/2, ..., n-1
      do k = n/2, n-1
        eigenval(k) = M_ONE
        eigenvec(:, k) = M_ZERO
        eigenvec(k, k) = M_ONE
      end do

      do k = 1, n-1
        neigenvec(:, k) = eigenvec(:, k)
        do j = 1, k-1
          neigenvec(:, k) = neigenvec(:, k) - & 
            (dot_product(eigenvec(:, k), neigenvec(:, j)) / dot_product(neigenvec(:, j), neigenvec(:, j))) * neigenvec(:, j)
        end do
      end do
      do k = 1, n-1
        neigenvec(:, k) = neigenvec(:, k)/sqrt(dot_product(neigenvec(:, k),neigenvec(:, k)))
      end do
      eigenvec = neigenvec

      par%ihypersphere_transform = eigenvec
      par%hypersphere_transform = transpose(eigenvec)

      forall(j=1:n-1, k=1:n-1) 
        par%ihypersphere_transform(j, k) = par%ihypersphere_transform(j, k) / sqrt(eigenval(k))
        par%hypersphere_transform(j, k) = par%hypersphere_transform(j, k) * sqrt(eigenval(k))
      end forall

      SAFE_DEALLOCATE_A(alpha)
      SAFE_DEALLOCATE_A(eigenvec)
      SAFE_DEALLOCATE_A(neigenvec)
      SAFE_DEALLOCATE_A(eigenval)

      SAFE_ALLOCATE(par%utransf (1:par%dim, 1:par%dim))
      SAFE_ALLOCATE(par%utransfi(1:par%dim, 1:par%dim))
      par%utransf  = M_ZERO
      par%utransfi = M_ZERO
      forall(mm = 1:par%dim) par%utransf(mm, mm) = M_ONE
      forall(mm = 1:par%dim) par%utransfi(mm, mm) = M_ONE

      ! WARNING: Here we are missing the cases in which cf_common%representation is
      ! either ctr_fourier_series or ctr_zero_fourier_series.
      if( cf_common%mode .eq. controlfunction_mode_f ) then
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
        SAFE_ALLOCATE(eigenvec(1:par%dim, 1:par%dim))
        SAFE_ALLOCATE(eigenval(1:par%dim))
        par%utransf  = M_ZERO
        par%utransfi = M_ZERO

        do mm = 1, par%dim
          select case(cf_common%representation)
          case(ctr_fourier_series, ctr_fourier_series_h)
            call tdf_init_numerical(fm, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
              cf_common%omegamax, rep = TDF_FOURIER_SERIES)
          case(ctr_zero_fourier_series_h, ctr_zero_fourier_series)
            call tdf_init_numerical(fm, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
              cf_common%omegamax, rep = TDF_ZERO_FOURIER)
          end select
          call tdf_set_numerical(fm, mm, M_ONE)
          select case(cf_common%representation)
            case(ctr_fourier_series_h, ctr_fourier_series)
              call tdf_fourier_to_numerical(fm)
            case(ctr_zero_fourier_series_h, ctr_zero_fourier_series)
              call tdf_zerofourier_to_numerical(fm)
          end select
          do i = 1, tdf_niter(fm) + 1
            t = (i-1)*tdf_dt(fm)
            call tdf_set_numerical(fm, i, tdf(fm, i)*cos(par%w0*t))
          end do

          do nn = mm, par%dim
            select case(cf_common%representation)
            case(ctr_fourier_series_h, ctr_fourier_series)
              call tdf_init_numerical(fn, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
                cf_common%omegamax, rep = TDF_FOURIER_SERIES)
            case(ctr_zero_fourier_series_h, ctr_zero_fourier_series)
              call tdf_init_numerical(fn, tdf_niter(par%f(1)), tdf_dt(par%f(1)), &
                cf_common%omegamax, rep = TDF_ZERO_FOURIER)
            end select
            call tdf_set_numerical(fn, nn, M_ONE)
            select case(cf_common%representation)
              case(ctr_fourier_series_h, ctr_fourier_series)
                call tdf_fourier_to_numerical(fn)
              case(ctr_zero_fourier_series_h, ctr_zero_fourier_series)
                call tdf_zerofourier_to_numerical(fn)
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

        eigenvec = par%utransf
        call lalg_eigensolve(par%dim, eigenvec, eigenval)

        ! We need to make sure that eigenvectors have the same sign on all machines, which is not guaranteed
        ! by LAPACK. So, we will use the following criterion: the sign of the first non-null component should be
        ! positive.
        do nn = 1, par%dim
          do mm = 1, par%dim
            if( eigenvec(mm, nn)*eigenvec(mm, nn) > CNST(1.0e-20) ) then
              !eigenvec(1:par%dim, nn) = sign(eigenvec(mm, nn), M_ONE) * eigenvec(1:par%dim, nn)
              if(eigenvec(mm, nn) < M_ZERO) eigenvec(1:par%dim, nn) = - eigenvec(1:par%dim, nn)
              exit
            end if
          end do
        end do

        do mm = 1, par%dim
          do nn = 1, par%dim
            eigenvec(mm, nn) = eigenvec(mm, nn) * sqrt(eigenval(nn))
          end do
        end do
        par%utransf = transpose(eigenvec)
        par%utransfi = par%utransf

      end if

    case default

      SAFE_ALLOCATE(par%utransf (1:par%dim, 1:par%dim))
      SAFE_ALLOCATE(par%utransfi(1:par%dim, 1:par%dim))
      par%utransf  = M_ZERO
      par%utransfi = M_ZERO
      forall(mm = 1:par%dim) par%utransf(mm, mm) = M_ONE
      forall(mm = 1:par%dim) par%utransfi(mm, mm) = M_ONE

    end select

    POP_SUB(controlfunction_trans_matrix)
  end subroutine controlfunction_trans_matrix
  ! ---------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
