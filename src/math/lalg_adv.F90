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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module lalg_adv_m
  use global_m
  use lapack_m
  use math_m
  use messages_m
  use mpi_m
  use profiling_m
  use blacs_proc_grid_m
  use scalapack_m
  use utils_m
  
  implicit none

  private
  public ::                       &
    lalg_cholesky,                &
    lalg_geneigensolve,           &
    lalg_eigensolve,              &
    lalg_eigensolve_nonh,         &
    lalg_determinant,             &
    lalg_inverter,                &
    lalg_sym_inverter,            &
    lalg_linsyssolve,             &
    lalg_singular_value_decomp,   &
    lalg_svd_inverse,             &
    lalg_invert_upper_triangular, &
    lalg_invert_lower_triangular, &
    lalg_lowest_geneigensolve,    &
    lalg_lowest_eigensolve,       &
    zlalg_exp,                    &
    zlalg_phi,                    &
    lalg_zpseudoinverse,          &
    lalg_zeigenderivatives,       &
    lalg_check_zeigenderivatives, &
    lalg_zdni,                    &
    lalg_zduialpha,               &
    lalg_zd2ni

  type(profile_t), save :: cholesky_prof, eigensolver_prof

  interface lalg_cholesky
    module procedure dcholesky, zcholesky
  end interface lalg_cholesky

  interface lalg_geneigensolve
    module procedure dgeneigensolve, zgeneigensolve
  end interface lalg_geneigensolve

  interface lalg_eigensolve_nonh
    module procedure zeigensolve_nonh, deigensolve_nonh
  end interface lalg_eigensolve_nonh

  interface lalg_eigensolve
    module procedure deigensolve, zeigensolve
#ifdef HAVE_SCALAPACK
    module procedure deigensolve_scalapack, zeigensolve_scalapack
#endif
  end interface lalg_eigensolve

  !> Note that lalg_determinant and lalg_inverter are just wrappers
  !! over the same routine.
  interface lalg_determinant
    module procedure ddeterminant, zdeterminant
  end interface lalg_determinant

  interface lalg_inverter
    module procedure ddeterminant, zdeterminant
  end interface lalg_inverter

  interface lalg_sym_inverter
    module procedure dsym_inverter, zsym_inverter
  end interface lalg_sym_inverter

  interface lalg_linsyssolve
    module procedure dlinsyssolve, zlinsyssolve
  end interface lalg_linsyssolve

  interface lalg_singular_value_decomp
    module procedure zsingular_value_decomp
  end interface lalg_singular_value_decomp

  interface lalg_svd_inverse
    module procedure zsvd_inverse
  end interface lalg_svd_inverse

  interface lalg_invert_upper_triangular
    module procedure dinvert_upper_triangular, zinvert_upper_triangular
  end interface lalg_invert_upper_triangular
  
  interface lalg_invert_lower_triangular
    module procedure dinvert_lower_triangular, zinvert_lower_triangular
  end interface lalg_invert_lower_triangular
  
  interface lalg_lowest_geneigensolve
    module procedure dlowest_geneigensolve, zlowest_geneigensolve
  end interface lalg_lowest_geneigensolve

  interface lalg_lowest_eigensolve
    module procedure dlowest_eigensolve, zlowest_eigensolve
  end interface lalg_lowest_eigensolve

contains

  !>-------------------------------------------------
  !!
  !! This routine calculates the exponential of a matrix by using an
  !! eigenvalue decomposition.
  !!
  !! For the hermitian case:
  !!
  !!   A = V D V^T => exp(A) = V exp(D) V^T
  !!
  !! and in general
  !!
  !!   A = V D V^-1 => exp(A) = V exp(D) V^-1
  !!
  !! This is slow but it is simple to implement, and for the moment it
  !! does not affect performance.
  !!
  !!---------------------------------------------
  subroutine zlalg_exp(nn, pp, aa, ex, hermitian)
    integer,           intent(in)      :: nn
    CMPLX,             intent(in)      :: pp
    CMPLX,             intent(in)      :: aa(:, :)
    CMPLX,             intent(inout)   :: ex(:, :)
    logical,           intent(in)      :: hermitian

    CMPLX, allocatable :: evectors(:, :), zevalues(:)
    FLOAT, allocatable :: evalues(:)
    CMPLX :: deter
    
    integer :: ii

    PUSH_SUB(zlalg_exp)

    SAFE_ALLOCATE(evectors(1:nn, 1:nn))

    if(hermitian) then
      SAFE_ALLOCATE(evalues(1:nn))
      SAFE_ALLOCATE(zevalues(1:nn))

      evectors(1:nn, 1:nn) = aa(1:nn, 1:nn)

      call lalg_eigensolve(nn, evectors, evalues)

      zevalues(1:nn) = exp(pp*evalues(1:nn))

      do ii = 1, nn
        ex(1:nn, ii) = zevalues(1:nn)*conjg(evectors(ii, 1:nn))
      end do

      ex(1:nn, 1:nn) = matmul(evectors(1:nn, 1:nn), ex(1:nn, 1:nn))

      SAFE_DEALLOCATE_A(evalues)
      SAFE_DEALLOCATE_A(zevalues)
    else
      SAFE_ALLOCATE(zevalues(1:nn))

      evectors(1:nn, 1:nn) = aa(1:nn, 1:nn)

      call lalg_eigensolve_nonh(nn, evectors, zevalues)
      
      zevalues(1:nn) = exp(pp*zevalues(1:nn))

      ex(1:nn, 1:nn) = evectors(1:nn, 1:nn)

      deter = lalg_inverter(nn, evectors)
      
      do ii = 1, nn
        evectors(1:nn, ii) = zevalues(1:nn)*evectors(1:nn, ii)
      end do
      
      ex(1:nn, 1:nn) = matmul(ex(1:nn, 1:nn), evectors(1:nn, 1:nn))

      SAFE_DEALLOCATE_A(zevalues)
    end if

    POP_SUB(zlalg_exp)
  end subroutine zlalg_exp


  !>-------------------------------------------------
  !!
  !! This routine calculates phi(pp*A), where A is a matrix,
  !! pp is any complex number, and phi is the function:
  !! 
  !! phi(x) = (e^x - 1)/x
  !!
  !! For the Hermitian case, for any function f:
  !!
  !!   A = V D V^T => f(A) = V f(D) V^T
  !!
  !! and in general
  !!
  !!   A = V D V^-1 => f(A) = V f(D) V^-1
  !!
  !!---------------------------------------------
  subroutine zlalg_phi(nn, pp, aa, ex, hermitian)
    integer,           intent(in)      :: nn
    CMPLX,             intent(in)      :: pp
    CMPLX,             intent(in)      :: aa(:, :)
    CMPLX,             intent(inout)   :: ex(:, :)
    logical,           intent(in)      :: hermitian

    CMPLX, allocatable :: evectors(:, :), zevalues(:)
    FLOAT, allocatable :: evalues(:)
    CMPLX :: deter
    
    integer :: ii

    PUSH_SUB(zlalg_phi)

    SAFE_ALLOCATE(evectors(1:nn, 1:nn))

    if(hermitian) then
      SAFE_ALLOCATE(evalues(1:nn))
      SAFE_ALLOCATE(zevalues(1:nn))

      evectors = aa

      call lalg_eigensolve(nn, evectors, evalues)

      forall(ii = 1:nn) zevalues(ii) = (exp(pp*evalues(ii)) - M_z1) / (pp*evalues(ii))

      do ii = 1, nn
        ex(1:nn, ii) = zevalues(1:nn)*conjg(evectors(ii, 1:nn))
      end do

      ex(1:nn, 1:nn) = matmul(evectors(1:nn, 1:nn), ex(1:nn, 1:nn))

      SAFE_DEALLOCATE_A(evalues)
      SAFE_DEALLOCATE_A(zevalues)
    else
      SAFE_ALLOCATE(zevalues(1:nn))

      evectors(1:nn, 1:nn) = aa(1:nn, 1:nn)

      call lalg_eigensolve_nonh(nn, evectors, zevalues)
      
      forall(ii = 1:nn) zevalues(ii) = (exp(pp*zevalues(ii)) - M_z1) / (pp*zevalues(ii))

      ex(1:nn, 1:nn) = evectors(1:nn, 1:nn)

      deter = lalg_inverter(nn, evectors)
      
      do ii = 1, nn
        evectors(1:nn, ii) = zevalues(1:nn)*evectors(1:nn, ii)
      end do
      
      ex(1:nn, 1:nn) = matmul(ex(1:nn, 1:nn), evectors(1:nn, 1:nn))

      SAFE_DEALLOCATE_A(zevalues)
    end if

    POP_SUB(zlalg_phi)
  end subroutine zlalg_phi


  !>-------------------------------------------------
  !! Computes the necessary ingredients to obtain,
  !! later, the first and second derivatives of the
  !! eigenvalues of a Hermitean complex matrix zmat,
  !! and the first derivatives of the eigenvectors.
  !!
  !! This follows the scheme of J. R. Magnus,
  !! Econometric Theory 1, 179 (1985), restricted to
  !! Hermitean matrices, although probably this can be
  !! found in other sources. 
  !!---------------------------------------------
  subroutine lalg_zeigenderivatives(n, mat, zeigenvec, zeigenval, zmat)
    integer, intent(in) :: n
    CMPLX, intent(in) :: mat(:, :)
    CMPLX, intent(out) :: zeigenvec(:, :)
    CMPLX, intent(out) :: zeigenval(:)
    CMPLX, intent(out) :: zmat(:, :, :)

    integer :: i, alpha, beta
    CMPLX, allocatable :: knaught(:, :)
    CMPLX, allocatable :: lambdaminusdm(:, :)
    CMPLX, allocatable :: ilambdaminusdm(:, :)
    CMPLX, allocatable :: unit(:, :)

    PUSH_SUB(lalg_zeigenderivatives)

    SAFE_ALLOCATE(unit(n, n))
    SAFE_ALLOCATE(knaught(n, n))
    SAFE_ALLOCATE(lambdaminusdm(n, n))
    SAFE_ALLOCATE(ilambdaminusdm(n, n))

    zeigenvec = mat
    call lalg_eigensolve_nonh(n, zeigenvec, zeigenval)

    unit = M_z0
    forall(i = 1:n) unit(i, i) = M_z1

    do i = 1, n
      do alpha = 1, n
        do beta = 1, n
          knaught(alpha, beta) = - zeigenvec(alpha, i)*conjg(zeigenvec(beta, i))
        end do
      end do
      knaught = knaught + unit
      lambdaminusdm = zeigenval(i)*unit - mat
      call lalg_zpseudoinverse(n, lambdaminusdm, ilambdaminusdm)
      zmat(:, :, i) = matmul(ilambdaminusdm, knaught)
    end do

    SAFE_DEALLOCATE_A(unit)
    SAFE_DEALLOCATE_A(knaught)
    SAFE_DEALLOCATE_A(lambdaminusdm)
    SAFE_DEALLOCATE_A(ilambdaminusdm)
    POP_SUB(lalg_zeigenderivatives)
  end subroutine lalg_zeigenderivatives


  !>-------------------------------------------------
  !! Computes the Moore-Penrose pseudoinverse of a
  !! complex matrix.
  !!---------------------------------------------
  subroutine lalg_zpseudoinverse(n, mat, imat)
    integer, intent(in) :: n
    CMPLX, intent(in) :: mat(:, :)
    CMPLX, intent(out) :: imat(:, :)

    integer :: i
    CMPLX, allocatable :: u(:, :), vt(:, :), sigma(:, :)
    FLOAT, allocatable :: sg_values(:)

    PUSH_SUB(lalg_zpseudoinverse)

    SAFE_ALLOCATE(u(n, n)) 
    SAFE_ALLOCATE(vt(n, n)) 
    SAFE_ALLOCATE(sigma(n, n)) 
    SAFE_ALLOCATE(sg_values(n))

    imat = mat
    call  lalg_singular_value_decomp(n, imat, u, vt, sg_values)

    sigma = M_z0
    do i = 1, n
      if(abs(sg_values(i)) <= CNST(1.0e-12) * maxval(abs(sg_values))) then
        sigma(i, i) = M_z0
      else
        sigma(i, i) = M_z1 / sg_values(i)
      end if
    end do

    vt = conjg(transpose(vt))
    u = conjg(transpose(u))
    imat = matmul(vt, matmul(sigma, u))

    ! Check if we truly have a pseudoinverse
    vt = matmul(mat, matmul(imat, mat)) - mat
    if( maxval(abs(vt)) > CNST(1.0e-10) * maxval(abs(mat))) then
      write(*, *) maxval(abs(vt))
      write(*, *) vt
      write(*, *)
      write(*, *) 1.0e-10 * maxval(abs(mat))
      write(*, *) maxval(abs(vt)) > CNST(1.0e-10) * maxval(abs(mat))
      write(*, *) mat
      write(message(1), '(a)') 'Pseudoinverse failed.'
      call messages_fatal(1)
    end if

    SAFE_DEALLOCATE_A(u) 
    SAFE_DEALLOCATE_A(vt) 
    SAFE_DEALLOCATE_A(sigma)
    SAFE_DEALLOCATE_A(sg_values)
    POP_SUB(lalg_zpseudoinverse)
  end subroutine lalg_zpseudoinverse


  !>-------------------------------------------------
  !! The purpose of this routine is to check that "lalg_zeigenderivatives"
  !! is working properly, and therefore, it is not really called anywhere
  !! in the code. It is here only for debugging purposes (perhaps it will
  !! disappear in the future...)
  !!-------------------------------------------------
  subroutine lalg_check_zeigenderivatives(n, mat)
    integer, intent(in) :: n
    CMPLX, intent(in) :: mat(:, :)

    integer :: alpha, beta, gamma, delta
    CMPLX :: deltah, zuder_direct, zder_direct, zuder_directplus, zuder_directminus
    CMPLX, allocatable :: zeigenvec(:, :), dm(:, :), zeigref_(:, :), zeigenval(:), mmatrix(:, :, :)
    CMPLX, allocatable :: zeigplus(:), zeigminus(:), zeig0(:), zeigplusminus(:), zeigminusplus(:)

    PUSH_SUB(lalg_check_zeigenderivatives)

    SAFE_ALLOCATE(zeigenvec(n, n))
    SAFE_ALLOCATE(dm(n, n))
    SAFE_ALLOCATE(zeigref_(n, n))
    SAFE_ALLOCATE(zeigenval(n))
    SAFE_ALLOCATE(mmatrix(n, n, n))
    SAFE_ALLOCATE(zeigplus(n))
    SAFE_ALLOCATE(zeigminus(n))
    SAFE_ALLOCATE(zeig0(n))
    SAFE_ALLOCATE(zeigplusminus(n))
    SAFE_ALLOCATE(zeigminusplus(n))

    ASSERT(n  ==  2)

    dm = mat
    call lalg_zeigenderivatives(2, dm, zeigref_, zeigenval, mmatrix)


    deltah = (CNST(0.000001), CNST(0.000))
    !deltah = M_z1 * maxval(abs(dm)) * CNST(0.001)
    do alpha = 1, n
      do beta = 1, n
        zder_direct = lalg_zdni(zeigref_(:, 2), alpha, beta)
        zuder_direct = lalg_zduialpha(zeigref_(:, 2), mmatrix(:, :, 2), 2, alpha, beta)

        zeigenvec = dm
        zeigenvec(alpha, beta) = zeigenvec(alpha, beta) + deltah
        call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigplus)
        zeigenvec(:, 1) = (M_ONE/sum(conjg(zeigref_(1:2, 1))*zeigenvec(1:2, 1))) * zeigenvec(:, 1)
        zeigenvec(:, 2) = (M_ONE/sum(conjg(zeigref_(1:2, 2))*zeigenvec(1:2, 2))) * zeigenvec(:, 2)
        zuder_directplus = zeigenvec(2, 2)

        zeigenvec = dm
        zeigenvec(alpha, beta) = zeigenvec(alpha, beta) - deltah
        call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigminus)
        zeigenvec(:, 1) = (M_ONE/sum(conjg(zeigref_(1:2, 1))*zeigenvec(1:2, 1))) * zeigenvec(:, 1)
        zeigenvec(:, 2) = (M_ONE/sum(conjg(zeigref_(1:2, 2))*zeigenvec(1:2, 2))) * zeigenvec(:, 2)
        zuder_directminus = zeigenvec(2, 2)

        write(*, '(2i1,4f24.12)') alpha, beta, zder_direct, (zeigplus(2) - zeigminus(2))/(M_TWO * deltah)
        write(*, '(2i1,4f24.12)') alpha, beta, &
            zuder_direct, (zuder_directplus - zuder_directminus) / (M_TWO * deltah)

        do gamma = 1, n
          do delta = 1, n
            if(alpha == gamma .and. beta == delta) then
               zder_direct = lalg_zd2ni(zeigref_(:, 1), mmatrix(:, :, 1), alpha, beta, gamma, delta)

               zeigenvec = dm
               call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeig0)

               zeigenvec = dm
               zeigenvec(alpha, beta) = zeigenvec(alpha, beta) + deltah
               call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigplus)

               zeigenvec = dm
               zeigenvec(alpha, beta) = zeigenvec(alpha, beta) - deltah
               call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigminus)

               write(*, '(4i1,4f24.12)') alpha, beta, gamma, delta, &
                 zder_direct, (zeigplus(1) + zeigminus(1) - M_TWO*zeig0(1))/(deltah**2)
            else
               zder_direct = lalg_zd2ni(zeigref_(:, 1), mmatrix(:, :, 1), alpha, beta, gamma, delta)

               zeigenvec = dm
               zeigenvec(alpha, beta) = zeigenvec(alpha, beta) + deltah
               zeigenvec(gamma, delta) = zeigenvec(gamma, delta) + deltah
               call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigplus)

               zeigenvec = dm
               zeigenvec(alpha, beta) = zeigenvec(alpha, beta) - deltah
               zeigenvec(gamma, delta) = zeigenvec(gamma, delta) - deltah
               call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigminus)

               zeigenvec = dm
               zeigenvec(alpha, beta) = zeigenvec(alpha, beta) + deltah
               zeigenvec(gamma, delta) = zeigenvec(gamma, delta) - deltah
               call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigplusminus)

               zeigenvec = dm
               zeigenvec(alpha, beta) = zeigenvec(alpha, beta) - deltah
               zeigenvec(gamma, delta) = zeigenvec(gamma, delta) + deltah
               call lalg_eigensolve_nonh(2, zeigenvec(:, :), zeigminusplus)

               write(*, '(4i1,4f24.12)') alpha, beta, gamma, delta, &
                 zder_direct, (zeigplus(1) + zeigminus(1) - zeigplusminus(1) - zeigminusplus(1))/(M_FOUR*deltah**2)


            end if
          end do
        end do

      end do
    end do

    SAFE_DEALLOCATE_A(zeigenval)
    SAFE_DEALLOCATE_A(mmatrix)
    SAFE_DEALLOCATE_A(zeigenvec)
    SAFE_DEALLOCATE_A(dm)
    SAFE_DEALLOCATE_A(zeigref_)
    SAFE_DEALLOCATE_A(zeigplus)
    SAFE_DEALLOCATE_A(zeigminus)
    SAFE_DEALLOCATE_A(zeig0)
    SAFE_DEALLOCATE_A(zeigplusminus)
    SAFE_DEALLOCATE_A(zeigminusplus)
    POP_SUB(lalg_check_zeigenderivatives)
  end subroutine lalg_check_zeigenderivatives

  CMPLX function lalg_zdni(eigenvec, alpha, beta)
    integer, intent(in) :: alpha, beta
    CMPLX :: eigenvec(2)
    lalg_zdni = conjg(eigenvec(alpha)) * eigenvec(beta)
  end function lalg_zdni

  CMPLX function lalg_zduialpha(eigenvec, mmatrix, alpha, gamma, delta)
    integer, intent(in) :: alpha, gamma, delta
    CMPLX :: eigenvec(2), mmatrix(2, 2)
    lalg_zduialpha = mmatrix(alpha, gamma) * eigenvec(delta)
  end function lalg_zduialpha

  CMPLX function lalg_zd2ni(eigenvec, mmatrix, alpha, beta, gamma, delta)
    integer, intent(in) :: alpha, beta, gamma, delta
    CMPLX :: eigenvec(2), mmatrix(2, 2)
    lalg_zd2ni  = conjg(mmatrix(alpha, delta) * eigenvec(gamma)) * eigenvec(beta) + &
                  conjg(eigenvec(alpha)) * mmatrix(beta, gamma) * eigenvec(delta)
  end function lalg_zd2ni

#ifdef HAVE_LAPACK
#include "lalg_adv_lapack_inc.F90"
#endif

#ifdef HAVE_SCALAPACK
#include "undef.F90"
#include "complex.F90"
#include "lalg_adv_scalapack_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "lalg_adv_scalapack_inc.F90"
#endif

end module lalg_adv_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
