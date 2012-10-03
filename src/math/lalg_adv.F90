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

module lalg_adv_m
  use global_m
  use lapack_m
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
    zlalg_phi

  type(profile_t), save :: cholesky_prof, eigensolver_prof

  interface lalg_cholesky
    module procedure dcholesky, zcholesky
  end interface

  interface lalg_geneigensolve
    module procedure dgeneigensolve, zgeneigensolve
  end interface

  interface lalg_eigensolve_nonh
    module procedure zeigensolve_nonh, deigensolve_nonh
  end interface

  interface lalg_eigensolve
    module procedure deigensolve, zeigensolve
#ifdef HAVE_SCALAPACK
    module procedure deigensolve_scalapack, zeigensolve_scalapack
#endif
  end interface

  ! Note that lalg_determinant and lalg_inverter are just wrappers
  ! over the same routine.
  interface lalg_determinant
    module procedure ddeterminant, zdeterminant
  end interface

  interface lalg_inverter
    module procedure ddeterminant, zdeterminant
  end interface

  interface lalg_sym_inverter
    module procedure dsym_inverter, zsym_inverter
  end interface

  interface lalg_linsyssolve
    module procedure dlinsyssolve, zlinsyssolve
  end interface

  interface lalg_singular_value_decomp
    module procedure zsingular_value_decomp
  end interface

  interface lalg_svd_inverse
    module procedure zsvd_inverse
  end interface

  interface lalg_invert_upper_triangular
    module procedure dinvert_upper_triangular, zinvert_upper_triangular
  end interface
  
  interface lalg_invert_lower_triangular
    module procedure dinvert_lower_triangular, zinvert_lower_triangular
  end interface
  
  interface lalg_lowest_geneigensolve
    module procedure dlowest_geneigensolve, zlowest_geneigensolve
  end interface

  interface lalg_lowest_eigensolve
    module procedure dlowest_eigensolve, zlowest_eigensolve
  end interface

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
