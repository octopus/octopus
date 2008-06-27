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

! --------------------------------------------------------- 
subroutine X(preconditioner_apply)(pre, gr, h, a, b, omega)
  type(preconditioner_t), intent(in)    :: pre
  type(grid_t),           intent(inout) :: gr
  type(hamiltonian_t),    intent(inout) :: h
  R_TYPE,                 intent(inout) :: a(:,:)
  R_TYPE,                 intent(out)   :: b(:,:)
  R_TYPE,       optional, intent(in)    :: omega
  
  integer :: idim
  FLOAT   :: omega_
  type(profile_t), save :: preconditioner_prof

  call profiling_in(preconditioner_prof, "PRECONDITIONER")

  omega_ = M_ZERO
  if(present(omega)) omega_ = omega

  select case(pre%which)
  case(PRE_NONE)
    do idim = 1, h%d%dim
      call lalg_copy(NP, a(:,idim), b(:,idim))
    end do

  case(PRE_SMOOTHING)
    do idim = 1, h%d%dim
      call X(derivatives_oper)(pre%op, gr%f_der%der_discr, a(:, idim), b(:, idim))
    end do

  case(PRE_JACOBI)
    call apply_D_inverse(a, b)

  case(PRE_POISSON)
    do idim = 1, h%d%dim
      call X(poisson_solve) (gr, b(:, idim), a(:, idim), all_nodes=.false.)
      call lalg_scal(NP, R_TOTYPE(M_ONE/(M_TWO*M_PI)), b(:,idim))
    end do

  case default
   write(message(1), '(a,i4,a)') "Error: unknown preconditioner ", pre%which, "."
   call write_fatal(1)

  end select

  call profiling_out(preconditioner_prof)

contains

  subroutine apply_D_inverse(a, b)
    R_TYPE, intent(in)  :: a(:,:)
    R_TYPE, intent(out) :: b(:,:)

    FLOAT, allocatable :: diag(:)

    ALLOCATE(diag(NP), NP)

    do idim = 1, h%d%dim
      diag(:) = pre%diag_lapl(1:NP) + h%ep%vpsl(1:NP) + h%vhxc(1:NP, idim)

      b(1:NP,idim) = a(1:NP,idim)/(diag(1:NP) + omega_)
    end do

    deallocate(diag)
  end subroutine apply_D_inverse

end subroutine X(preconditioner_apply)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
