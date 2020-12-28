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

  ! ---------------------------------------------------------
  subroutine X(eigensolver_run)(eigens, namespace, gr, st, hm, iter, ik)
    type(eigensolver_t),      intent(inout) :: eigens
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    integer,                  intent(in)    :: iter, ik

    integer :: maxiter

    PUSH_SUB(X(eigensolver_run)) 

    maxiter = eigens%es_maxiter

    select case(eigens%es_type)
    case(RS_CG_NEW)
      call X(eigensolver_cg2_new)(namespace, gr, st, hm, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
    case(RS_CG)
      call X(eigensolver_cg)(namespace, gr, st, hm, hm%xc, eigens%pre, eigens%tolerance, maxiter, &
        eigens%converged(ik), ik, eigens%diff(:, ik), eigens%orthogonalize_to_all, &
        eigens%conjugate_direction, eigens%additional_terms)
    case(RS_PLAN)
      call X(eigensolver_plan)(namespace, gr, st, hm, eigens%pre, eigens%tolerance, maxiter, eigens%converged(ik), ik, &
        eigens%diff(:, ik))
    case(RS_EVO)
      maxiter = 1
      call X(eigensolver_evolution)(namespace, gr%mesh, st, hm, eigens%exponential_operator, eigens%tolerance, maxiter, &
        eigens%converged(ik), ik, eigens%diff(:, ik), tau = eigens%imag_time)
    case(RS_LOBPCG)
      call X(eigensolver_lobpcg)(namespace, gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
        eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size)
    case(RS_RMMDIIS)
      if(iter <= eigens%rmmdiis_minimization_iter) then
        maxiter = 2
        call X(eigensolver_rmmdiis_min)(namespace, gr, st, hm, eigens%pre, maxiter, eigens%converged(ik), ik)
      else
        call X(eigensolver_rmmdiis)(namespace, gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
          eigens%converged(ik), ik, eigens%diff(:, ik))
      end if
    case(RS_PSD)
      call X(eigensolver_rmmdiis_min)(namespace, gr, st, hm, eigens%pre, maxiter, eigens%converged(ik), ik)
    end select

    ! Do subspace diagonalization always after the eigensolver
    ! Tests have shown that this leads to better convergence
    if(st%calc_eigenval) then
      call X(subspace_diag)(eigens%sdiag, namespace, gr%mesh, st, hm, ik, st%eigenval(:, ik), eigens%diff(:, ik))
    end if


    eigens%matvec = eigens%matvec + maxiter

    POP_SUB(X(eigensolver_run))
  end subroutine X(eigensolver_run)
