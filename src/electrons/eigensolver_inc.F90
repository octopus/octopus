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
  subroutine X(eigensolver_run)(eigens, namespace, mesh, st, hm, iter, ik)
    type(eigensolver_t),      intent(inout) :: eigens
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    integer,                  intent(in)    :: iter, ik

    integer :: maxiter

    PUSH_SUB(X(eigensolver_run)) 

    maxiter = eigens%es_maxiter

    select case(eigens%es_type)
    case(RS_CG_NEW)
      call X(eigensolver_cg_jiang)(namespace, mesh, st, hm, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
    case(RS_CG)
      call X(eigensolver_cg2)(namespace, mesh, st, hm, hm%xc, eigens%pre, eigens%tolerance, maxiter, &
        eigens%converged(ik), ik, eigens%diff(:, ik), eigens%orthogonalize_to_all, &
        eigens%conjugate_direction, eigens%additional_terms, eigens%energy_change_threshold)
    case(RS_PLAN)
      call X(eigensolver_plan)(namespace, mesh, st, hm, eigens%pre, eigens%tolerance, maxiter, eigens%converged(ik), ik, &
        eigens%diff(:, ik))
    case(RS_EVO)
      maxiter = 1
      call X(eigensolver_evolution)(namespace, mesh, st, hm, eigens%exponential_operator, eigens%tolerance, maxiter, &
        eigens%converged(ik), ik, eigens%diff(:, ik), tau = eigens%imag_time)
    case(RS_RMMDIIS)
      if(iter <= eigens%rmmdiis_minimization_iter) then
        maxiter = 2
        call X(eigensolver_rmmdiis_min)(namespace, mesh, st, hm, eigens%pre, maxiter, eigens%converged(ik), ik)
      else
        call X(eigensolver_rmmdiis)(namespace, mesh, st, hm, eigens%pre, eigens%tolerance, maxiter, &
          eigens%converged(ik), ik, eigens%diff(:, ik))
      end if
    end select

    ! Do subspace diagonalization always after the eigensolver
    ! Tests have shown that this leads to better convergence
    if(st%calc_eigenval) then
      call X(subspace_diag)(eigens%sdiag, namespace, mesh, st, hm, ik, st%eigenval(:, ik), eigens%diff(:, ik))
    end if


    eigens%matvec = eigens%matvec + maxiter

    POP_SUB(X(eigensolver_run))
  end subroutine X(eigensolver_run)
