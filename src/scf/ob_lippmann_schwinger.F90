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
!! $Id: scf.F90 4182 2008-05-14 14:02:30Z acastro $

! This module solves the Schroedinger equation for a system with open
! boundaries for a prescribed energy.

#include "global.h"

module ob_lippmann_schwinger_m
  use blas_m
  use eigensolver_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use lalg_basic_m
  use math_m
  use messages_m
  use ob_interface_m
  use ob_lead_m
  use profiling_m
  use simul_box_m
  use solvers_m
  use states_m

  implicit none

  private
  public :: &
    lippmann_schwinger

  ! Pointers to communicate with iterative linear solver.
  integer, pointer             :: ist_p, ik_p
  FLOAT, pointer               :: energy_p
  CMPLX, pointer               :: green_p(:, :, :, :)
  type(grid_t), pointer        :: gr_p
  type(hamiltonian_t), pointer :: hm_p
  type(states_t), pointer      :: st_p

contains

  ! ---------------------------------------------------------
  ! Solve the Lippmann-Schwinger equation for the open boundary
  ! system. Use convergence criteria in eigens.
  subroutine lippmann_schwinger(eigens, hm, gr, st)
    type(eigensolver_t),        intent(out)   :: eigens
    type(hamiltonian_t), target, intent(inout) :: hm
    type(grid_t), target,        intent(inout) :: gr
    type(states_t), target,      intent(inout) :: st

    integer                    :: idim, ip, ip_lead, il, np_intf, iter
    integer, target            :: ist, ik
    FLOAT, target              :: energy
    FLOAT                      :: tol, res
    CMPLX, allocatable         :: rhs(:, :)
    CMPLX, allocatable, target :: green(:, :, :, :)
    logical                    :: conv
    
    call push_sub('ob_lippmann_schwinger.lippmann_schwinger')

    ALLOCATE(rhs(gr%mesh%np_part, st%d%dim), gr%mesh%np_part*st%d%dim)
    ALLOCATE(green(gr%intf(LEFT)%np, gr%intf(LEFT)%np, st%d%dim, NLEADS), gr%intf(LEFT)%np**2*st%d%dim**NLEADS)

    eigens%converged = 0
    eigens%matvec    = 0

    ist_p    => ist
    ik_p     => ik
    green_p  => green
    gr_p     => gr
    hm_p      => hm
    st_p     => st
    energy_p => energy

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        ! Solve Schroedinger equation for this energy.
        energy = st%ob_eigenval(ist, ik)

        ASSERT(ubound(st%zphi, dim = 1) >= gr%mesh%np_part)

        ! Calculate right hand side e-T-V0-sum(a)[H_ca*g_a*H_ac].
        rhs(:, :) = M_z0
        call zhamiltonian_apply(hm, gr, st%zphi(:, :, ist, ik), rhs(:, :), ist, ik, kinetic_only=.true.)

        ! Apply lead potential.
        do idim = 1, st%d%dim
          do ip = 1, gr%mesh%np
            ip_lead = mod(ip-1, gr%intf(LEFT)%np) + 1
            rhs(ip, idim) = rhs(ip, idim) + hm%lead_vks(ip_lead, idim, LEFT)*st%zphi(ip, idim, ist, ik)
          end do
        end do

        ! Add energy.
        do idim = 1, st%d%dim
          call lalg_scal(gr%mesh%np, -M_z1, rhs(:, idim))
          call lalg_axpy(gr%mesh%np, energy, st%zphi(:, idim, ist, ik), rhs(:, idim))
        end do

        ! Apply term with lead Green functions.
        do il = 1, NLEADS
          do idim = 1, st%d%dim
            np_intf = gr%intf(il)%np
            call lalg_copy(np_intf**2, hm%lead_green(:, 1, idim, ist, ik, il), green(:, 1, idim, il))
            call apply_coupling(green(:, :, idim, il), hm%lead_h_offdiag(:, :, il), &
              green(:, :, idim, il), np_intf, il)
          end do
        end do
        do il = 1, NLEADS
          do idim = 1, st%d%dim
            call interface_apply_sym_op(gr%intf(il), -M_z1, green(:, :, idim, il), &
              st%zphi(:, idim, ist, ik), M_z1, rhs(:, idim))
          end do
        end do

        ! Solve linear system lhs psi = rhs.
        iter = eigens%es_maxiter
        tol  = eigens%final_tol

        conv = .false.
        call zqmr_sym(gr%mesh%np_part*st%d%dim, st%zpsi(:, 1, ist, ik), rhs(:, 1), lhs, dotu, nrm2, precond, &
          iter, residue=res, threshold=tol, converged=conv)

        eigens%matvec = eigens%matvec + iter + 1 + 2
        if(conv) then
          eigens%converged = eigens%converged + 1
        end if
        eigens%diff(ist, ik) = res
      end do
    end do

    call pop_sub()
  end subroutine lippmann_schwinger


  ! ---------------------------------------------------------
  ! Dot product for QMR solver, works for x being a spinor.
  CMPLX function dotu(x, y)
    CMPLX, intent(in) :: x(:)
    CMPLX, intent(in) :: y(:)

    integer :: np_part, np, dim, idim
    CMPLX   :: dot

    call push_sub('ob_lippmann_schwinger.dotu')

    np_part = gr_p%mesh%np_part
    np      = gr_p%mesh%np
    dim     = st_p%d%dim

    dot = M_ZERO
    do idim = 1, dim
      dot = dot + blas_dotu(np, x(l(idim)), 1, y(l(idim)), 1)
    end do
    dotu = dot

    call pop_sub()

  contains

    integer function l(idim)
      integer, intent(in) :: idim

      l = (idim-1)*np_part+1
    end function l
  end function dotu


  ! ---------------------------------------------------------
  ! Norm for QMR solver. This routine works for x being a spinor
  ! considered as one vector.
  FLOAT function nrm2(x)
    CMPLX, intent(in) :: x(:)

    integer :: np_part, np, dim, idim
    FLOAT   :: nrm

    call push_sub('ob_lippmann_schwinger.nrm2')

    np_part = gr_p%mesh%np_part
    np      = gr_p%mesh%np
    dim     = st_p%d%dim

    nrm = M_ZERO
    do idim = 1, dim
      nrm = nrm + lalg_nrm2(np, x(l(idim):u(idim)))
    end do
    nrm2 = nrm

    call pop_sub()

  contains

    integer function l(idim)
      integer, intent(in) :: idim

      l = (idim-1)*np_part+1
    end function l

    integer function u(idim)
      integer, intent(in) :: idim

      u = l(idim)+np_part-1
    end function u
  end function nrm2

  
  ! ---------------------------------------------------------
  ! The left hand side of the Lippmann-Schwinger equation
  ! e-H-sum(a)[H_ca*g_a*H_ac].
  ! Used by the iterative linear solver.
  subroutine lhs(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp_x(:, :)
    CMPLX, allocatable :: tmp_y(:, :)
    integer            :: np, np_part, idim, il, dim

    call push_sub('ob_lippmann_schwinger.lhs')

    np      = gr_p%mesh%np
    np_part = gr_p%mesh%np_part
    dim     = st_p%d%dim

    ALLOCATE(tmp_x(np_part, dim), np_part*dim)
    ALLOCATE(tmp_y(np_part, dim), np_part*dim)

    do idim = 1, dim
      call lalg_copy(np, x(l(idim):u(idim)), tmp_x(:, idim))
    end do
    call zhamiltonian_apply(hm_p, gr_p, tmp_x, tmp_y, ist_p, ik_p)

    ! y <- e x - tmp_y
    do idim = 1, dim
      call lalg_copy(np, tmp_y(:, idim), y(l(idim):u(idim)))
      call lalg_scal(np, -M_z1, y(l(idim):u(idim)))
      call lalg_axpy(np, energy_p, x(l(idim):u(idim)), y(l(idim):u(idim)))
    end do

    do il = 1, NLEADS
      do idim = 1, dim
        call interface_apply_sym_op(gr_p%intf(il), &
          -M_z1, green_p(:, :, idim, il), tmp_x(:, idim), M_z1, y(l(idim):u(idim)))
      end do
    end do

    deallocate(tmp_x, tmp_y)
    call pop_sub()

  contains

    integer function l(idim)
      integer, intent(in) :: idim

      l = (idim-1)*np_part+1
    end function l

    integer function u(idim)
      integer, intent(in) :: idim

      u = l(idim)+np_part-1
    end function u
  end subroutine lhs


  ! ---------------------------------------------------------
  ! Identity preconditioner. Since preconditioning with the inverse of
  ! the diagonal did not improve the convergence we put identity here
  ! until we have something better.
  subroutine precond(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    call push_sub('ob_lippmann_schwinger.precond')

    y = x

    call pop_sub()
  end subroutine precond
end module ob_lippmann_schwinger_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
