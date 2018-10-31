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
subroutine X(preconditioner_apply)(pre, gr, hm, ik, a, b, omega, energy, ist)
  type(preconditioner_t), intent(in)    :: pre
  type(grid_t), target,   intent(in)    :: gr
  type(hamiltonian_t),    intent(in)    :: hm
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: a(:,:)
  R_TYPE,                 intent(inout) :: b(:,:)
  R_TYPE,       optional, intent(in)    :: omega
  FLOAT,        optional, intent(in)    :: energy
  integer,      optional, intent(in)    :: ist

  integer :: idim, iter
  R_TYPE  :: omega_
  type(profile_t), save :: preconditioner_prof
  FLOAT :: kinetic_energy, energy_
  R_TYPE, allocatable :: xx(:,:), xxsq(:,:), tmp(:,:)
  R_TYPE, allocatable :: a1(:,:), a2(:,:)
  FLOAT :: jacobi_relax

  call profiling_in(preconditioner_prof, "PRECONDITIONER")
  PUSH_SUB(X(preconditioner_apply))

  omega_ = M_ZERO
  if(present(omega)) omega_ = omega

  select case(pre%which)
  case(PRE_NONE)
    do idim = 1, hm%d%dim
      call lalg_copy(gr%mesh%np, a(:,idim), b(:,idim))
    end do

  case(PRE_FILTER)
    do idim = 1, hm%d%dim
      call X(derivatives_perform)(pre%op, gr%der, a(:, idim), b(:, idim))
    end do

  case(PRE_JACOBI)
    call apply_D_inverse(a, b)

  case(PRE_POISSON)
    do idim = 1, hm%d%dim
      call X(poisson_solve)(psolver, b(:, idim), a(:, idim), all_nodes=.false.)
      call lalg_scal(gr%mesh%np, R_TOTYPE(M_ONE/(M_TWO*M_PI)), b(:,idim))
    end do

  case(PRE_MULTIGRID)
    call multigrid()

  case(PRE_ENERGY)
    ! these arguments are necessary for the preconditioner to work
    ASSERT(present(energy))
    ASSERT(present(ist))
    SAFE_ALLOCATE(xx(1:gr%mesh%np, 1:hm%d%dim))
    SAFE_ALLOCATE(xxsq(1:gr%mesh%np, 1:hm%d%dim))
    SAFE_ALLOCATE(tmp(1:gr%mesh%np, 1:hm%d%dim))

    ! Calculate kinetic energy
    call X(hamiltonian_apply)(hm, gr%der, a, tmp, ist, ik, terms=TERM_KINETIC)
    kinetic_energy = R_REAL(X(mf_dotp) (gr%mesh, hm%d%dim, a, tmp))

    do idim = 1, hm%d%dim
      ! is the indexing of the potential correct for spinor cases?
      xx(:,idim) = pre%energy_ratio_factor * abs(energy - hm%vhxc(:,idim)) /&
        kinetic_energy
      xxsq(:,idim) = xx(:,idim) * xx(:,idim)
      tmp(:,idim) = CNST(27.0) + CNST(18.0)*xx(:,idim) + CNST(12.0)*xxsq(:,idim) &
        + CNST(8.0)*xx(:,idim)*xxsq(:,idim)
      b(1:gr%mesh%np, idim) = tmp(:,idim) / (tmp(:,idim) + CNST(16.0)*xxsq(:,idim)*xxsq(:,idim)) &
        * a(1:gr%mesh%np, idim)
    end do

    SAFE_DEALLOCATE_A(xx)
    SAFE_DEALLOCATE_A(xxsq)
    SAFE_DEALLOCATE_A(tmp)

  case(PRE_JACOBI2)
    SAFE_ALLOCATE(a1(1:gr%mesh%np_part, hm%d%dim))
    SAFE_ALLOCATE(a2(1:gr%mesh%np_part, hm%d%dim))

    omega_ = M_ZERO
    jacobi_relax = pre%relax_factor

    do idim = 1, hm%d%dim
      b(1:gr%mesh%np,idim) = jacobi_relax*a(1:gr%mesh%np,idim)/pre%diag_lapl(1:gr%mesh%np)
    end do
    do iter = 1, pre%jacobi_iterations
      do idim = 1, hm%d%dim
        call X(derivatives_lapl)(gr%der, b(:, idim), a2(:, idim))
        b(1:gr%mesh%np,idim) = b(1:gr%mesh%np,idim) + jacobi_relax/pre%diag_lapl(1:gr%mesh%np) &
          *(CNST(-0.5)*a2(1:gr%mesh%np, idim) + a(1:gr%mesh%np, idim))
      end do
    end do

    SAFE_DEALLOCATE_A(a1)
    SAFE_DEALLOCATE_A(a2)

  case default
   write(message(1), '(a,i4,a)') "Unknown preconditioner ", pre%which, "."
   call messages_fatal(1)

  end select

  POP_SUB(X(preconditioner_apply))
  call profiling_out(preconditioner_prof)
contains

  subroutine apply_D_inverse(a, b)
    R_TYPE, intent(in)  :: a(:,:)
    R_TYPE, intent(out) :: b(:,:)

    FLOAT, allocatable :: diag(:)

    PUSH_SUB(X(preconditioner_apply).apply_D_inverse)

    SAFE_ALLOCATE(diag(1:gr%mesh%np))

    do idim = 1, hm%d%dim
      diag(:) = pre%diag_lapl(1:gr%mesh%np) + hm%ep%vpsl(1:gr%mesh%np) + hm%vhxc(1:gr%mesh%np, idim)

      b(1:gr%mesh%np,idim) = a(1:gr%mesh%np,idim)/(diag(1:gr%mesh%np) + omega_)
    end do

    SAFE_DEALLOCATE_A(diag)
    POP_SUB(X(preconditioner_apply).apply_D_inverse)
  end subroutine apply_D_inverse

  ! -----------------------------------------------------

  subroutine multigrid
    FLOAT :: step

    R_TYPE, allocatable :: d0(:), q0(:)
    R_TYPE, allocatable :: r1(:), d1(:), q1(:), t1(:)
    R_TYPE, allocatable :: r2(:), d2(:), q2(:)

    type(mesh_t), pointer :: mesh0, mesh1, mesh2
    integer :: idim, ip

    PUSH_SUB(X(preconditioner_apply).multigrid)

    mesh0 => gr%mgrid_prec%level(0)%mesh
    mesh1 => gr%mgrid_prec%level(1)%mesh
    mesh2 => gr%mgrid_prec%level(2)%mesh

    SAFE_ALLOCATE(d0(1:mesh0%np_part))
    SAFE_ALLOCATE(q0(1:mesh0%np))

    SAFE_ALLOCATE(r1(1:mesh1%np))
    SAFE_ALLOCATE(d1(1:mesh1%np_part))
    SAFE_ALLOCATE(q1(1:mesh1%np_part))
    SAFE_ALLOCATE(t1(1:mesh1%np_part))

    SAFE_ALLOCATE(r2(1:mesh2%np))
    SAFE_ALLOCATE(d2(1:mesh2%np_part))
    SAFE_ALLOCATE(q2(1:mesh2%np))

    step = CNST(0.66666666)/pre%diag_lapl(1)

    do idim = 1, hm%d%dim

      d0 = M_ZERO
      q0 = M_ZERO
      r1 = M_ZERO
      d1 = M_ZERO
      q1 = M_ZERO
      r2 = M_ZERO
      d2 = M_ZERO
      q2 = M_ZERO

      ! move to level  1
      call X(multigrid_fine2coarse)(gr%mgrid_prec%level(1)%tt, gr%mgrid_prec%level(0)%der, &
        gr%mgrid_prec%level(1)%mesh, a(:, idim), r1, FULLWEIGHT)

      forall (ip = 1:mesh1%np)
        r1(ip) = -r1(ip)
        d1(ip) = CNST(4.0)*step*r1(ip)
      end forall

      call X(derivatives_lapl)(gr%mgrid_prec%level(1)%der, d1, q1)

      forall (ip = 1:mesh1%np) q1(ip) = CNST(-0.5)*q1(ip) - r1(ip)

      ! move to level  2

      call X(multigrid_fine2coarse)(gr%mgrid_prec%level(2)%tt, gr%mgrid_prec%level(1)%der, &
        gr%mgrid_prec%level(2)%mesh, q1, r2, FULLWEIGHT)

      forall (ip = 1:mesh2%np) d2(ip) = CNST(16.0)*step*r2(ip)

      call X(derivatives_lapl)(gr%mgrid_prec%level(2)%der, d2, q2)

      forall (ip = 1:mesh2%np)
        q2(ip) = CNST(-0.5)*q2(ip) - r2(ip)
        d2(ip) = d2(ip) - CNST(16.0)*step*q2(ip)
      end forall

      ! back to level 1

      call X(multigrid_coarse2fine)(gr%mgrid_prec%level(2)%tt, gr%mgrid_prec%level(2)%der, &
        gr%mgrid_prec%level(1)%mesh, d2, t1)

      forall (ip = 1:mesh1%np)
        q1(ip) = q1(ip) + t1(ip)
        d1(ip) = d1(ip) - q1(ip)
      end forall

      call X(derivatives_lapl)(gr%mgrid_prec%level(1)%der, d1, q1)

      forall (ip = 1:mesh1%np)
        q1(ip) = CNST(-0.5)*q1(ip) - r1(ip)
        d1(ip) = d1(ip) - CNST(4.0)*step*q1(ip)
      end forall

      ! and finally back to level 0

      call X(multigrid_coarse2fine)(gr%mgrid_prec%level(1)%tt ,gr%mgrid_prec%level(1)%der, &
        gr%mgrid_prec%level(0)%mesh, d1, d0)

      forall (ip = 1:mesh0%np) d0(ip) = -d0(ip)

      call X(derivatives_lapl)(gr%mgrid_prec%level(0)%der, d0, q0)

      forall (ip = 1:mesh0%np)
        q0(ip) = CNST(-0.5)*q0(ip) - a(ip, idim)
        d0(ip) = d0(ip) - step*q0(ip)
        b(ip, idim) = -d0(ip)
      end forall

    end do

    POP_SUB(X(preconditioner_apply).multigrid)
  end subroutine multigrid


end subroutine X(preconditioner_apply)

! ----------------------------------------

subroutine X(preconditioner_apply_batch)(pre, gr, hm, ik, aa, bb, omega)
  type(preconditioner_t), intent(in)    :: pre
  type(grid_t),           intent(in)    :: gr
  type(hamiltonian_t),    intent(in)    :: hm
  integer,                intent(in)    :: ik
  type(batch_t),          intent(inout) :: aa
  type(batch_t),          intent(inout) :: bb
  R_TYPE,       optional, intent(in)    :: omega(:)

  integer :: ii
  type(profile_t), save :: prof

  PUSH_SUB(X(preconditioner_apply_batch))
  call profiling_in(prof, 'PRECONDITIONER_BATCH')

  if(pre%which == PRE_FILTER) then

    call X(derivatives_batch_perform)(pre%op, gr%der, aa, bb)

  else if(pre%which == PRE_NONE) then

    call batch_copy_data(gr%der%mesh%np, aa, bb)

  else
    ASSERT(.not. batch_is_packed(aa))
    ASSERT(.not. batch_is_packed(bb))
    do ii = 1, aa%nst
      if (present(omega)) then
        call X(preconditioner_apply)(pre, gr, hm, ik, aa%states(ii)%X(psi), bb%states(ii)%X(psi), omega(ii))
      else
        call X(preconditioner_apply)(pre, gr, hm, ik, aa%states(ii)%X(psi), bb%states(ii)%X(psi))
      end if
    end do
  end if

  call profiling_out(prof)
  POP_SUB(X(preconditioner_apply_batch))
end subroutine X(preconditioner_apply_batch)
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
