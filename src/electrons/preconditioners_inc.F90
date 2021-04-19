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
subroutine X(preconditioner_apply)(pre, namespace, mesh, hm, a, b, ik, omega)
  type(preconditioner_t),   intent(in)    :: pre
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t), target,     intent(in)    :: mesh
  type(hamiltonian_elec_t), intent(in)    :: hm
  R_TYPE, contiguous,       intent(inout) :: a(:,:)
  R_TYPE, contiguous,       intent(inout) :: b(:,:)
  integer,                  intent(in)    :: ik
  R_TYPE,         optional, intent(in)    :: omega

  integer :: idim
  R_TYPE  :: omega_
  type(profile_t), save :: preconditioner_prof

  type(wfs_elec_t), target :: batch_a
  type(wfs_elec_t) :: batch_b
  type(wfs_elec_t), pointer :: batch_ea

  call profiling_in(preconditioner_prof, TOSTRING(X(PRECONDITIONER)))
  PUSH_SUB(X(preconditioner_apply))

  omega_ = M_ZERO
  if(present(omega)) omega_ = omega

  select case(pre%which)
  case(PRE_NONE)
    call lalg_copy(mesh%np, hm%d%dim, a, b)

  case(PRE_FILTER)
    call wfs_elec_init(batch_a, hm%d%dim, 1, 1, a, ik)
    call wfs_elec_init(batch_b, hm%d%dim, 1, 1, b, ik)
    call boundaries_set(pre%der%boundaries, batch_a)
    if (allocated(hm%hm_base%phase)) then
      SAFE_ALLOCATE(batch_ea)
      call batch_a%copy_to(batch_ea)
      call hamiltonian_elec_base_phase(hm%hm_base, mesh, mesh%np_part, .false., batch_ea, src = batch_a)
      batch_b%has_phase = .true.
    else
      batch_ea => batch_a
    end if

     call X(derivatives_batch_perform)(pre%op, pre%der, batch_ea, batch_b, set_bc = .false.)

    if (allocated(hm%hm_base%phase)) then
      call hamiltonian_elec_base_phase(hm%hm_base, mesh, mesh%np, .true., batch_b)
      call batch_ea%end(copy = .false.)
      SAFE_DEALLOCATE_P(batch_ea)
    end if
    call batch_a%end()
    call batch_b%end()


  case(PRE_JACOBI)
    call apply_D_inverse(a, b)

  case(PRE_POISSON)
    do idim = 1, hm%d%dim
      call X(poisson_solve)(hm%psolver, b(:, idim), a(:, idim), all_nodes=.false.)
    end do
    call lalg_scal(mesh%np, hm%d%dim, R_TOTYPE(M_ONE/(M_TWO*M_PI)), b)

  case(PRE_MULTIGRID)
    call multigrid()

  case default
   write(message(1), '(a,i4,a)') "Unknown preconditioner ", pre%which, "."
   call messages_fatal(1, namespace=namespace)

  end select

  POP_SUB(X(preconditioner_apply))
  call profiling_out(preconditioner_prof)
contains

  subroutine apply_D_inverse(a, b)
    R_TYPE, intent(in)  :: a(:,:)
    R_TYPE, intent(out) :: b(:,:)

    FLOAT, allocatable :: diag(:)

    PUSH_SUB(X(preconditioner_apply).apply_D_inverse)

    SAFE_ALLOCATE(diag(1:mesh%np))

    do idim = 1, hm%d%dim
      diag(:) = pre%diag_lapl(1:mesh%np) + hm%ep%vpsl(1:mesh%np) + hm%vhxc(1:mesh%np, idim)

      b(1:mesh%np,idim) = a(1:mesh%np,idim)/(diag(1:mesh%np) + omega_)
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
    integer :: idim, ip, j

    PUSH_SUB(X(preconditioner_apply).multigrid)

    mesh0 => pre%mgrid%level(0)%mesh
    mesh1 => pre%mgrid%level(1)%mesh
    mesh2 => pre%mgrid%level(2)%mesh

    SAFE_ALLOCATE(d0(1:mesh0%np_part))
    SAFE_ALLOCATE(q0(1:mesh0%np_part))

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
      q1 = M_ZERO
      r2 = M_ZERO
      d2 = M_ZERO
      q2 = M_ZERO


      ! move to level  1
      call X(multigrid_fine2coarse)(pre%mgrid%level(1)%tt, pre%mgrid%level(0)%der, &
        pre%mgrid%level(1)%mesh, a(:,idim), r1, FULLWEIGHT)
      ! r1 has the opposite sign of r2 to avoid an unnecessary operation in the first step

      do ip = 1, mesh1%np
        d1(ip) = -CNST(4.0)*step*r1(ip)
      end do

      ! pre-smoothing
      do j = 1, pre%npre
        call X(derivatives_lapl)(pre%mgrid%level(1)%der, d1, q1, factor = -M_HALF)

        do ip = 1, mesh1%np
          q1(ip) = q1(ip) + r1(ip)
          d1(ip) = d1(ip) - CNST(4.0)*step*q1(ip)
        end do
      end do

      call X(derivatives_lapl)(pre%mgrid%level(1)%der, d1, q1, factor = -M_HALF)

      call lalg_axpy(mesh1%np, M_ONE, q1, r1)


      ! move to level  2
      call X(multigrid_fine2coarse)(pre%mgrid%level(2)%tt, pre%mgrid%level(1)%der, &
        pre%mgrid%level(2)%mesh, q1, r2, FULLWEIGHT)

      do ip = 1, mesh2%np
        d2(ip) = CNST(16.0)*step*r2(ip)
      end do

      ! Jacobi steps on coarsest grid
      do j = 1, pre%nmiddle
        call X(derivatives_lapl)(pre%mgrid%level(2)%der, d2, q2, factor = -M_HALF)

        do ip = 1, mesh2%np
          q2(ip) = q2(ip) - r2(ip)
          d2(ip) = d2(ip) - CNST(16.0)*step*q2(ip)
        end do
      end do

      ! back to level 1
      call X(multigrid_coarse2fine)(pre%mgrid%level(2)%tt, pre%mgrid%level(2)%der, &
        pre%mgrid%level(1)%mesh, d2, t1)

      do ip = 1, mesh1%np
        d1(ip) = d1(ip) - t1(ip)
      end do

      ! post-smoothing
      do j = 1, pre%npost
        call X(derivatives_lapl)(pre%mgrid%level(1)%der, d1, q1, factor = -M_HALF)

        do ip = 1, mesh1%np
          q1(ip) = q1(ip) + r1(ip)
          d1(ip) = d1(ip) - CNST(4.0)*step*q1(ip)
        end do
      end do

      ! and finally back to level 0
      call X(multigrid_coarse2fine)(pre%mgrid%level(1)%tt ,pre%mgrid%level(1)%der, &
        pre%mgrid%level(0)%mesh, d1, q0)

      do ip = 1, mesh0%np
         d0(ip) = - q0(ip)
      end do

      ! post-smoothing
      do j = 1, pre%npost
        call X(derivatives_lapl)(pre%mgrid%level(0)%der, d0, q0, factor = -M_HALF)

        do ip = 1, mesh0%np
          q0(ip) = q0(ip) - a(ip, idim)
          d0(ip) = d0(ip) - step*q0(ip)
        end do
      end do

      do ip = 1, mesh0%np
        b(ip, idim) = -d0(ip)
      end do

    end do

    POP_SUB(X(preconditioner_apply).multigrid)
  end subroutine multigrid


end subroutine X(preconditioner_apply)

! ----------------------------------------

subroutine X(preconditioner_apply_batch)(pre, namespace, mesh, hm, aa, bb, ik, omega)
  type(preconditioner_t),   intent(in)    :: pre
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  type(hamiltonian_elec_t), intent(in)    :: hm
  class(batch_t),           intent(inout) :: aa
  class(batch_t),           intent(inout) :: bb
  integer,                  intent(in)    :: ik
  R_TYPE,         optional, intent(in)    :: omega(:)

  integer :: ii
  type(profile_t), save :: prof
  R_TYPE, allocatable :: psia(:, :), psib(:, :)

  PUSH_SUB(X(preconditioner_apply_batch))
  call profiling_in(prof, TOSTRING(X(PRECONDITIONER_BATCH)))

  call aa%check_compatibility_with(bb)

  if(pre%which == PRE_FILTER) then

    call X(derivatives_batch_perform)(pre%op, pre%der, aa, bb)

  else if(pre%which == PRE_NONE) then

    call aa%copy_data_to(mesh%np, bb)

  else if(pre%which == PRE_MULTIGRID) then

    call multigrid_batch()

  else
    SAFE_ALLOCATE(psia(1:mesh%np_part, 1:hm%d%dim))
    SAFE_ALLOCATE(psib(1:mesh%np, 1:hm%d%dim))
    do ii = 1, aa%nst
      call batch_get_state(aa, ii, mesh%np, psia)
      if(present(omega)) then
        call X(preconditioner_apply)(pre, namespace, mesh, hm, psia, psib, ik, omega(ii))
      else
        call X(preconditioner_apply)(pre, namespace, mesh, hm, psia, psib, ik)
      end if
      call batch_set_state(bb, ii, mesh%np, psib)
    end do
    SAFE_DEALLOCATE_A(psia)
    SAFE_DEALLOCATE_A(psib)
  end if

  call profiling_out(prof)
  POP_SUB(X(preconditioner_apply_batch))

  contains
   ! -----------------------------------------------------

  subroutine multigrid_batch
    FLOAT :: step

    class(batch_t), allocatable :: d0, q0
    class(batch_t), allocatable :: r1, d1, q1, t1
    class(batch_t), allocatable :: r2, d2, q2

    type(mesh_t), pointer :: mesh0, mesh1, mesh2
    integer :: j

    PUSH_SUB(X(preconditioner_apply_batch).multigrid_batch)

    mesh0 => pre%mgrid%level(0)%mesh
    mesh1 => pre%mgrid%level(1)%mesh
    mesh2 => pre%mgrid%level(2)%mesh

    step = CNST(0.66666666)/pre%diag_lapl(1)


    ! move to level  1
    call aa%clone_to(r1, new_np = mesh1%np_part)
    call X(multigrid_fine2coarse_batch)(pre%mgrid%level(1)%tt, pre%mgrid%level(0)%der, &
        pre%mgrid%level(1)%mesh, aa, r1, FULLWEIGHT)
    ! r1 has the opposite sign of r2 to avoid an unnecessary operation in the first step

    call r1%clone_to(d1, copy_data=.true.)
    call batch_scal(mesh1%np, -CNST(4.0)*step, d1)
    call r1%clone_to(q1)

    ! pre-smoothing
    do j = 1, pre%npre
      call X(derivatives_batch_perform)(pre%mgrid%level(1)%der%lapl, pre%mgrid%level(1)%der, d1, q1, factor = -M_HALF)
      
      call batch_axpy(mesh1%np, M_ONE, r1, q1)
      call batch_axpy(mesh1%np, -CNST(4.0)*step, q1, d1)
    end do

    call X(derivatives_batch_perform)(pre%mgrid%level(1)%der%lapl, pre%mgrid%level(1)%der, d1, q1, factor = -M_HALF)

    ! move to level  2
    call q1%clone_to(r2, new_np = mesh2%np_part)
    call X(multigrid_fine2coarse_batch)(pre%mgrid%level(2)%tt, pre%mgrid%level(1)%der, &
        pre%mgrid%level(2)%mesh, q1, r2, FULLWEIGHT)

    call r2%clone_to(d2, copy_data=.true., new_np = mesh2%np_part) 
    call batch_scal(mesh2%np, CNST(16.0)*step, d2)
    
    call r2%clone_to(q2, new_np = mesh2%np)

    ! Jacobi steps on coarsest grid
    do j = 1, pre%nmiddle
      call X(derivatives_batch_perform)(pre%mgrid%level(2)%der%lapl, pre%mgrid%level(2)%der, d2, q2, factor = M_HALF)

      call batch_axpy(mesh2%np, M_ONE, r2, q2)
      call batch_axpy(mesh2%np, CNST(16.0)*step, q2, d2)
    end do

    call r2%end()
    call q2%end()
    SAFE_DEALLOCATE_A(r2)
    SAFE_DEALLOCATE_A(q2)

    ! back to level 1
    call d1%clone_to(t1)
    call X(multigrid_coarse2fine_batch)(pre%mgrid%level(2)%tt, pre%mgrid%level(2)%der, &
        pre%mgrid%level(1)%mesh, d2, t1)

    call d2%end()
    SAFE_DEALLOCATE_A(d2)

    call batch_axpy(mesh1%np, -M_ONE, t1, d1)

    call t1%end()
    SAFE_DEALLOCATE_A(t1)

    ! post-smoothing
    do j = 1, pre%npost
      call X(derivatives_batch_perform)(pre%mgrid%level(1)%der%lapl, pre%mgrid%level(1)%der, d1, q1, factor = -M_HALF)

      call batch_axpy(mesh1%np, M_ONE, r1, q1)
      call batch_axpy(mesh1%np, -CNST(4.0)*step, q1, d1)
    end do

    call q1%end()
    call r1%end()
    SAFE_DEALLOCATE_A(q1)
    SAFE_DEALLOCATE_A(r1)

    ! and finally back to level 0
    call aa%clone_to(q0)
    call X(multigrid_coarse2fine_batch)(pre%mgrid%level(1)%tt ,pre%mgrid%level(1)%der, &
        pre%mgrid%level(0)%mesh, d1, q0)
 
    call d1%end()
    SAFE_DEALLOCATE_A(d1)

    call q0%clone_to(d0, copy_data=.true., new_np = mesh0%np_part)
    call batch_scal(mesh0%np, -M_ONE, d0)

    ! post-smoothing
    do j = 1, pre%npost
      call X(derivatives_batch_perform)(pre%mgrid%level(0)%der%lapl, pre%mgrid%level(0)%der, d0, q0, factor = M_HALF)

      call batch_axpy(mesh0%np, M_ONE, aa, q0)
      call batch_axpy(mesh0%np, step, q0, d0)
    end do

    call q0%end()
    SAFE_DEALLOCATE_A(q0)

    call batch_set_zero(bb)
    call batch_axpy(mesh0%np, -M_ONE, d0, bb)

    call d0%end()
    SAFE_DEALLOCATE_A(d0)

    POP_SUB(X(preconditioner_apply_batch).multigrid_batch)
  end subroutine multigrid_batch

end subroutine X(preconditioner_apply_batch)
!! Local Variables:
!! mode: f90
!! coding: utf-8

