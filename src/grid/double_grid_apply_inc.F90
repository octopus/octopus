!! Copyright (C) 2002-2006 X. Andrade
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
!! $Id: double_grid_apply.F90 2711 2007-02-13 17:36:18Z xavier $

subroutine double_grid_apply (this, spec, mesh, sm, x_atom, vl, l, lm, ic)
  type(double_grid_t),     intent(in)    :: this
  type(species_t), target, intent(in)    :: spec
  type(mesh_t),            intent(in)    :: mesh
  type(submesh_t),         intent(in)    :: sm
  FLOAT,                   intent(in)    :: x_atom(1:MAX_DIM)
  FLOAT,                   intent(out)   :: vl(:)
  integer, optional,       intent(in)    :: l
  integer, optional,       intent(in)    :: lm
  integer, optional,       intent(in)    :: ic

  FLOAT :: r, xx(1:MAX_DIM)
  FLOAT :: vv
  FLOAT :: tmp(1:MAX_DIM)

  integer :: is
  integer :: ii, jj, kk
  integer, allocatable :: map_inv(:)
  FLOAT,   allocatable :: vs(:)
  type(ps_t), pointer :: ps

#ifdef HAVE_OPENMP
  integer(kind=omp_lock_kind) :: lock
#endif

  ! inclusion from double_grid.F90 results in local and non-local versions
  PUSH_SUB(double_grid_apply)

  ps => species_ps(spec)
  xx = M_ZERO !must be initialized so xx(dim + 1:MAX_DIM) = 0.0

  if (.not. this%use_double_grid) then 

    do is = 1, sm%np
      r = sm%x(is, 0)
      xx(1:3) = sm%x(is, 1:3)
      calc_pot(vl(is))
    end do

  else

    call profiling_in(profiler, profiler_label)

    ASSERT(.not. simul_box_is_periodic(mesh%sb))

    SAFE_ALLOCATE(map_inv(1:sm%mesh%np_part))
    call submesh_get_inv(sm, map_inv)

#ifdef HAVE_OPENMP
    vl(1:sm%np) = M_ZERO
    call omp_init_lock(lock)
    !$omp parallel private(vs)
#endif

    !$omp critical
    SAFE_ALLOCATE(vs(0:sm%np_part))
    !$omp end critical

    vs = M_ZERO

    !for each grid point
    !$omp do private(ii, jj, kk, vv, tmp, r, xx)
    do is = 1, sm%np_part

      ! iterate over the fine grid
      do ii = -this%nn, this%nn
        do jj = -this%nn, this%nn
          do kk = -this%nn, this%nn

            xx(1:3) = mesh%x(sm%map(is), 1:3) + mesh%spacing(1:3)/this%spacing_divisor * (/ii, jj, kk/) - x_atom(1:3)
            r = sqrt(sum(xx(1:3)**2))

            ! calculate the value of the potential at that point
            calc_pot(vv)

            ! and apply the corresponding term to all neighbouring grid points
            call apply_to_nb(this, mesh, is, ii, jj, kk, sm%map, map_inv, vv, vs)

          end do !kk
        end do !jj
      end do !ii

    end do !is
    !$omp end do nowait

#ifndef HAVE_OPENMP
    vl(1:sm%np) = vs(1:sm%np)/(this%spacing_divisor**3)
#else
    call omp_set_lock(lock)
    vl(1:sm%np) = vl(1:sm%np) + vs(1:sm%np)/(this%spacing_divisor**3)
    call omp_unset_lock(lock)
#endif

    SAFE_DEALLOCATE_A(vs)

#ifdef HAVE_OPENMP
    !$omp end parallel
    call omp_destroy_lock(lock)
#endif

    call profiling_out(profiler)

  end if

  nullify(ps)

  POP_SUB(double_grid_apply)
end subroutine double_grid_apply

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
