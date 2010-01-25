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
  type(double_grid_t),    intent(in)    :: this
  type(species_t),        intent(in)    :: spec
  type(mesh_t),           intent(in)    :: mesh
  type(submesh_t),        intent(in)    :: sm
  FLOAT,                  intent(in)    :: x_atom(1:MAX_DIM)
  FLOAT,                  intent(out)   :: vl(:)

  integer, optional,      intent(in)    :: l
  integer, optional,      intent(in)    :: lm
  integer, optional,      intent(in)    :: ic

  FLOAT :: r, xx(1:MAX_DIM)
  FLOAT :: vv, tmp(1:MAX_DIM)

  integer :: is, is2, ip
  integer :: ii, jj, kk, ll, mm, nn
  integer :: start(1:3), pp, qq, rr
  integer, allocatable :: jxyz_inv(:)
  FLOAT, allocatable :: vs(:)
  type(ps_t), pointer :: ps

#ifdef USE_OMP
  integer(kind=omp_lock_kind) :: lock
#endif

  ! inclusion from double_grid.F90 results in local and non-local versions
  call push_sub('double_grid_apply.double_grid_apply')

  ps => species_ps(spec)
  xx = M_ZERO !must be initialized so xx(dim + 1:MAX_DIM) = 0.0

  if (.not. this%use_double_grid) then 

    do is = 1, sm%ns
      rr = sm%x(is, 0)
      xx(1:3) = sm%x(is, 1:3)
      calc_pot(vl(is))
    end do

  else

    call profiling_in(profiler, profiler_label)

    ASSERT(.not. simul_box_is_periodic(mesh%sb))

    SAFE_ALLOCATE(jxyz_inv(1:sm%np_part))
    call submesh_get_inv(sm, jxyz_inv)

#ifdef USE_OMP
    vl(1:sm%ns) = M_ZERO
    call omp_init_lock(lock)
    !$omp parallel private(vs)
#endif

    !$omp critical
    SAFE_ALLOCATE(vs(0:sm%ns_part))
    !$omp end critical

    vs = M_ZERO

    !for each grid point
    !$omp do private(ii, jj, kk, ip, ll, mm, nn, pp, qq, rr, is2, start, vv, tmp, rr, xx)
    do is = 1, sm%ns_part

      do ii = -this%nn, this%nn
        do jj = -this%nn, this%nn
          do kk = -this%nn, this%nn

            xx(1:3) = mesh%x(sm%jxyz(is), 1:3) + mesh%spacing(1:3)/this%spacing_divisor * (/ii, jj, kk/) - x_atom(1:3)
            rr = sqrt(sum(xx(1:3)**2))

            calc_pot(vv)

            ip = sm%jxyz(is)
#ifdef HAVE_MPI                    
            if (mesh%parallel_in_domains) then
              !map the local point to a global point
              if (ip <= mesh%np) then
                !inner points
                ip = ip - 1 + mesh%vp%xlocal(mesh%vp%partno)
                ip = mesh%vp%local(ip)
              else if (ip <= mesh%np + mesh%vp%np_ghost(mesh%vp%partno)) then
                !ghost points
                ip = ip - 1 - mesh%np + mesh%vp%xghost(mesh%vp%partno) 
                ip = mesh%vp%ghost(ip)
              else
                !boundary points
                ip = ip - 1 - (mesh%np + mesh%vp%np_ghost(mesh%vp%partno)) + mesh%vp%xbndry(mesh%vp%partno)
                ip = mesh%vp%bndry(ip)
              end if
            end if
#endif
            start(1:3) = mesh%idx%Lxyz(ip, 1:3) + this%interpolation_min * (/ii, jj, kk/)
            
            pp = start(1)
            do ll = this%interpolation_min, this%interpolation_max
              
              qq = start(2)
              do mm = this%interpolation_min, this%interpolation_max
                
                rr = start(3)
                do nn = this%interpolation_min, this%interpolation_max

                    ip = mesh%idx%Lxyz_inv(pp, qq, rr)
#ifdef HAVE_MPI      
                    !map the global point to a local point
                    if (mesh%parallel_in_domains) ip = vec_global2local(mesh%vp, ip, mesh%vp%partno)
#endif
                    is2 = jxyz_inv(ip)
                    vs(is2) = vs(is2)  + this%co(ll)*this%co(mm)*this%co(nn)*vv

                    rr = rr + kk
                  end do
                  qq = qq + jj
                end do
                pp = pp + ii
              end do


          end do !kk
        end do !jj
      end do !ii

    end do !is
    !$omp end do nowait

#ifndef USE_OMP
    vl(1:sm%ns) = vs(1:sm%ns)/(this%spacing_divisor**3)
#else
    call omp_set_lock(lock)
    vl(1:sm%ns) = vl(1:sm%ns) + vs(1:sm%ns)/(this%spacing_divisor**3)
    call omp_unset_lock(lock)
#endif

    SAFE_DEALLOCATE_A(vs)

#ifdef USE_OMP
    !$omp end parallel
    call omp_destroy_lock(lock)
#endif

    call profiling_out(profiler)

  end if

  nullify(ps)

  call pop_sub()

end subroutine double_grid_apply

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
