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

subroutine double_grid_apply (this, s, m, sm, x_atom, vl, l, lm, ic)
  type(double_grid_t),    intent(in)    :: this
  type(species_t),        intent(in)    :: s
  type(mesh_t),           intent(in)    :: m
  type(submesh_t),        intent(in)    :: sm
  FLOAT,                  intent(in)    :: x_atom(1:MAX_DIM)
  FLOAT,                  intent(out)   :: vl(:)

  integer, optional,      intent(in)    :: l
  integer, optional,      intent(in)    :: lm
  integer, optional,      intent(in)    :: ic

  FLOAT :: r, x(1:MAX_DIM)
  FLOAT :: vv, tmp(1:MAX_DIM)

  integer :: is, is2, ip
  integer :: ii, jj, kk, ll, mm, nn
  integer :: start(1:3), pp, qq, rr
  integer, allocatable :: jxyz_inv(:)

  FLOAT, allocatable :: vs(:)

#ifdef USE_OMP
  integer(kind=omp_lock_kind) :: lock
#endif

  call push_sub('double_grid_apply.double_grid_apply_(non_)local')

  if (.not. this%use_double_grid) then 
    
    !$omp parallel do private(r, x)
    do is = 1, sm%ns
      r = sm%x(is, 0)
      x(1:3) = sm%x(is, 1:3)
      calc_pot(vl(is))
    end do
    !$omp end parallel do 
    
  else

#ifdef HAVE_MPI    
    if (m%parallel_in_domains .and. .not. conf%devel_version) then
      message(1) = "Double grid is not implemented for parallelization in domains (yet)."
      call write_fatal(1)
    end if
#endif

    call profiling_in(profiler, profiler_label)

    ASSERT(.not. simul_box_is_periodic(m%sb))

    ALLOCATE(jxyz_inv(1:sm%np_part), sm%np_part)
    call submesh_get_inv(sm, jxyz_inv)

#ifdef USE_OMP
    vl(1:sm%ns) = M_ZERO
    call omp_init_lock(lock)
    !$omp parallel private(vs)
#endif

    !$omp critical
    ALLOCATE(vs(0:sm%ns_part), (sm%ns_part+1))
    !$omp end critical

    vs = M_ZERO

    !for each grid point
    !$omp do private(ii, jj, kk, ip, ll, mm, nn, pp, qq, rr, is2, start, vv, tmp, r, x)
    do is = 1, sm%ns_part

      do ii = -this%nn, this%nn
        do jj = -this%nn, this%nn
          do kk = -this%nn, this%nn

            x(1:3) = m%x(sm%jxyz(is), 1:3) + m%h(1:3)/this%spacing_divisor * (/ii, jj, kk/) - x_atom(1:3)
            r = sqrt(sum(x(1:3)**2))

            calc_pot(vv)

            ip = sm%jxyz(is)
#ifdef HAVE_MPI                    
            if (m%parallel_in_domains) then
              !map the local point to a global point
              if (ip <= m%np) then
                !inner points
                ip = ip - 1 + m%vp%xlocal(m%vp%partno)
                ip = m%vp%local(ip)
              else if (ip <= m%np + m%vp%np_ghost(m%vp%partno)) then
                !ghost points
                ip = ip - 1 - m%np + m%vp%xghost(m%vp%partno) 
                ip = m%vp%ghost(ip)
              else
                !boundary points
                ip = ip - 1 - (m%np + m%vp%np_ghost(m%vp%partno)) + m%vp%xbndry(m%vp%partno)
                ip = m%vp%bndry(ip)
              end if
            end if
#endif
            start(1:3) = m%idx%Lxyz(ip, 1:3) + this%interpolation_min * (/ii, jj, kk/)
            
            pp = start(1)
            do ll = this%interpolation_min, this%interpolation_max
              
              qq = start(2)
              do mm = this%interpolation_min, this%interpolation_max
                
                rr = start(3)
                do nn = this%interpolation_min, this%interpolation_max

                    ip = m%idx%Lxyz_inv(pp, qq, rr)
#ifdef HAVE_MPI      
                    !map the global point to a local point
                    if (m%parallel_in_domains) ip = vec_global2local(m%vp, ip, m%vp%partno)
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

    deallocate(vs)

#ifdef USE_OMP
    !$omp end parallel
    call omp_destroy_lock(lock)
#endif

    call profiling_out(profiler)

  end if

  call pop_sub()

end subroutine double_grid_apply

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
