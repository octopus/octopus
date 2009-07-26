!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

  ! ---------------------------------------------------------
  subroutine X(multigrid_coarse2fine)(tt, coarse_der, coarse_mesh, f_coarse, f_fine)
    type(transfer_table_t),  intent(in)    :: tt
    type(derivatives_t),     intent(in)    :: coarse_der
    type(mesh_t),            intent(in)    :: coarse_mesh
    R_TYPE,                  intent(inout) :: f_coarse(:)
    R_TYPE,                  intent(out)   :: f_fine(:)

    FLOAT, pointer :: vol_pp(:)

    integer :: i, i1, i2, i4, i8
    integer :: jj, j(8)
    FLOAT   :: vol_total

    call push_sub('multigrid.Xmultigrid_coarse2fine')

    call profiling_in(interp_prof, "MG_INTERPOLATION")

    call X(set_bc)(coarse_der, f_coarse)

#ifdef HAVE_MPI
    if(coarse_mesh%parallel_in_domains) call X(vec_ghost_update)(coarse_mesh%vp, f_coarse)
#endif

    vol_pp => coarse_mesh%vol_pp

    i1 = 0;  i2 = 0;  i4 = 0;  i8 = 0;
    do i = 1, tt%n_fine
      select case(tt%fine_i(i))
      case(1)
        i1 = i1 + 1
        j(1:1) = tt%to_fine1(1:1, i1)
      case(2)
        i2 = i2 + 1
        j(1:2) = tt%to_fine2(1:2, i2)
      case(4)
        i4 = i4 + 1
        j(1:4) = tt%to_fine4(1:4, i4)
      case(8)
        i8 = i8 + 1
        j(1:8) = tt%to_fine8(1:8, i8)
      end select

      if(coarse_mesh%use_curvilinear) then
        f_fine(i) = M_ZERO
        vol_total = M_ZERO
        do jj = 1, tt%fine_i(i)
          f_fine(i) = f_fine(i) + vol_pp(j(jj))*f_coarse(j(jj))
          vol_total = vol_total + vol_pp(j(jj))
        end do
        f_fine(i) = f_fine(i)/vol_total
      else
        f_fine(i) = M_ZERO
        vol_total = M_ZERO
        do jj = 1, tt%fine_i(i)
          f_fine(i) = f_fine(i) + vol_pp(1)*f_coarse(j(jj))
          vol_total = vol_total + vol_pp(1)
        end do
        f_fine(i) = f_fine(i)/vol_total
      end if

    end do
    call profiling_out(interp_prof)
    call pop_sub()
  end subroutine X(multigrid_coarse2fine)

  ! ---------------------------------------------------------
  subroutine X(multigrid_fine2coarse)(mgrid, ilevel, f_fine, f_coarse, method_p)
    type(multigrid_t),    intent(in)    :: mgrid
    integer,              intent(in)    :: ilevel
    R_TYPE,               intent(inout) :: f_fine(:)
    R_TYPE,               intent(out)   :: f_coarse(:)
    integer, optional,    intent(in)    :: method_p

    integer :: method

    call push_sub('multigrid.multigrid_fine2coarse')

    if(present(method_p)) then
      method=method_p
    else
      method=FULLWEIGHT
    end if

    select case(method)
    case(FULLWEIGHT)
      call X(multigrid_restriction)(mgrid%level(ilevel)%tt, mgrid%level(ilevel - 1)%der, &
        mgrid%level(ilevel - 1)%mesh, mgrid%level(ilevel)%mesh, f_fine, f_coarse)
    case(INJECTION)
      call X(multigrid_injection)(mgrid%level(ilevel)%tt, f_fine, f_coarse)
    case default
      write(message(1), '(a,i2,a)') 'Multigrid: Restriction method  = ', method, ' is not valid.'
      call write_fatal(1)
    end select

    call pop_sub()
  end subroutine X(multigrid_fine2coarse)


  ! ---------------------------------------------------------
  subroutine X(multigrid_injection)(tt, f_fine, f_coarse)
    type(transfer_table_t), intent(in)  :: tt
    R_TYPE,                 intent(in)  :: f_fine(:)
    R_TYPE,                 intent(out) :: f_coarse(:)

    integer :: i

    call push_sub('multigrid.multigrid_injection')
    call profiling_in(injection_prof, "MG_INJECTION")

    do i = 1, tt%n_coarse
      f_coarse(i) = f_fine(tt%to_coarse(i))
    end do

    call profiling_out(injection_prof)
    call pop_sub()
  end subroutine X(multigrid_injection)

  ! ---------------------------------------------------------
  subroutine X(multigrid_restriction)(tt, fine_der, fine_mesh, coarse_mesh, f_fine, f_coarse)
    type(transfer_table_t), intent(in)    :: tt
    type(derivatives_t),    intent(in)    :: fine_der
    type(mesh_t),           intent(in)    :: fine_mesh
    type(mesh_t),           intent(in)    :: coarse_mesh
    R_TYPE,                 intent(inout) :: f_fine(:)
    R_TYPE,                 intent(out)   :: f_coarse(:)

    FLOAT :: weight(-1:1,-1:1,-1:1)

    integer :: n, fn, di, dj, dk, d, fi(MAX_DIM)

    call push_sub('multigrid.multigrid_restriction')
    call profiling_in(restrict_prof, "MG_RESTRICTION")

    do di = -1, 1
      do dj = -1, 1
        do dk = -1, 1
          d = abs(di) + abs(dj) + abs(dk)
          weight(di, dj, dk) = CNST(0.5)**d
        end do
      end do
    end do

    call X(set_bc)(fine_der, f_fine)

#ifdef HAVE_MPI
    if(fine_mesh%parallel_in_domains) call X(vec_ghost_update)(fine_mesh%vp, f_fine)
#endif

    do n = 1, tt%n_coarse
      fn = tt%to_coarse(n)
#ifdef HAVE_MPI
      ! translate to a global index
      if(fine_mesh%parallel_in_domains) fn = fine_mesh%vp%local(fn - 1 + fine_mesh%vp%xlocal(fine_mesh%vp%partno))
#endif
      fi(:) = fine_mesh%idx%Lxyz(fn, :)

      f_coarse(n) = M_ZERO

      do di = -1, 1
        do dj = -1, 1
          do dk = -1, 1
            fn = fine_mesh%idx%Lxyz_inv(fi(1) + di, fi(2) + dj, fi(3) + dk)

#ifdef HAVE_MPI
            ! translate to a local index
            if(fine_mesh%parallel_in_domains) fn = vec_global2local(fine_mesh%vp, fn, fine_mesh%vp%partno)
#endif
            if(fine_mesh%use_curvilinear) then
              f_coarse(n) = f_coarse(n) + weight(di, dj, dk)*f_fine(fn)*fine_mesh%vol_pp(fn)
            else
              f_coarse(n) = f_coarse(n) + weight(di, dj, dk)*f_fine(fn)*fine_mesh%vol_pp(1)
            end if

          end do
        end do
      end do

      if(fine_mesh%use_curvilinear) then
        f_coarse(n) = f_coarse(n)/coarse_mesh%vol_pp(n)
      else
        f_coarse(n) = f_coarse(n)/coarse_mesh%vol_pp(1)
      end if
    end do

    call profiling_count_operations(tt%n_coarse*(27*3 + 1))
    call profiling_out(restrict_prof)
    call pop_sub()
  end subroutine X(multigrid_restriction)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
