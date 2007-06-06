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
!! $Id: specie.F90 2711 2007-02-13 17:36:18Z xavier $

#ifdef DG_VECTORIAL
#define SECIND , 1:3
#define VECDIM 3
#else
#define SECIND
#define VECDIM 1
#endif

subroutine double_grid_apply (this, s, m, sm, x_atom, vl, l, lm, ic)
  type(double_grid_t),    intent(in)    :: this
  type(specie_t),         intent(in)    :: s
  type(mesh_t),           intent(in)    :: m
  type(submesh_t),        intent(in)    :: sm
  FLOAT,                  intent(in)    :: x_atom(1:MAX_DIM)
#ifndef DG_VECTORIAL
  FLOAT,                  intent(out)   :: vl(:)
#else
  FLOAT,                  intent(out)   :: vl(:, :)
#endif
  integer, optional,      intent(in)    :: l
  integer, optional,      intent(in)    :: lm
  integer, optional,      intent(in)    :: ic



  FLOAT :: r, x(1:3)
#ifndef DG_VECTORIAL
  FLOAT :: vv, tmp(1:3)
#else 
  FLOAT :: vv(1:3), tmp
#endif
  integer :: is, is2
  integer :: ii, jj, kk, ll, mm, nn
  integer :: start(1:3), pp, qq, rr

#ifndef DG_VECTORIAL 
  FLOAT, allocatable :: vs(:)
#else 
  FLOAT, allocatable :: vs(:,:)
#endif

  if (.not. this%use_double_grid) then 

    do is = 1, sm%ns
      x(1:3) = m%x(sm%jxyz(is), 1:3) - x_atom(1:3)
      r = sqrt(sum(x(1:3)**2))
      calc_pot(vl(is SECIND))
    end do

  else

    ALLOCATE(vs(0:sm%ns_part SECIND), VECDIM*(sm%ns_part+1))

    vs = M_ZERO

    !for each grid point
    do is = 1, sm%ns

      do ii = -this%nn, this%nn
        do jj = -this%nn, this%nn
          do kk = -this%nn, this%nn

            x(1:3) = m%x(sm%jxyz(is), 1:3) + m%h(1:3)/this%spacing_divisor * (/ii, jj, kk/) - x_atom(1:3)
            r = sqrt(sum(x(1:3)**2))

            calc_pot(vv)

            start(1:3) = m%Lxyz(sm%jxyz(is), 1:3) + this%interpolation_min * (/ii, jj, kk/)
            
            pp = start(1)
            do ll = this%interpolation_min, this%interpolation_max
              
              qq = start(2)
              do mm = this%interpolation_min, this%interpolation_max
                
                rr = start(3)
                do nn = this%interpolation_min, this%interpolation_max

                    is2 = sm%jxyz_inv(m%Lxyz_inv(pp, qq, rr))
                    vs(is2 SECIND) = vs(is2 SECIND)  + this%co(ll)*this%co(mm)*this%co(nn)*vv

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

    vl(1:sm%ns SECIND) = vs(1:sm%ns SECIND)/(this%spacing_divisor**3)

    deallocate(vs)

  end if

end subroutine double_grid_apply

#undef SECIND
#undef VECDIM

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
