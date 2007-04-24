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

#include "global.h"

module double_grid_m

  use curvlinear_m
  use datasets_m
  use geometry_m
  use global_m
  use math_m
  use mesh_m
  use messages_m
  use functions_m
  use mesh_function_m
  use lib_oct_m
  use lib_oct_gsl_spline_m
  use lib_oct_parser_m
  use simul_box_m
  use specie_m
  use submesh_m

  implicit none
  
  private
  
  public ::              &
       double_grid_t,    &
       double_grid_init, &
       double_grid_end,  &
       double_grid_apply_local,     &
       double_grid_apply_glocal,    &
       double_grid_apply_non_local, &
       double_grid_get_rmax,        &
       double_grid_get_hmax

  type double_grid_t
     private
     FLOAT   :: h_fine(MAX_DIM)
     integer :: spacing_divisor
     integer :: interpolation_min
     integer :: interpolation_max
     integer :: nn
     logical :: use_double_grid
  end type double_grid_t

  FLOAT, parameter :: co (-4:5) = (/ &
       CNST(	3449600.0	)/CNST(	7142567040.0	),& !-4
       CNST(	-4484480.0	)/CNST(	793618560.0	),& !-3
       CNST(	6406400.0	)/CNST(	198404640.0	),& !-2
       CNST(	-11211200.0	)/CNST(	85030560.0	),& !-1
       CNST(	44844800.0	)/CNST(	56687040.0	),& ! 0
       CNST(	22422400.0	)/CNST(	56687040.0	),& ! 1
       CNST(	-8968960.0	)/CNST(	85030560.0	),& ! 2
       CNST(	5605600.0	)/CNST(	198404640.0	),& ! 3
       CNST(	-4076800.0	)/CNST(	793618560.0	),& ! 4
       CNST(	3203200.0	)/CNST(	7142567040.0	)/) ! 5

contains
  
  subroutine double_grid_init(this, m)
    type(double_grid_t), intent(out) :: this
    type(mesh_t),        intent(in)  :: m

    this%spacing_divisor = 3
    this%interpolation_min = -4
    this%interpolation_max = 5
    this%nn = (this%spacing_divisor - 1)/2
    this%h_fine(1:MAX_DIM) = m%h(1:MAX_DIM)/this%spacing_divisor
   
    !%Variable DoubleGrid
    !%Type logical
    !%Default no
    !%Section Mesh
    !%Description
    !% Enables or disables the use of a double grid technique to
    !% increase the precision of the application of the
    !% pseudopotentials.
    !%End

    call loct_parse_logical(check_inp('DoubleGrid'), .false., this%use_double_grid)

  end subroutine double_grid_init

  subroutine double_grid_end(this)
    type(double_grid_t), intent(inout) :: this
 
  end subroutine double_grid_end
  
  subroutine double_grid_apply_local(this, s, m, sm, x_atom, vl)
    type(double_grid_t), intent(in)  :: this
    type(specie_t),      intent(in)  :: s
    type(mesh_t),        intent(in)  :: m
    type(submesh_t),     intent(in)  :: sm
    FLOAT,               intent(in)  :: x_atom(1:MAX_DIM)
    FLOAT,               intent(out) :: vl(:)
    
    FLOAT :: r, vv
    FLOAT, allocatable :: vls(:)
    integer :: is
    integer :: ii, jj, kk
    integer :: start(1:3), pp, qq, rr, is2
    integer :: ll, mm, nn

    call push_sub('double_grid.double_grid_apply_local')

    if (.not. this%use_double_grid) then 

      do is = 1, sm%ns
        r = sqrt(sum((m%x(sm%jxyz(is), 1:3) - x_atom(1:3))**2))
        vl(is) = loct_splint(s%ps%vl, r)
      end do

    else

      ALLOCATE(vls(0:sm%ns), sm%ns)

      vls(0:sm%ns) = M_ZERO

      !for each grid point
      do is = 1, sm%ns

        !for each point in the fine mesh around the grid point
        do ii = -this%nn, this%nn
          do jj = -this%nn, this%nn
            do kk = -this%nn, this%nn

              r = sqrt(sum((m%x(sm%jxyz(is), 1:3) + this%h_fine(1:3) * (/ii, jj, kk/) - x_atom(1:3))**2))

              !calculate the potential
              vv = loct_splint(s%ps%vl, r)

              start(1:3) = m%Lxyz(sm%jxyz(is), 1:3) + this%interpolation_min * (/ii, jj, kk/)
              
              pp = start(1)
              do ll = this%interpolation_min, this%interpolation_max
                
                qq = start(2)
                do mm = this%interpolation_min, this%interpolation_max
                  
                  rr = start(3)
                  do nn = this%interpolation_min, this%interpolation_max
                    
                    is2 = sm%jxyz_inv(m%Lxyz_inv(pp, qq, rr))
                    vls(is2) = vls(is2) + co(ll)*co(mm)*co(nn)*vv
                    
                    rr = rr + kk
                  end do
                  qq = qq + jj
                end do
                pp = pp + ii
              end do

            end do
          end do
        end do
        
      end do

      vl(1:sm%ns) = vls(1:sm%ns)/(this%spacing_divisor**3)

      deallocate(vls)

    end if

    call pop_sub()
  end subroutine double_grid_apply_local

  subroutine double_grid_apply_glocal(this, s, m, sm, x_atom, dvl)
    type(double_grid_t),    intent(in)    :: this
    type(specie_t),         intent(in)    :: s
    type(mesh_t),           intent(in)    :: m
    type(submesh_t),        intent(in)    :: sm
    FLOAT,                  intent(in)    :: x_atom(1:MAX_DIM)
    FLOAT,                  intent(out)   :: dvl(:, :)
    
    FLOAT :: r, r2, x(1:3), vv(1:3)
    integer :: is, is2
    integer :: ii, jj, kk, ll, mm, nn
    integer :: start(1:3), pp, qq, rr

    FLOAT, allocatable :: dvs(:,:)


    if (.not. this%use_double_grid) then 

      do is = 1, sm%ns
        x(1:3) = m%x(sm%jxyz(is), 1:3) - x_atom(1:3)
        r = sqrt(sum(x(1:3)**2))
        if ( r > CNST(1e-5) ) then 
          dvl(is, 1:3) = -loct_splint(s%ps%dvl, r)*x(1:3)/r
        end if
      end do
 
    else

      ALLOCATE(dvs(0:sm%ns, 1:3), 3*(sm%ns+1))

      dvs(0:sm%ns, 1:3) = M_ZERO

      !for each grid point
      do is = 1, sm%ns
        
        do ii = -this%nn, this%nn
          do jj = -this%nn, this%nn
            do kk = -this%nn, this%nn
              
              x(1:3) = m%x(sm%jxyz(is), 1:3) + this%h_fine(1:3) * (/ii, jj, kk/) - x_atom(1:3)

              r2 = sum(x(1:3)**2)

              if ( r > CNST(1e-7)) then 

                vv(1:3) = -loct_splint(s%ps%dvl, sqrt(r2))*x(1:3)/sqrt(r2)

                start(1:3) = m%Lxyz(sm%jxyz(is), 1:3) + this%interpolation_min * (/ii, jj, kk/)

                pp = start(1)
                do ll = this%interpolation_min, this%interpolation_max
                  
                  qq = start(2)
                  do mm = this%interpolation_min, this%interpolation_max

                    rr = start(3)
                    do nn = this%interpolation_min, this%interpolation_max

                      is2 = sm%jxyz_inv(m%Lxyz_inv(pp, qq, rr))
                      dvs(is2, 1:3) = dvs(is2, 1:3)  + co(ll)*co(mm)*co(nn)*vv(1:3)
                      
                      rr = rr + kk
                    end do
                    qq = qq + jj
                  end do
                  pp = pp + ii
               end do
               
             end if

            end do !kk
          end do !jj
        end do !ii

      end do !is

      dvl(1:sm%ns, 1:3) = dvs(1:sm%ns, 1:3)/(this%spacing_divisor**3)
      
      deallocate(dvs)

    end if

  end subroutine double_grid_apply_glocal

  subroutine double_grid_apply_non_local(this, s, m, x_atom, sm, l, lm, ic, vnl, dvnl)
    type(double_grid_t),    intent(in)  :: this
    type(specie_t), target, intent(in)  :: s
    type(mesh_t),           intent(in)  :: m
    FLOAT,                  intent(in)  :: x_atom(1:MAX_DIM)
    type(submesh_t),        intent(in)  :: sm
    integer,                intent(in)  :: l
    integer,                intent(in)  :: lm
    integer,                intent(in)  :: ic
    FLOAT,                  intent(out) :: vnl(:)
    FLOAT,                  intent(out) :: dvnl(:, :)
    
    FLOAT :: x(MAX_DIM), vv, dvv(1:MAX_DIM)
    integer :: is, ii, jj, kk, idir

    integer :: start(1:3), pp, qq, rr, is2
    integer :: ll, mm, nn
    
    FLOAT, allocatable :: vls(:), dvls(:,:)

    call push_sub('double_grid.double_grid_apply_local')

    if(.not. this%use_double_grid) then 

      do is = 1, sm%ns
        call specie_real_nl_projector(s, x_atom, m%x(sm%jxyz(is), :), l, lm, ic, vnl(is), dvnl(is, :))
      end do

    else
      
      ALLOCATE(vls(0:sm%ns), sm%ns+1)
      ALLOCATE(dvls(0:sm%ns, 1:3), 3*(sm%ns+1))

      vls(0:sm%ns) = M_ZERO
      dvls(0:sm%ns, 1:3) = M_ZERO
      
      do is = 1, sm%ns

        do ii = -this%nn, this%nn
          do jj = -this%nn, this%nn
            do kk = -this%nn, this%nn
              
              x(1:3) = m%x(sm%jxyz(is), 1:3) + this%h_fine(1:3) * (/ii, jj, kk/) 

              call specie_real_nl_projector(s, x_atom, x, l, lm, ic, vv, dvv(1:3))
              
              start(1:3) = m%Lxyz(sm%jxyz(is), 1:3) + this%interpolation_min * (/ii, jj, kk/)
              
              pp = start(1)
              do ll = this%interpolation_min, this%interpolation_max
                
                qq = start(2)
                do mm = this%interpolation_min, this%interpolation_max
                  
                  rr = start(3)
                  do nn = this%interpolation_min, this%interpolation_max
                    
                    is2 = sm%jxyz_inv(m%Lxyz_inv(pp, qq, rr))

                    vls(is2) = vls(is2) + co(ll)*co(mm)*co(nn)*vv
                    dvls(is2, 1:3) = dvls(is2, 1:3) + co(ll)*co(mm)*co(nn)*dvv(1:3)
                    
                    rr = rr + kk
                  end do
                  qq = qq + jj
                end do
                pp = pp + ii
              end do
              
            end do
          end do
        end do

      end do

      vnl(1:sm%ns) = vls(1:sm%ns)/(this%spacing_divisor**3)
      dvnl(1:sm%ns, 1:3) = dvls(1:sm%ns, 1:3)/(this%spacing_divisor**3)

      deallocate(vls, dvls)

    end if

    call pop_sub()

  end subroutine double_grid_apply_non_local

  FLOAT function double_grid_get_hmax(this, mesh) result(hmax)
    type(double_grid_t), intent(in) :: this
    type(mesh_t),        intent(in) :: mesh

    if(this%use_double_grid) then 
      hmax = maxval(this%h_fine(1:MAX_DIM))
    else
      hmax = maxval(mesh%h(1:MAX_DIM))
    end if
    
  end function double_grid_get_hmax

  FLOAT pure function double_grid_get_rmax(this, s, m) result(rmax)
    type(double_grid_t),     intent(in) :: this
    type(specie_t), target,  intent(in) :: s
    type(mesh_t),            intent(in) :: m
    
    rmax = loct_spline_cutoff_radius(s%ps%dvl, s%ps%projectors_sphere_threshold)

    if(this%use_double_grid) then 
      rmax = rmax + this%interpolation_max * maxval(m%h(1:3))
    end if
  end function double_grid_get_rmax

end module double_grid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
