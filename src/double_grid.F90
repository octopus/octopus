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
  use specie_m

  implicit none
  
  private
  
  public ::              &
       double_grid_t,    &
       double_grid_init, &
       double_grid_end,  &
       double_grid_apply_local,     &
       double_grid_apply_glocal,    &
       double_grid_apply_non_local, &
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
  
  subroutine double_grid_apply_local(this, s, m, x_atom, vl)
    type(double_grid_t), intent(in)  :: this
    type(specie_t), target, intent(in)  :: s
    type(mesh_t),        intent(in)  :: m
    FLOAT,               intent(in)  :: x_atom(1:MAX_DIM)
    FLOAT,               intent(out) :: vl(:)
    
    type(loct_spline_t), pointer  :: ps_spline

    FLOAT :: r2, vv
    integer :: ip, ip2
    integer :: ii, jj, kk, ll, mm, nn
    integer :: start(1:3), pp, qq, rr

    ps_spline => s%ps%vl

    if ( .not. this%use_double_grid) then 

      call dmf_put_radial_spline(m, s%ps%vl, x_atom, vl)
      
    else

      vl(1:m%np) = M_ZERO

      !for each grid point
      do ip = 1, m%np

        !for each point in the fine mesh around the grid point
        do ii = -this%nn, this%nn
          do jj = -this%nn, this%nn
            do kk = -this%nn, this%nn
              
              r2 = sum( (m%x(ip, 1:3) + this%h_fine(1:3)*(/ ii, jj, kk /) - x_atom(1:3))**2 )

              if ( r2 < specie_local_cutoff_radius(s)**2 ) then 

                !calculate the potential
                vv = loct_splint(ps_spline, sqrt(r2)) 
                
                !add the value of the potential to all the points that depend on this point

                start(1:3) = m%Lxyz(ip, 1:3) + this%interpolation_min * (/ ii, jj, kk /)

                pp = start(1)
                do ll = this%interpolation_min, this%interpolation_max

                  qq = start(2)
                  do mm = this%interpolation_min, this%interpolation_max

                    rr = start(3)
                    do nn = this%interpolation_min, this%interpolation_max

                      ! here: (/ pp, qq, rr /) = m%Lxyz(ip, 1:3) + (/ ll*ii, mm*jj, nn*kk/)
                      ip2 = m%Lxyz_inv(pp, qq, rr)
                      if(ip2 > 0 .and. ip2 <= m%np) &
                           vl(ip2) = vl(ip2) + co(nn)*co(mm)*co(ll)*vv

                      rr = rr + kk
                    end do
                    qq = qq + jj
                  end do
                 pp = pp + ii
                end do
                
              end if

            end do
          end do
        end do
        
      end do

      vl(1:m%np) = vl(1:m%np)/(this%spacing_divisor**3)

    end if

  end subroutine double_grid_apply_local

  subroutine double_grid_apply_glocal(this, s, m, x_atom, dvl)
    type(double_grid_t),    intent(in)    :: this
    type(specie_t), target, intent(in)    :: s
     type(mesh_t),          intent(in)    :: m
    FLOAT,                  intent(in)    :: x_atom(1:MAX_DIM)
    FLOAT,                  intent(out)   :: dvl(:, :)
    
    type(loct_spline_t), pointer  :: ps_spline
    FLOAT :: r, r2, x(1:3), vv(1:3)
    integer :: ip, ip2
    integer :: ii, jj, kk, ll, mm, nn
    integer :: start(1:3), pp, qq, rr

    ps_spline => s%ps%dvl

    if (.not. this%use_double_grid) then 
      
      do ip = 1, m%np
        x(1:3) = m%x(ip, 1:3) - x_atom(1:3)
        r = sqrt(sum(x(1:3)**2))
        if ( r > CNST(1e-5) ) then 
          dvl(ip, 1:3) = -loct_splint(ps_spline, r)*x(1:3)/r
        else 
          dvl(ip, 1:3) = M_ZERO
        end if
      end do
 
    else

      dvl(1:m%np, 1:3) = M_ZERO

      do ip = 1, m%np
        
        do ii = -this%nn, this%nn
          do jj = -this%nn, this%nn
            do kk = -this%nn, this%nn
              
              x(1:3) = m%x(ip, 1:3) + this%h_fine(1:3)*(/ ii, jj, kk /) - x_atom(1:3)
              
              r2 = sum(x(1:3)**2)

              if ( r2 > CNST(1e-10) .and. r2 < specie_local_cutoff_radius(s)**2 ) then 

                vv(1:3) = -loct_splint(ps_spline, sqrt(r2))*x(1:3)/sqrt(r2)

                start(1:3) = m%Lxyz(ip, 1:3) + this%interpolation_min * (/ ii, jj, kk /)

                pp = start(1)
                do ll = this%interpolation_min, this%interpolation_max
                  
                  qq = start(2)
                  do mm = this%interpolation_min, this%interpolation_max

                    rr = start(3)
                    do nn = this%interpolation_min, this%interpolation_max

                      ip2 = m%Lxyz_inv(pp, qq, rr)
                      if(ip2 > 0 .and. ip2 <= m%np) then 
                        dvl(ip2, 1:3) = dvl(ip2, 1:3) + co(ll)*co(mm)*co(nn)*vv(1:3)
                      end if

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

      end do !ip

      dvl(1:m%np, 1:3) = dvl(1:m%np, 1:3)/(this%spacing_divisor**3)

    end if

  end subroutine double_grid_apply_glocal

  subroutine double_grid_apply_non_local(this, s, m, x_atom, ns, jxyz, l, lm, ic, vnl, dvnl)
    type(double_grid_t),    intent(in)  :: this
    type(specie_t), target, intent(in)  :: s
    type(mesh_t),           intent(in)  :: m
    FLOAT,                  intent(in)  :: x_atom(1:MAX_DIM)
    integer,                intent(in)  :: ns
    integer,                intent(in)  :: jxyz(:)
    integer,                intent(in)  :: l
    integer,                intent(in)  :: lm
    integer,                intent(in)  :: ic
    FLOAT,                  intent(out) :: vnl(:)
    FLOAT,                  intent(out) :: dvnl(:, :)
    
    FLOAT, allocatable :: jxyz_inv(:)

    FLOAT :: x(MAX_DIM), vv, dvv(1:MAX_DIM)
    integer :: is, is2
    integer :: ii, jj, kk, ll, mm, nn
    integer :: start(1:3), pp, qq, rr

    if(.not. this%use_double_grid) then 

      do is = 1, ns
        call specie_real_nl_projector(s, x_atom, m%x(jxyz(is), :), l, lm, ic, vnl(is), dvnl(is, :))
      end do

    else

      !build the inverse of jxyz, points not in the sphere go to 0
      ALLOCATE(jxyz_inv(1:m%np), m%np)
      jxyz_inv(1:m%np) = 0
      do is = 1, ns
       jxyz_inv(jxyz(is)) = is
      end do

      vnl(1:ns) = M_ZERO
      dvnl(1:ns, 1:3) = M_ZERO

      do is = 1, ns

        do ii = -this%nn, this%nn
          do jj = -this%nn, this%nn
            do kk = -this%nn, this%nn
              
              x(1:3) = m%x(jxyz(is), 1:3) + this%h_fine(1:3) * (/ii, jj, kk/) 

              call specie_real_nl_projector(s, x_atom, x, l, lm, ic, vv, dvv(1:3))
              
              start(1:3) = m%Lxyz(jxyz(is), 1:3) + this%interpolation_min * (/ ii, jj, kk /)
              
              pp = start(1)
              do ll = this%interpolation_min, this%interpolation_max
                
                qq = start(2)
                do mm = this%interpolation_min, this%interpolation_max
                  
                  rr = start(3)
                  do nn = this%interpolation_min, this%interpolation_max

                    is2 = jxyz_inv(m%Lxyz_inv(pp, qq, rr))

                    if(is2 /= 0 ) then 
                      vnl(is2) = vnl(is2) + co(ll)*co(mm)*co(nn)*vv
                      dvnl(is2, 1:3) = dvnl(is2, 1:3) + co(ll)*co(mm)*co(nn)*dvv(1:3)
                    end if

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

      deallocate(jxyz_inv)
      
      vnl(1:ns) = vnl(1:ns)/(this%spacing_divisor**3)
      dvnl(1:ns, 1:3) = dvnl(1:ns, 1:3)/(this%spacing_divisor**3)

    end if

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

end module double_grid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
