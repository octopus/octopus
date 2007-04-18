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
       dg_add_localization_density, &
       dg_get_density_correction,   &
       dg_get_potential_correction, &
       dg_filter_potential,         &
       dg_get_hmax

  type double_grid_t
     private
     FLOAT   :: h_fine(MAX_DIM)
     integer :: spacing_divisor
     integer :: interpolation_order
     integer :: interpolation_min
     integer :: interpolation_max
     integer :: nn
     integer :: loc_function
     logical :: use_double_grid
     logical :: filter
  end type double_grid_t

  FLOAT,   parameter :: rmax = CNST(30.0), sigma = CNST(0.625)
  integer, parameter :: npoints = 3000
  
  integer, parameter :: &
       F_NONE = 0,      &
       F_GAUSSIAN = 1 
  
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

  FLOAT, allocatable :: cc(:,:,:)

contains
  
  subroutine double_grid_init(this, m)
    type(double_grid_t), intent(out) :: this
    type(mesh_t),        intent(in)  :: m

    integer :: ii, jj, kk

    this%spacing_divisor = 3
    this%interpolation_order = 9
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

    !%Variable DoubleGridFilter
    !%Type logical
    !%Default no
    !%Section Mesh
    !%Description
    !% When this is enabled a filter is applied to the pseudopotential.
    !%End

    call loct_parse_logical(check_inp('DoubleGridFilter'), .false., this%filter)

    !%Variable LocalizationDensity
    !%Type integer
    !%Default none
    !%Section Mesh
    !%Description
    !% When a localization density is used, the local part of the
    !% pseudopotential is separated in a localized part and a long
    !% range term. The long range term is calculated by solving the
    !% poisson equation of a localized charge density. This separation
    !% is required to apply to use the double grid technique or to
    !% apply a filter to the local part of the pseudopotential. This
    !% variable selects which kind of function is used for the
    !% localized charge density. By default this technique is not
    !% used.
    !%Option none 0 
    !% The potential is not separated. This is the default.
    !%Option gaussian 3
    !% A gaussian charge distribution is used, the correspoding
    !% potential is erf(r)/r.
    !%End

    call loct_parse_int(check_inp('LocalizationDensity'), F_GAUSSIAN, this%loc_function)

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
    type(loct_spline_t)  :: ft
    FLOAT :: r, vv
    integer :: ip, ip2
    integer :: ii, jj, kk, idx(1:3)
    integer :: ll, mm, nn

    if(dg_add_localization_density(this)) then 
      ps_spline => s%ps%vll
    else
      ps_spline => s%ps%vl
    end if

    if(this%use_double_grid ) then 
      write(message(1),'(a)')         'Info: Integral of the local part of the pseudopotential'
      call write_info(1)
    end if

    if(this%use_double_grid) then
      call loct_spline_init(ft)
      call loct_spline_3dft(ps_spline, ft, CNST(1.0))
      write(message(1),'(a, f14.6)')  '      exact       = ', loct_splint(ft, CNST(0.0))
      call write_info(1)
      call loct_spline_end(ft)  
    end if

    do ip = 1, m%np
      r = sqrt(sum( (m%x(ip, :) - x_atom(:))**2 ))
      vl(ip) = loct_splint(ps_spline, r)
!      write(12,*) r, vl(ip)
    end do

    if (this%use_double_grid) then 

      write(message(1),'(a, f14.6)')  '      coarse grid = ', dmf_integrate(m, vl)
      call write_info(1)

      vl(1:m%np) = M_ZERO

      !for each grid point
      do ip = 1, m%np

        !for each point in the fine mesh around the grid point
        do ii = -this%nn, this%nn
          do jj = -this%nn, this%nn
            do kk = -this%nn, this%nn
              
              idx = (/ ii, jj, kk /)

              !calculate the potential
              r = sqrt(sum( (m%x(ip, 1:3) + this%h_fine(1:3)*idx(1:3) - x_atom(1:3))**2 ))
              vv = loct_splint(ps_spline, r) 

              !add thsi value to all the points that depend on this point
              if ( r < specie_local_cutoff_radius(s) ) then 
                do ll = this%interpolation_min, this%interpolation_max
                  do mm = this%interpolation_min, this%interpolation_max
                    do nn = this%interpolation_min, this%interpolation_max
                      
                      idx = m%Lxyz(ip, 1:3) + (/ ii, jj, kk /)*(/ ll , mm, nn /)
                      ip2 = m%Lxyz_inv(idx(1), idx(2), idx(3))
                      vl(ip2) = vl(ip2) + co(ll)*co(mm)*co(nn)*vv
                      
                    end do
                  end do
                end do
                
              end if

            end do
          end do
        end do
        
      end do

      !divide by the factor between the spacing of the coarse and the fine grid
      vl(1:m%np) = vl(1:m%np)/(this%spacing_divisor**3)
      write(message(1),'(a, f14.6)')  '      double grid = ', dmf_integrate(m, vl)
      call write_info(1)

!      do ip = 1, m%np
!        write(13,*) vl(ip)
!      end do
!      stop

    end if

  end subroutine double_grid_apply_local

  subroutine double_grid_apply_glocal(this, s, m, x_atom, dvl)
    type(double_grid_t),    intent(in)    :: this
    type(specie_t), target, intent(in)    :: s
     type(mesh_t),          intent(in)    :: m
    FLOAT,                  intent(in)    :: x_atom(1:MAX_DIM)
    FLOAT,                  intent(out)   :: dvl(:, :)
    
    type(loct_spline_t), pointer  :: ps_spline
    FLOAT :: r, x(1:3), vv(1:3)
    integer :: ip
    integer :: ii, jj, kk, idx(1:3)
    integer :: iw, ipnb, idx_ip(1:3)
    integer :: ll, mm, nn, ip2
    
    if(dg_add_localization_density(this)) then 
      ps_spline => s%ps%dvll
    else
      ps_spline => s%ps%dvl
    end if

    if (.not. this%use_double_grid) then 
      
      do ip = 1, m%np
        x(1:3) = m%x(ip, 1:3) - x_atom(1:3)
        r = sqrt(sum(x(1:3)**2))
        if ( r > CNST(1e-7) ) then 
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
              
              idx = (/ ii, jj, kk /)

              x(1:3) = m%x(ip, 1:3) + this%h_fine(1:3)*idx(1:3) - x_atom(1:3)

              
              r = sqrt(sum(x(1:3)**2))
              if ( r > CNST(1e-7) .and. r < specie_local_cutoff_radius(s) ) then 

                vv(1:3) = -loct_splint(ps_spline, r)*x(1:3)/r
                
                do ll = this%interpolation_min, this%interpolation_max
                  do mm = this%interpolation_min, this%interpolation_max
                    do nn = this%interpolation_min, this%interpolation_max
                      
                      idx = m%Lxyz(ip, 1:3) + (/ ii, jj, kk /)*(/ ll , mm, nn /)
                      ip2 = m%Lxyz_inv(idx(1), idx(2), idx(3))
                      dvl(ip2, 1:3) = dvl(ip2, 1:3) + co(ll)*co(mm)*co(nn)*vv(1:3)
                      
                    end do
                  end do
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
    integer :: is, is2, idx_ip(1:MAX_DIM)
    integer :: ii, jj, kk, idx(1:3)
    integer :: ll, mm, nn

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
              
              idx = (/ii, jj, kk/)
              x(1:3) = m%x(jxyz(is), 1:3) + this%h_fine(1:3) * idx(1:3)
              call specie_real_nl_projector(s, x_atom, x, l, lm, ic, vv, dvv(1:3))

              do ll = this%interpolation_min, this%interpolation_max
                do mm = this%interpolation_min, this%interpolation_max
                  do nn = this%interpolation_min, this%interpolation_max
                    
                    idx = m%Lxyz(jxyz(is), 1:3) + (/ ii, jj, kk /)*(/ ll , mm, nn /)
                    is2 = jxyz_inv(m%Lxyz_inv(idx(1), idx(2), idx(3)))

                    if(is2 > 0 ) then 
                      vnl(is2) = vnl(is2) + co(ll)*co(mm)*co(nn)*vv
                      dvnl(is2, 1:3) = dvnl(is2, 1:3) + co(ll)*co(mm)*co(nn)*dvv(1:3)
                    end if

                  end do
                end do
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

  logical function dg_add_localization_density(this) 
    type(double_grid_t), intent(in) :: this
    
    dg_add_localization_density = ( this%loc_function /= F_NONE )
    
  end function dg_add_localization_density


  subroutine dg_get_potential_correction(this, vlc)
    type(double_grid_t), intent(in) :: this
    type(loct_spline_t), intent(out) :: vlc

    FLOAT :: dr, r(1:npoints), pot(1:npoints)
    integer :: ii
    
    dr = rmax/(npoints-M_ONE)
        
    r(1) = M_ZERO

    pot(1) = M_TWO/(sqrt(M_TWO*M_PI)*sigma)

    do ii = 1, npoints
      if ( ii /= 1 ) pot(ii) = loct_erf(r(ii)/(sigma*sqrt(M_TWO)))/r(ii)
      if ( ii /= npoints ) r(ii+1) = r(ii) + dr
    end do
    
    call loct_spline_fit(npoints, r, pot, vlc)

  end subroutine dg_get_potential_correction

  subroutine dg_get_density_correction(this, rc)
    type(double_grid_t), intent(in) :: this
    type(loct_spline_t), intent(out) :: rc

    integer :: ii
    FLOAT :: dr, r(1:npoints), rho(1:npoints)

    dr = rmax/(npoints-M_ONE)
    
    r(1) = M_ZERO
    do ii = 1, npoints
      rho(ii) = M_ONE/(sigma*sqrt(M_TWO*M_PI))**3*exp(-M_HALF*r(ii)**2/sigma**2)
      if ( ii /= npoints ) r(ii+1) = r(ii) + dr
    end do
    
    call loct_spline_fit(npoints, r, rho, rc)

  end subroutine dg_get_density_correction


  logical function dg_filter_potential(this)
    type(double_grid_t), intent(in) :: this

    dg_filter_potential = this%filter
  end function dg_filter_potential

  FLOAT function dg_get_hmax(this, mesh) result(hmax)
    type(double_grid_t), intent(in) :: this
    type(mesh_t),        intent(in) :: mesh

    if(this%use_double_grid) then 
      hmax = maxval(this%h_fine(1:MAX_DIM))
    else
      hmax = maxval(mesh%h(1:MAX_DIM))
    end if
    
  end function dg_get_hmax

end module double_grid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
