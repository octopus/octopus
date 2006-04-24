!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module specie_grid_m
  use global_m
  use string_m
  use mpi_m
  use messages_m
  use datasets_m
  use lib_oct_gsl_spline_m
  use io_m
  use units_m
  use atomic_m
  use ps_m
  use math_m
  use specie_m
  use mesh_m
  use mesh_function_m

  implicit none

  private
  public ::                     &
       specie_get_density

contains

  ! ---------------------------------------------------------
  subroutine specie_get_density(s, pos, m, rho)
    type(specie_t), intent(in) :: s
    FLOAT,          intent(in) :: pos(MAX_DIM)
    type(mesh_t),   intent(in) :: m
    FLOAT,          intent(out) :: rho(:)

    FLOAT :: d, dmin
    integer :: i, imin

    call push_sub('specie.specie_get_density')

    select case(s%type)

    case(SPEC_ALL_E)

      dmin=M_ZERO

      !find the point of the grid that is closer to the atom
      do i=1,m%np
        d = sum( ( pos(1:MAX_DIM)-m%x(i,1:MAX_DIM) )**2 )

        if ( ( d < dmin ) .or. ( i == 1 ) ) then 
          imin = i
          dmin = d 

        end if

      end do

      if (dmin > CNST(1e-5)) then 

        write(message(1), '(a,f12.2,a)') "Atom displaced ", sqrt(dmin), " [b]"
        write(message(2), '(a,3f12.2)') "Original position ", pos
        write(message(3), '(a,3f12.2)') "Displaced position ", m%x(imin,:) 

        call write_warning(3)

      endif

      rho(1:m%np) = M_ZERO
      rho(imin) = -s%Z/m%vol_pp(imin)


    end select

    call pop_sub()

  end subroutine specie_get_density

end module specie_grid_m
