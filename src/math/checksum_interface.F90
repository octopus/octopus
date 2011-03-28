!! Copyright (C) 2010 X. Andrade
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
!! $Id: checksum_interface.F90 4852 2009-01-08 23:31:14Z dstrubbe $

module checksum_m
  
  public ::                 &
    checksum_calculate

  interface
    subroutine checksum_calculate(algorithm, narray, array, checksum)
      integer,    intent(in)  :: algorithm
      integer,    intent(in)  :: narray
      integer,    intent(in)  :: array
      integer(8), intent(out) :: checksum
    end subroutine checksum_calculate
  end interface

end module checksum_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
