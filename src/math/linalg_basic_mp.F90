!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: linal_basic_mp.F90 3154 2007-09-01 15:14:11Z xavier $

real(8) function nrm2_mp_r(n1, vec1) result(nrm2)
  integer, intent(in) :: n1
  real(4), intent(in) :: vec1(:)
  
  integer :: ip
  real(8) :: tmp
  nrm2 = 0.0_8
  
  do ip = 1, n1
    tmp = real(vec1(ip), 8)
    nrm2 = nrm2 + tmp**2
  end do
  
  nrm2 = sqrt(nrm2)
  
end function nrm2_mp_r

real(8) function nrm2_mp_c(n1, vec1) result(nrm2)
  integer,    intent(in) :: n1
  complex(4), intent(in) :: vec1(:)
  
  integer :: ip
  real(8) :: re
  real(8) :: im
  nrm2 = 0.0_8
  
  do ip = 1, n1
    re = real(vec1(ip), 8)
    im = real(aimag(vec1(ip)), 8)
    nrm2 = nrm2 + re**2 + im**2
  end do
  
  nrm2 = sqrt(nrm2)
  
end function nrm2_mp_c

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
