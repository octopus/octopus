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

#include "config_F90.h"

program make_st
  use global
  use liboct
  use states
  use system

  implicit none

  type(system_type) :: sys
  character(len=100) :: str
  integer :: i, n, ik, ist, idim, type

  ! Initialize stuff
  call global_init()
  call units_init()
  call system_init(sys)

#ifdef T_REAL
  deallocate(sys%st%dpsi)
#endif

  allocate(sys%st%zpsi (0:sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik), &
       sys%st%eigenval(sys%st%nst, sys%st%nik))
  
  if(.not.zstates_load_restart ("tmp/restart.static", sys%m, sys%st)) then
    message(1) = "Error opening 'restart.static' file"
    call write_fatal(1)
  endif

  str = C_string("MakeStates")
  n = oct_parse_block_n(str)
  do i = 1, n
    call oct_parse_block_int(str, i-1, 0, ik)
    call oct_parse_block_int(str, i-1, 1, ist)
    call oct_parse_block_int(str, i-1, 2, idim)
    call oct_parse_block_int(str, i-1, 3, type)
    select case(type)
    case(1)
      call wf_gaussian(i-1)
    end select
  end do

  ! renormalize functions
  call wf_renormalize()

  ! save wfs in a new static file
  call zstates_write_restart("tmp/restart.static.new", sys%m, sys%st)

contains
  subroutine wf_gaussian(line)
    integer, intent(in) :: line
    
    real(r8) :: x1(3), x(3), s, k(3)
    integer :: i, j
    
    ! read gaussian parameters
    x1 = M_ZERO; k = M_ZERO
    call oct_parse_block_double(str, line, 4,  s)
    s = s * units_inp%length%factor

    j = 5
    do i = 1, conf%dim
      call oct_parse_block_double(str, line, j,  x1(i))
      x1(i) = x1(i) * units_inp%length%factor
      j = j + 1
    end do

    do i = 1, conf%dim
      call oct_parse_block_double(str, line, j,  k(i))
      k(i) = k(i) / units_inp%length%factor
      j = j + 1
    end do
    
    ! build a gaussian
    do i = 1, sys%m%np
      call mesh_xyz(sys%m, i, x(1:conf%dim))
      sys%st%zpsi(i, idim, ist, ik) = exp(-sum((x(1:conf%dim)-x1(1:conf%dim))**2)/(2*s*s) + &
           M_zI*sum(k(1:conf%dim)*(x(1:conf%dim)-x1(1:conf%dim))))
    end do
    
  end subroutine wf_gaussian
  
  subroutine wf_renormalize()
    integer :: ik, ist
    real(r8) :: nrm2
    
    do ik = 1, sys%st%nik
      do ist = 1, sys%st%nst
        nrm2 = sum(abs(sys%st%zpsi(1:sys%m%np, 1:sys%st%dim, ist, ik))**2)*sys%m%vol_pp
        sys%st%zpsi(1:sys%m%np, 1:sys%st%dim, ist, ik) = sys%st%zpsi(1:sys%m%np, 1:sys%st%dim, ist, ik)/ sqrt(nrm2)
      end do
    end do
  end subroutine wf_renormalize

end program make_st
