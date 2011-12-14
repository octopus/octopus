!! Copyright (C) 2007 Xavier Andrade
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
!! $Id: help.F90 $

#include "global.h"

program photoelectron_spectrum
  
  use command_line_m
  use datasets_m
  use global_m
  use io_m
  use messages_m
  use parser_m
  use pes_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  integer              :: argc, ierr, mode, interp

  integer              :: dim, ll(MAX_DIM), ii
  FLOAT                :: Emax, Emin,Estep, uEstep,uEspan(2), pol(3)
  FLOAT, pointer       :: lk(:),RR(:)
  FLOAT, allocatable   :: PESK(:,:,:)
  logical              :: interpolate
  type(block_t) :: blk
  
  character(len=80) :: filename

  !Initial values
  ll = 1 
  mode = 1
  interpolate = .true. 

  call global_init()
  call parser_init()
  if(parse_block('CalculationMode', blk) == 0) then
    call datasets_init(3, blk)
  else
    call datasets_init(3)
  end if
  call io_init()
  call io_init_datasets()

  call getopt_init(ierr)
  if(ierr.ne.0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-photoelectron-spectrum command is not available."
    call messages_fatal(2)
  end if
  
  call getopt_photoelectron_spectrum(mode,interp,uEstep, uEspan,pol)
  if(interp .eq. 0) interpolate = .false.

  call PES_mask_read_info(tmpdir, dim, Emax, Estep, ll(1), Lk,RR)

  write(message(1), '(a)') 'Read PES info file.'
  call messages_info(1)

  do ii=2, dim
    ll(ii) = ll(1)
  end do    
  
  SAFE_ALLOCATE(PESK(1:ll(1),1:ll(2),1:ll(3)))

  filename=io_workpath('td.general/PESM_map.obf')
  call io_binary_read(trim(filename),ll(1)**dim,PESK, ierr) 
  if(ierr > 0) then
    message(1) = "Failed to read file "//trim(filename)
    call messages_fatal(1)
  end if


  write(message(1), '(a)') 'Read PES restart file.'
  call messages_info(1)

  !! set user values
  if(uEstep >  0 .and. uEstep > Estep)    Estep = uEstep
  
  if(uEspan(1) > 0 ) Emin = uEspan(1)

  if(uEspan(2) > 0 ) Emax = uEspan(2)


  call unit_system_init()
 

  select case(mode)
  case(1) ! Energy-resolved
    write(message(1), '(a)') 'Calculating energy-resolved PES'
    call messages_info(1)
    call PES_mask_dump_power_totalM(PESK,'td.general/PES_power.sum', Lk, dim, Emax, Estep, interpolate)
 
 
  case(2) ! Angle and energy resolved
    write(message(1), '(a)') 'Not yet implemented'
    call messages_warning(1)

  case(3) ! On a plane
    write(message(1), '(a)') 'Calculating momentum-resolved PES on plane z=0'
    call messages_info(1)
    call PES_mask_dump_full_mapM(PESK, 'td.general/PES_map.z=0', Lk, dim, dir = 3)    

  case(4) ! Angle energy resolved on plane 
    write(message(1), '(a)') 'Calculating angle and energy-resolved PES (cartesian) with polarization axis (1,0,0)'
    call messages_info(1)
    if(uEstep >  0 .and. uEstep > Estep) then
      Estep = uEstep
    else
      Estep = Emax/size(Lk,1)
    end if

    call PES_mask_dump_ar_plane_M(PESK,'td.general/PES_energy_map.sum', Lk, dim, Emax, Estep)

  end select


  write(message(1), '(a)') 'Done'
  call messages_info(1)

  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()
  
  SAFE_DEALLOCATE_A(PESK)    

end program photoelectron_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
