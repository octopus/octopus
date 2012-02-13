!! Copyright (C) 2011 U. De Giovannini
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
  use utils_m
  use varinfo_m

  implicit none

  integer              :: argc, ierr, mode, interp

  integer              :: dim, ll(MAX_DIM), ii, dir
  FLOAT                :: Emax, Emin,Estep, uEstep,uEspan(2), pol(3)
  FLOAT                :: uThstep,uThspan(2),uPhstep,uPhspan(2) 
  FLOAT                :: ThBatch(3),PhBatch(3), center(3)
  FLOAT, pointer       :: lk(:),RR(:)
  FLOAT, allocatable   :: PESK(:,:,:)
  logical              :: interpol
  type(block_t)        :: blk
  
  character(len=80) :: filename

  !Initial values
  ll = 1 
  mode = 1
  interpol = .true. 

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

  !set default values
  mode = 1
  interp = 1
  uEstep = -1
  uEspan = (/-1,-1/)
  uThstep = -1
  uThspan = (/-1,-1/)
  uPhstep = -1
  uPhspan = (/-1,-1/)
  center = (/0,0,0/)
  
  call get_laser_polarizaion(pol)
  
  call getopt_photoelectron_spectrum(mode,interp,uEstep, uEspan,&
                                     uThstep, uThspan, uPhstep, uPhspan, pol, center)
  if(interp .eq. 0) interpol = .false.

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
 
  write(message(1),'(a,f10.2,a2,f10.2,a2,f10.2,a1)') &
                   "Zenith axis: (",pol(1),", ",pol(2),", ",pol(3),")"
  call messages_info(1)


  ! choose what to calculate
  ! these functions are defined in pes_mask_out_inc.F90
  select case(mode)
  case(1) ! Energy-resolved
    write(message(1), '(a)') 'Calculating energy-resolved PES'
    call messages_info(1)
    call PES_mask_dump_power_totalM(PESK,'td.general/PES_power.sum', Lk, dim, Emax, Estep, interpol)
 
 
  case(2) ! Angle and energy resolved
    write(message(1), '(a)') 'Calculating angle- and energy-resolved PES'
    call messages_info(1)
    call PES_mask_dump_ar_polar_M(PESK,'td.general/PES_angle_energy.map', Lk, dim, pol, Emax, Estep)


  case(3) ! On a plane
    
    dir = -1
    if(sum((pol-(/1 ,0 ,0/))**2)  <= 1E-14  )  dir = 1
    if(sum((pol-(/0 ,1 ,0/))**2)  <= 1E-14  )  dir = 2
    if(sum((pol-(/0 ,0 ,1/))**2)  <= 1E-14  )  dir = 3

    filename = "td.general/PES_velocity.map."//index2axis(dir)//"=0"


    if (dir == -1) then
        write(message(1), '(a)') 'Unrecognized plane. Use -V to change.'
        call messages_fatal(1)
      else
        write(message(1), '(a)') 'Calculating velocity map on a plane '//index2axis(dir)//"=0"
        call messages_info(1)
    end if 
    
    call PES_mask_dump_full_mapM_cut(PESK, filename, Lk, dim, dir)    

  case(4) ! Angle energy resolved on plane 
    write(message(1), '(a)') 'Calculating angle and energy-resolved PES'
    call messages_info(1)
    if(uEstep >  0 .and. uEstep > Estep) then
      Estep = uEstep
    else
      Estep = Emax/size(Lk,1)
    end if

    call PES_mask_dump_ar_plane_M(PESK,'td.general/PES_energy.map', Lk, dim, pol, Emax, Estep)

  case(5) ! Full momentum resolved matrix 
    write(message(1), '(a)') 'Calculating full momentum-resolved PES'
    call messages_info(1)
    call PES_mask_dump_full_mapM(PESK, 'td.general/PES_fullmap', Lk)        

  end select


  write(message(1), '(a)') 'Done'
  call messages_info(1)

  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()
  
  SAFE_DEALLOCATE_A(PESK)    

  contains



    ! ===============================
    ! = get the laser polarization  =
    ! ===============================
    subroutine get_laser_polarizaion(lPol)
       FLOAT,   intent(out) :: lPol(:) 
       
        type(block_t)       :: blk
        integer             :: no_l
        
        PUSH_SUB(get_laser_polarization)
        
        no_l = 0
        if(parse_block(datasets_check('TDExternalFields'), blk) == 0) then
          no_l = parse_block_n(blk)

          call parse_block_float(blk, 0, 1, lPol(1))
          call parse_block_float(blk, 0, 2, lPol(2))
          call parse_block_float(blk, 0, 3, lPol(3))


          call parse_block_end(blk)
        end if
        
        if(no_l > 1) then
          message(1)="There are more the one external field. Polarization will be selected "
          message(2)="from the first field. Use -V to change axis."
          call messages_info(2)
        end if

        POP_SUB(get_laser_polarization)
    end subroutine get_laser_polarizaion




end program photoelectron_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
