#include "config_F90.h"

program hs_from_acc
  use global
  use spectrum
  use liboct

  integer :: ierr
  character(len=100) :: txt
  type(spec_type) :: s
  type(spec_sh) :: sh

  ! init liboct
  ierr = oct_parse_init(C_string('inp'), C_string('out.oct'))
  if(ierr .ne. 0) then
    message(1) = "Error initializing liboct"
    call write_fatal(1)
  end if
  conf%verbose = 30
  call units_init()

  call oct_parse_double(C_string("SpecStartTime"), 0._r8, s%start_time)
  call oct_parse_double(C_string("SpecEndTime"), -1._r8, s%end_time)
  call oct_parse_double(C_string("SpecEnergyStep"), 0.05_r8, s%energy_step)
  call oct_parse_double(C_string("SpecMaxEnergy"), 20._r8, s%max_energy)

  ! adjust units
  s%start_time      = s%start_time      * units_inp%time%factor
  s%end_time        = s%end_time        * units_inp%time%factor
  s%energy_step     = s%energy_step     * units_inp%energy%factor
  s%max_energy      = s%max_energy      * units_inp%energy%factor

  call oct_parse_str('HSPolarization', 'z', txt)
  sh%pol = txt(1:1)
  if(sh%pol.ne.'x' .and. sh%pol.ne.'y' .and. sh%pol.ne.'z' .and. &
       sh%pol.ne.'+' .and. sh%pol.ne.'-') then
    message(1) = "HSPolarization has an invalid value"
    message(2) = "Valid values are ('x' | 'y' | 'z' | '+' | '-')"
    call write_fatal(2)
  end if

  call oct_parse_str('SystemName', 'system', txt)
  call spectrum_hs_from_acc(trim(txt), trim(txt)//'.acc-hs', s, sh, .true.)

  deallocate(sh%sp)
  stop  
end program hs_from_acc
