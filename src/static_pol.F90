#include "config.h"

module static_pol
use liboct
use io
use units
use system
use hamiltonian
use scf

implicit none

contains

subroutine static_pol_run(scf, sys, h)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  type(scf_type), intent(inout) :: scf
  
  integer :: iunit, ios, i_start, i, j, is
  real(r8) :: e_field
  real(r8), allocatable :: Vpsl_save(:)
  real(r8) :: dipole(0:3, 3, sys%st%nspin)
  logical :: out_pol

  sub_name = 'static_pol_run'; call push_sub()

  ! read in e_field value
  call oct_parse_double(C_string('POLStaticField'), &
       0.001_r8/units_inp%energy%factor*units_inp%length%factor, e_field)
  e_field = e_field * units_inp%energy%factor / units_inp%length%factor
  if (e_field <= 0._r8) then
    write(message(1), '(a,e14.6,a)') "Input: '", e_field, "' is not a valid POLStaticField"
    message(2) = '(0 < POLStaticField)'
    call write_fatal(2)
  end if

  ! open restart file
  call io_assign(iunit)
  open(iunit, file=trim(sys%sysname)//'.pol_restart', status='unknown')
  do i_start = 0, 3
    read(iunit, fmt=*, iostat=ios) dipole(i_start, :, :)
    if(ios.ne.0) exit
  end do
  if(ios.eq.0) then ! all information is already calculated
    i_start = 4 ! do not do the do loop
    out_pol = .true.
  else
    out_pol = .false.
  end if
  call io_close(iunit)

  ! save local pseudopot
  allocate(Vpsl_save(sys%m%np))
  Vpsl_save = h%Vpsl

  do i = i_start, 3
    if(clean_stop()) exit

    write(message(1), '(a, i2)')'Info: Calculating dipole moment for field ', i
    call write_info(1)

    h%Vpsl = Vpsl_save

    select case(i)
    case(1)
      h%Vpsl = h%Vpsl - sys%m%Lx(:)*sys%m%h*e_field
    case(2)
      h%Vpsl = h%Vpsl - sys%m%Ly(:)*sys%m%h*e_field
    case(3)
      h%Vpsl = h%Vpsl - sys%m%Lz(:)*sys%m%h*e_field
    end select

    call scf_run(scf, sys, h)
    
    ! calculate dipole
    do is = 1, sys%st%nspin
      do j = 1, 3
        select case (j)
        case(1) 
          dipole(i, j, is) = sum(sys%st%rho(:, is)*sys%m%Lx(:))*sys%m%H*sys%m%vol_pp
        case(2)
          dipole(i, j, is) = sum(sys%st%rho(:, is)*sys%m%Ly(:))*sys%m%H*sys%m%vol_pp
        case(3)
          dipole(i, j, is) = sum(sys%st%rho(:, is)*sys%m%Lz(:))*sys%m%H*sys%m%vol_pp
        end select
      end do
    end do

    ! output dipole
    call io_assign(iunit)
    open(iunit, file=trim(sys%sysname)//'.pol_restart', status='unknown')
    do j = 0, i
      write(iunit, fmt=*) dipole(j, :, :)
    end do
    call io_close(iunit)

    if(i == 3) then 
      out_pol = .true.
    end if
  end do

  deallocate(Vpsl_save)

  if(out_pol) then ! output pol file
    call io_assign(iunit)
    open(iunit, file=trim(sys%sysname)//'.pol', status='unknown')
    write(iunit, '(3a)') '# Static polarizability tensor [', &
         trim(units_out%length%abbrev), '^3]'
         
    do is = 1, sys%st%nspin
      do j = 1, 3
        write(iunit, '(3f12.6)') (dipole(j, :, is) - dipole(0, :, is))/e_field &
             / units_out%length%factor**3
      end do
    end do
    call io_close(iunit)
  end if

  call pop_sub()
  return
end subroutine static_pol_run

end module static_pol
