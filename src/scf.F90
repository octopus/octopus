#include "config.h"

module scf
use units
use system
use states
use hamiltonian
use eigen_solver
use mix

implicit none

type scf_type  ! some variables used for the scf cycle
  integer :: max_iter ! maximum number of scf iterations
  real(r8) :: conv_abs_dens, conv_rel_dens, &
       conv_abs_ener, conv_rel_ener ! several convergence criteria
  
  real(r8) :: abs_dens, rel_dens, abs_ener, rel_ener
  
  type(mix_type) :: smix
  type(eigen_solver_type) :: eigens
end type scf_type
  
contains

subroutine scf_init(scf, sys)
  type(scf_type), intent(inout) :: scf
  type(system_type), intent(IN) :: sys

  sub_name = 'systm_scf_init'; call push_sub()

  scf%max_iter = fdf_integer("MaximumIter", 200)
  scf%conv_abs_dens = fdf_double("ConvAbsDens", 1e-5_r8)
  scf%conv_rel_dens = fdf_double("ConvRelDens", 0.0_r8)
  scf%conv_abs_ener = fdf_double("ConvAbsEnergy", 0.0_r8) * units_inp%energy%factor
  scf%conv_rel_ener = fdf_double("ConvRelEnergy", 0.0_r8) * units_inp%energy%factor

  if(scf%max_iter <= 0 .and. &
      scf%conv_abs_dens <= 0.0_r8 .and. scf%conv_rel_dens <= 0.0_r8 .and. &
      scf%conv_abs_ener <= 0.0_r8 .and. scf%conv_rel_ener <= 0.0_r8) then
    message(1) = "Input: Not all convergence criteria can be <= 0"
    message(2) = "Please set one of the following:"
    message(3) = "MaximumIter | ConvAbsDens | ConvRelDens | ConvAbsEnergy | ConvRelEnergy"
    call write_fatal(3)
  end if

  if(scf%max_iter <= 0) scf%max_iter = huge(scf%max_iter)

  ! Handle mixing now...
  call mix_init(scf%smix, sys%m, sys%st)

  ! now the eigen solver stuff
  call eigen_solver_init(scf%eigens)

  call pop_sub()
  return
end subroutine scf_init

subroutine scf_end(scf)
  type(scf_type), intent(inout) :: scf

  call mix_end(scf%smix)

  return
end subroutine scf_end

subroutine scf_run(scf, sys, h)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  type(scf_type), intent(inout) :: scf

  integer :: iter, iunit, ik, ist, id
  real(r8) :: old_etot, diff(sys%st%nst, sys%st%nik)
  logical :: finish

  sub_name = 'scf_run'; call push_sub()

  do iter = 1, scf%max_iter
    if(clean_stop()) exit

    call eigen_solver_run(scf%eigens, sys, h, iter, diff)

    ! compute eigenvalues
    call R_FUNC(hamiltonian_eigenval) (h, sys, 1, sys%st%nst) ! eigenvalues
    call states_fermi(sys%st)                         ! occupations

    ! output eigenvalues
    call states_write_eigenvalues(stdout, sys%st%nst, sys%st, diff)

    ! now compute total energy
    old_etot = h%etot
    call hamiltonian_energy(h, sys, -1)
    scf%abs_ener = abs(old_etot - h%etot)
    scf%rel_ener = scf%abs_ener / abs(h%etot)

    ! compute new density
    call mix_dens(scf%smix, iter, sys%st, sys%m, scf%abs_dens)
    scf%rel_dens = scf%abs_dens / sys%st%qtot

    ! compute new potentials
    call hamiltonian_setup(h, sys)

    ! save restart information
    call R_FUNC(states_write_restart)(trim(sys%sysname)//".restart", sys%m, sys%st)

    ! are we finished?
    finish = &
        (scf%conv_abs_dens > 0.0_r8 .and. scf%abs_dens <= scf%conv_abs_dens) .or. &
        (scf%conv_rel_dens > 0.0_r8 .and. scf%rel_dens <= scf%conv_rel_dens) .or. &
        (scf%conv_abs_ener > 0.0_r8 .and. scf%abs_ener <= scf%conv_abs_ener) .or. &
        (scf%conv_rel_ener > 0.0_r8 .and. scf%rel_ener <= scf%conv_rel_ener)

    write(message(1), '(a,i4,a,e14.8,a,e14.8)') &
         'Info: iter = ', iter, ' abs_dens = ', scf%abs_dens, &
         ' abs_ener = ', scf%abs_ener
    call write_info(1)

    if(finish) then
      write(message(1), '(a, i4, a)')'Info: SCF converged in ', iter, ' iterations'
      call write_info(1)

!      call R_FUNC(ion_forces) (sys, rpsi, forces, &
!                               ion_ion_forces=ion_ion_forces, local_forces=local_forces, &
!                               non_local_forces=non_local_forces)
      exit
    end if
  end do

  if(.not.finish) then
    message(1) = 'SCF *not* converged!'
    call write_warning(1)
  end if

  ! output final information
  call  scf_write_static()

  call pop_sub()

contains

subroutine scf_write_static()
  integer iunit, i

  call io_assign(iunit)
  open(iunit, status='unknown', file=trim(sys%sysname)//'.static')

  ! mesh
  write(iunit, '(a,a)') 'System name: ', sys%sysname
  write(iunit, '(1x)')

  write(iunit, '(a)') 'Mesh:'
  call mesh_write_info(sys%m, iunit)
  write(iunit,'(1x)')

  !write(iunit, '(a)') 'Geometry [A]:'
  !call write_geom_info(sys, iunit)
  !write(iunit,'(1x)')

  if(.not. h%ip_app) then
    write(iunit, '(a)') 'Exchange and correlation functionals:'
    call xc_write_info(h%xc, iunit)
  else
    write(iunit, '(a)') 'Independent Particles'
  end if
  write(iunit,'(1x)')

  ! scf information
  if(finish) then
    write(iunit, '(a, i4, a)')'SCF converged in ', iter, ' iterations'
  else
    write(iunit, '(a)') 'SCF *not* converged!'
  end if
  write(iunit,'(1x)')

  call states_write_eigenvalues(iunit, sys%st%nst, sys%st, diff)
  write(iunit,'(1x)')

  write(iunit, '(a)') 'Energy [eV]:'
  call hamiltonian_energy(h, sys, iunit)
  write(iunit,'(1x)')

  write(iunit, '(a)') 'Convergence:'
  write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'abs_dens = ', scf%abs_dens, &
      ' (', scf%conv_abs_dens, ')'
  write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'rel_dens = ', scf%rel_dens, &
      ' (', scf%conv_rel_dens, ')'
  write(iunit, '(6x, a, es14.8,a,es14.8,4a)') 'abs_ener = ', scf%abs_ener, &
      ' (', scf%conv_abs_ener / units_out%energy%factor, ')', &
      ' [',  trim(units_out%energy%abbrev), ']'
  write(iunit, '(6x, a, es14.8,a,es14.8,4a)') 'rel_ener = ', scf%rel_ener, &
      ' (', scf%conv_rel_ener / units_out%energy%factor, ')', &
      ' [',  trim(units_out%energy%abbrev), ']'
  write(iunit,'(1x)') 

! TODO
!  write(iunit,'(a)') 'Forces on the ions [eV/A]'
!  write(iunit,'(a,10x,14x,a,14x,a,14x,a)') ' Ion','x','y','z'
!  do i=1,sys%nions
!     write(iunit,'(i4,a,3f15.6)') i,'  ion-ion:',ion_ion_forces(1:3,i)
!     write(iunit,'(4x,a,3f15.6)')   ' nonlocal:',non_local_forces(1:3,i)
!     write(iunit,'(4x,a,3f15.6)')   '    local:',local_forces(1:3,i)
!     write(iunit,'(4x,a,3f15.6)')   '    total:',ion_ion_forces(1:3,i) + &
!                                                 non_local_forces(1:3,i) + &
!                                                 local_forces(1:3,i)
!  enddo

  call io_close(iunit)
  return
end subroutine scf_write_static

end subroutine scf_run

! TODO use netcdf
subroutine scf_write_restart(filename, m, st)
  character(len=*), intent(in) :: filename
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st

  integer :: iunit, ik, ist, id

  call io_assign(iunit)
  open(iunit, status='unknown', file=trim(filename), form='unformatted')
    
  write(iunit) m%box_shape, m%h, m%rsize, m%zsize
  write(iunit) m%np, st%dim, 1, st%nst, st%nik, st%ispin
  do ik = 1, st%nik
    do ist = 1, st%nst
      do id = 1, st%dim
        write(iunit) st%R_FUNC(psi)(1:m%np, id, ist, ik)
      end do
    end do
  end do
  ! eigenvalues are also needed ;)
  do ik = 1, st%nik
    write(iunit) st%eigenval(:, ik)
  end do

  call io_close(iunit)

end subroutine scf_write_restart

end module scf
