#include "config.h"

module scf
use system
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
  scf%conv_abs_ener = fdf_double("ConvAbsEnergy", 0.0_r8)
  scf%conv_rel_ener = fdf_double("ConvRelEnergy", 0.0_r8)

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

subroutine scf_run(sys, h)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  
  type(scf_type) :: scf

  integer :: iter, iunit, ik, ist, id
  real(r8) :: old_etot, diff(sys%st%nst, sys%st%nik)
  logical :: finish

  sub_name = 'scf_run'; call push_sub()
  call scf_init(scf, sys)

  do iter = 1, scf%max_iter
    if(clean_stop()) exit

    call eigen_solver_run(scf%eigens, sys, h, iter, diff)

    ! compute eigenvalues
    call dhamiltonian_eigenval(h, sys, 1, sys%st%nst) ! eigenvalues
    call states_fermi(sys%st)                         ! occupations

    ! output eigenvalues
    call scf_write_eigenvalues(stdout, sys%st%nst, sys%st, diff)

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
    call io_assign(iunit)
    open(iunit, status='unknown', file=trim(sys%sysname)//".restart", form='unformatted')
    
    write(iunit) sys%m%box_shape, sys%m%h, sys%m%rsize, sys%m%zsize
    write(iunit) sys%m%np, sys%st%dim, 1, sys%st%nst, sys%st%nik, sys%st%ispin
    do ik = 1, sys%st%nik
      do ist = 1, sys%st%nst
        do id = 1, sys%st%dim
          write(iunit) sys%st%R_FUNC(psi)(1:sys%m%np, id, ist, ik)
        end do
      end do
    end do
    ! eigenvalues are also needed ;)
    do ik = 1, sys%st%nik
      write(iunit) sys%st%eigenval(:, ik)
    end do
    call io_close(iunit)

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
!  call  systm_write_scf(sys, scf, iter, diff, finish, ion_ion_forces, local_forces, & 
!                        non_local_forces)

  call scf_end(scf)
  call pop_sub()
end subroutine scf_run

subroutine scf_write_eigenvalues(iunit, nst, st, error)
  integer, intent(in) :: iunit, nst
  type(states_type), intent(IN) :: st
  real(r8), intent(in) :: error(nst, st%nik)

  integer ik, j
  real(r8) :: o(st%nik)

  if(iunit==stdout.and.conf%verbose<=20) return

  message(1) = 'Info: Eigenvalues ['//trim(units_out%mass%abbrev)//']'
  call write_info(1, iunit)

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    write(iunit, '(a4)', advance='no') '#'
    do ik = 1, st%nik
      write(iunit, '(1x,a12,1x,a12,2x,a10,i1,a1)', advance='no') &
           'Eigenvalue ', 'Occupation ', 'Error (', ik, ')'
    end do
    write(iunit, '(1x)', advance='yes')
    
    do j = 1, nst
      if(j > st%nst) then
        o = 0._r8
      else
        o = st%occ(j, :)
      end if
      
      write(iunit, '(i4)', advance='no') j
      do ik = 1, st%nik
        write(iunit, '(1x,f12.6,1x,f12.6,a2,f12.8,a1)', advance='no') &
             st%eigenval(j, ik), o(ik), ' (', error(j, ik), ')'
      end do
      write(iunit, '(1x)', advance='yes')
    end do

#ifdef HAVE_MPI
  end if
#endif
  
end subroutine scf_write_eigenvalues


end module scf
