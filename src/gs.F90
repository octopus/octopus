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

#include "global.h"

module ground_state
  use global
  use system
  use hamiltonian
  use lcao
  use states
  use restart
  use scf

  implicit none

  private
  public :: ground_state_run, &
            ground_state_init

contains

  ! ---------------------------------------------------------
  integer function ground_state_run(sys, h) result(ierr)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h

    type(scf_type) :: scfv

    call push_sub('ground_state_run')
    ierr = 0

    ! allocate wfs
    allocate(sys%st%X(psi)(sys%m%np, sys%st%d%dim, sys%st%nst, sys%st%d%nik))
    
    ! load wave-functions
    message(1) = 'Info: Loading wave-functions'
    call write_info(1)

    call X(restart_read) ("tmp/restart_gs", sys%st, sys%m, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not load wave-functions: Starting from scratch"
      call write_warning(1)
      ierr = 1
    else
      ! setup Hamiltonian
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)
      call X(system_h_setup) (sys, h)
    
      ! run self consistency
#ifdef COMPLEX_WFNS
      message(1) = 'Info: SCF using complex wavefunctions.'
#else
      message(1) = 'Info: SCF using real wavefunctions.'
#endif
      call write_info(1)

      call scf_init(scfv, sys%m, sys%st, h)
      call scf_run(scfv, sys%m, sys%f_der, sys%st, sys%geo, h, sys%outp)
      call scf_end(scfv)
    end if

    ! clean up
    deallocate(sys%st%X(psi)); nullify(sys%st%X(psi))
    call pop_sub()
  end function ground_state_run


  ! ---------------------------------------------------------
  subroutine ground_state_init(sys, h)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h

    integer :: err
    logical :: lcao_start 
    type(lcao_type) :: lcao_data

    call push_sub('ground_state_init')

    ! allocate wfs
    allocate(sys%st%X(psi)(sys%m%np, sys%st%d%dim, sys%st%nst, sys%st%d%nik))

    message(1) = 'Info: Random generating starting wavefunctions.'
    call write_info(1)
      
    ! wave functions are simply random gaussians
    call states_generate_random(sys%st, sys%m)

    ! this is certainly a better density
    call system_guess_density(sys%m, sys%geo, sys%st%qtot, sys%st%d%nspin, &
       sys%st%d%spin_channels, sys%st%rho)      

    call loct_parse_logical("LCAOStart", .true., lcao_start)
    if(lcao_start) then
      call X(h_calc_vhxc)(h, sys%m, sys%f_der, sys%st, calc_eigenval=.true.)

      message(1) = 'Info: Performing initial LCAO calculation.'
      call write_info(1)
     
      lcao_data%state = 0 ! Uninitialized here.
      call lcao_init(lcao_data, sys%m, sys%f_der, sys%st, sys%geo, h)
      if(lcao_data%state == 1) then
        call lcao_wf(lcao_data, sys%m, sys%st, h)
        call lcao_end(lcao_data)
        
        call states_fermi(sys%st, sys%m)                         ! occupations
        call states_write_eigenvalues(stdout, sys%st%nst, sys%st)
      end if
    end if
   
        ! write wave-functions to disk
    call X(restart_write)("tmp/restart_gs", sys%st, sys%m, err)
    if(err.ne.0) then
      message(1) = 'Unsuccesfull write of "tmp/restart_gs"'
      call write_fatal(1)
    end if

    ! clean up
    deallocate(sys%st%X(psi))
    nullify(sys%st%X(psi))

    call pop_sub()
  end subroutine ground_state_init

end module ground_state
