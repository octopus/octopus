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

module run_prog
use global
use oct_parser
use io
use units
use states
use system
use hamiltonian
use lcao
use scf
use unocc
use timedep
use static_pol
use geom_opt
use phonons
use opt_control
use pulpo

implicit none

type(system_type) :: sys
type(hamiltonian_type) :: h
type(scf_type) :: scfv
type(unocc_type) :: unoccv
integer :: calc_mode

! run stack
integer, private :: i_stack(100), instr
integer, private, parameter ::   &
     M_START_STATIC_CALC   = 1,  &
     M_RESUME_STATIC_CALC  = 2,  &
     M_START_UNOCC_STATES  = 3,  &
     M_RESUME_UNOCC_STATES = 4,  &
     M_START_TD            = 5,  &
     M_RESUME_TD           = 6,  &
     M_START_STATIC_POL    = 7,  &
     M_RESUME_STATIC_POL   = 8,  &
     M_BO_MD               = 9,  &
     M_GEOM_OPT            = 10, &
     M_PHONONS             = 11, &
     M_OPT_CONTROL         = 12, &
     M_PULPO_A_FEIRA       = 99

integer, private, parameter :: &
     I_SETUP_RPSI          =  1,  &
     I_END_RPSI            =  2,  &
     I_RANDOMIZE_RPSI      =  3,  &
     I_LOAD_RPSI           =  4,  &
     I_SETUP_HAMILTONIAN   =  5,  &
     I_SCF                 =  6,  &
     I_LCAO                =  7,  &
     I_SETUP_UNOCC         = 11,  &
     I_END_UNOCC           = 12,  &
     I_RANDOMIZE_UNOCC     = 13,  &
     I_LOAD_UNOCC          = 14,  &
     I_UNOCC_RUN           = 15,  &
     I_SETUP_TD            = 21,  &
     I_END_TD              = 22,  &
     I_INIT_ZPSI           = 23,  &
     I_LOAD_ZPSI           = 24,  &
     I_TD                  = 25,  &
     I_SETUP_OCC_AN        = 26,  &
     I_END_OCC_AN          = 27,  &
     I_START_POL           = 30,  &
     I_POL_SCF             = 31,  &
     I_GEOM_OPT            = 40,  &
     I_PHONONS             = 50,  &
     I_OPT_CONTROL         = 60,  &
     I_PULPO               = 99

contains

subroutine run()
  type(td_type), pointer :: td
  integer :: iunit, i
  real(r8) :: x
  logical :: log
  character(len=100) :: filename

  call push_sub('run')

  call run_init()

  instr = 0
  call define_run_modes()

  program: do
    if(instr <= 0) exit program

    select case(i_stack(instr))
    case(I_SETUP_RPSI)
      message(1) = 'Info: Allocating rpsi.'
      call write_info(1)
        
      allocate(sys%st%R_FUNC(psi)(sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))

    case(I_END_RPSI)
      message(1) = 'Info: Deallocating rpsi.'
      call write_info(1)

      if(associated(sys%st%R_FUNC(psi))) then
        deallocate(sys%st%R_FUNC(psi)); nullify(sys%st%R_FUNC(psi))
      end if

    case(I_RANDOMIZE_RPSI)
      message(1) = 'Info: Random generating starting wavefunctions.'
      call write_info(1)

      ! wave functions are simply random gaussians
      call states_generate_random(sys%st, sys%m)
      call R_FUNC(calcdens)(sys%st, sys%m%np, sys%st%rho)
      ! this is certainly a better density
      call lcao_dens(sys, sys%st%nspin, sys%st%rho)

    case(I_LOAD_RPSI)
      message(1) = 'Info: Loading rpsi.'
      call write_info(1)
      
      if(R_FUNC(states_load_restart)("tmp/restart.static", &
           sys%m, sys%st)) then
        call R_FUNC(calcdens)(sys%st, sys%m%np, sys%st%rho)
      else
        ! run scf unless it is already in the stack
        if(calc_mode .ne. M_RESUME_STATIC_CALC .and. &
             calc_mode .ne. M_START_STATIC_POL .and. calc_mode .ne. M_RESUME_STATIC_POL) then
          i_stack(instr) = I_SCF;               instr = instr + 1
          i_stack(instr) = I_LCAO;              instr = instr + 1
          i_stack(instr) = I_SETUP_HAMILTONIAN; instr = instr + 1
        end if
        i_stack(instr) = I_RANDOMIZE_RPSI
        cycle program         
      end if

    case(I_SETUP_HAMILTONIAN)
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)

      call R_FUNC(hamiltonian_setup)(h, sys%m, sys%st, sys) ! get potentials
      call states_fermi(sys%st, sys%m)                      ! occupations
      call hamiltonian_energy(h, sys%st, sys%eii, -1)       ! total energy

    case(I_SCF)
#ifdef COMPLEX_WFNS
      message(1) = 'Info: SCF using complex wavefunctions.'
#else
      message(1) = 'Info: SCF using real wavefunctions.'
#endif
      call write_info(1)
      call scf_init(scfv, sys)
      call scf_run(scfv, sys, h)
      call scf_end(scfv)

    case(I_LCAO)
      call oct_parse_logical("LCAOStart", .true., log)
      do i = 1, sys%nspecies
          log = log .and. (.not.sys%specie(i)%local)
      enddo
      if(log) then
        message(1) = 'Info: Performing LCAO calculation.'
        call write_info(1)
        call lcao_init(sys, h)
        call lcao_wf(sys, h)
        call lcao_end
        call states_fermi(sys%st, sys%m)                         ! occupations
        call states_write_eigenvalues(stdout, sys%st%nst, sys%st)
      end if

    case(I_SETUP_UNOCC)
      message(1) = 'Info: Initializing unoccupied states.'
      call write_info(1)

      call unocc_init(unoccv, sys%m, sys%st, sys)
      
    case(I_END_UNOCC)
      message(1) = 'Info: Finalizing unoccupied states.'
      call write_info(1)

      call unocc_end(unoccv)
      
    case(I_RANDOMIZE_UNOCC)
      message(1) = 'Info: Generating starting unoccupied states.'
      call write_info(1)
      
      ! first copy the occupied states
      unoccv%st%R_FUNC(psi)(:,:,1:sys%st%nst,:) = sys%st%R_FUNC(psi)(:,:,1:sys%st%nst,:)
      unoccv%st%occ(1:sys%st%nst,:) = sys%st%occ(1:sys%st%nst,:)
      call states_end(sys%st) ! to save memory

      call states_generate_random(unoccv%st, sys%m, sys%st%nst+1)

    case(I_LOAD_UNOCC)
      message(1) = 'Info: Loading unoccupied states.'
      call write_info(1)

      if(.not.R_FUNC(states_load_restart)("tmp/restart.occ", &
           sys%m, unoccv%st)) then
        if(calc_mode .ne. M_RESUME_UNOCC_STATES) then
          i_stack(instr) = I_UNOCC_RUN
          instr = instr + 1
        end if
        i_stack(instr) = I_RANDOMIZE_UNOCC
        cycle program 
      end if
      
    case(I_UNOCC_RUN)
      message(1) = 'Info: Calculation of unoccupied states.'
      call write_info(1)

      call unocc_run(unoccv, sys, h)

    case(I_SETUP_TD)
      message(1) = 'Info: Setting up TD.'
      call write_info(1)

      ! initialize structure
      allocate(td)
      ! this allocates zpsi
      call td_init(td, sys, sys%m, sys%st, h)

    case(I_END_TD)
      message(1) = 'Info: Cleaning up TD.'
      call write_info(1)

      call td_end(td)
      deallocate(td); nullify(td)
      deallocate(sys%st%zpsi); nullify(sys%st%zpsi)
      
    case(I_INIT_ZPSI)
      message(1) = 'Info: Initializing zpsi.'
      call write_info(1)

      ! load zpsi from static file.
      if(zstates_load_restart("tmp/restart.static", sys%m, sys%st)) then
        call zcalcdens(sys%st, sys%m%np, sys%st%rho, reduce=.true.)
        call zhamiltonian_setup(h, sys%m, sys%st, sys)
        x = minval(sys%st%eigenval(sys%st%st_start, :))
#ifdef HAVE_MPI
        call MPI_BCAST(x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, i)
#endif
        call hamiltonian_span(h, minval(sys%m%h), x)
        call hamiltonian_energy(h, sys%st, sys%eii, -1, reduce=.true.)        
      else
        i_stack(instr) = I_INIT_ZPSI
        instr = instr + 1; i_stack(instr) = I_SETUP_TD
        instr = instr + 1; i_stack(instr) = I_END_RPSI
        instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
        instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
        instr = instr + 1; i_stack(instr) = I_END_TD
        cycle program
      end if

    case(I_LOAD_ZPSI)
      message(1) = 'Info: Loading zpsi.'
      call write_info(1)

      write(filename, '(a,i3.3)') "tmp/restart.td.", mpiv%node
      if(zstates_load_restart(trim(filename), &
           sys%m, sys%st, iter=td%iter, v1=td%v_old(:, :, 1), v2=td%v_old(:, :, 2))) then
        ! define density and hamiltonian
        call zcalcdens(sys%st, sys%m%np, sys%st%rho, reduce=.true.)
        call zhamiltonian_setup(h, sys%m, sys%st, sys)
        call hamiltonian_span(h, minval(sys%m%h), minval(sys%st%eigenval(1,:)))
        call hamiltonian_energy(h, sys%st, sys%eii, -1, reduce=.true.)        
        
      else
        i_stack(instr) = I_INIT_ZPSI
        cycle program
      end if

    case(I_TD)
      message(1) = 'Info: Time-dependent simulation.'
      call write_info(1)

      call td_run(td, unoccv%st, sys, h)

    case(I_SETUP_OCC_AN)

      if(td%out_proj) then
        message(1) = 'Info: Seting up occupational analysis.'
        call write_info(1)
        instr = instr - 1
        call m_load_unocc()
        cycle program
      else
        allocate(unoccv%st) ! Otherwise td_run crashes
      end if

    case(I_END_OCC_AN)
      message(1) = 'Info: Cleaning up occupational analysis.'
      call write_info(1)

      if(td%out_proj) then
        i_stack(instr) = I_END_UNOCC
        cycle program
      else
        deallocate(unoccv%st)
      end if

    case(I_START_POL)
      ! just delete the pol file
      message(1) = 'Info: Starting static polarizability calculation'
      call write_info(1)

      call io_assign(iunit)
      open(iunit, file='tmp/restart.pol', status='unknown')
      write(iunit, '(a)') ' '
      call io_close(iunit)
        
    case(I_POL_SCF)
      message(1) = 'Info: Calculating static polarizability'
      call write_info(1)

      call scf_init(scfv, sys)
      call static_pol_run(scfv, sys, h)
      call scf_end(scfv)

    case(I_GEOM_OPT)
      message(1) = 'Info: Performing geometry optimization'
      call write_info(1)

      call scf_init(scfv, sys)
      call geom_opt_run(scfv, sys, h)
      call scf_end(scfv)

    case(I_PHONONS)
      message(1) = 'Info: Calculating phonon frequencies'
      call write_info(1)

      call scf_init(scfv, sys)
      call phonons_run(scfv, sys, h)
      call scf_end(scfv)

    case(I_OPT_CONTROL)
      message(1) = 'Info: Optimum control.'
      call write_info(1)

      call opt_control_run(td, sys, h)

    case(I_PULPO)
      call pulpo_print()

    case default
      write(message(1), '(a,i3,a)') "Instruction ", i_stack(instr), " not defined!"
      call write_warning(1)
    end select

    instr = instr - 1
  end do program
      
  call run_end()

  call pop_sub()
end subroutine run

subroutine run_init()
  ! initialize some stuff

  call oct_parse_int('CalculationMode', 1, calc_mode)
  if( (calc_mode < 1 .or. calc_mode > 12) .and. (calc_mode .ne. M_PULPO_A_FEIRA)) then
    write(message(1), '(a,i2,a)') "Input: '", calc_mode, "' is not a valid CalculationMode"
    message(2) = '  Calculation Mode = 1 <= start static calculation'
    message(3) = '                   = 2 <= resume static calculation'
    message(4) = '                   = 3 <= calculate unocuppied states'
    message(5) = '                   = 4 <= resume unocuppied states calculation'
    message(6) = '                   = 5 <= start td'
    message(7) = '                   = 6 <= resume td'
    message(8) = '                   = 7 <= start static polarizability'
    message(9) = '                   = 8 <= resume static polarizability'
    message(10)= '                   = 9 <= perform Born-Oppenheimer MD'
    message(11)= '                   =10 <= perform geometry minimization'
    message(12)= '                   =11 <= calculate phonon frequencies'
    message(13)= '                   =12 <= optimum control'
    message(14)= '                   =99 <= prints out the "Pulpo a Feira" recipe'
    call write_fatal(13)
  end if

  ! print dimension info
  write(message(1), '(a,i1,a)') 'Info: Octopus will run in ', conf%dim, ' dimension(s)'
  call write_info(1)
  write(message(1), '(a,i1,a)') 'Info: Octopus will treat system as periodic in ', conf%periodic_dim, ' dimension(s)'
  call write_info(1)

  message(1) = "Info: Boundary conditions:"
  if(conf%boundary_zero_derivative) then
    message(1) = trim(message(1)) + " zero derivatives"
  else
    message(1) = trim(message(1)) + " zero wave-functions"
  end if
  call write_info(1)


  if(calc_mode .ne. M_PULPO_A_FEIRA) then
    call units_init()
    call system_init(sys)
    call hamiltonian_init(h, sys)
    if(h%classic_pot > 0) call generate_classic_pot(h, sys)
    call generate_external_pot(h, sys)
  endif
  
end subroutine run_init

subroutine run_end()
  if(calc_mode .ne. M_PULPO_A_FEIRA) then
    call system_end(sys)
    call hamiltonian_end(h)
  endif
end subroutine run_end

subroutine m_load_unocc()
  instr = instr + 1; i_stack(instr) = I_END_RPSI
  instr = instr + 1; i_stack(instr) = I_LOAD_UNOCC   
  instr = instr + 1; i_stack(instr) = I_SETUP_UNOCC
  instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
  instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
  instr = instr + 1; i_stack(instr) = I_SETUP_RPSI

  return
end subroutine m_load_unocc

subroutine define_run_modes()
  select case(calc_mode)
  case(M_START_STATIC_CALC)
    instr = instr + 1; i_stack(instr) = I_END_RPSI
    instr = instr + 1; i_stack(instr) = I_SCF
    instr = instr + 1; i_stack(instr) = I_LCAO
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
    instr = instr + 1; i_stack(instr) = I_RANDOMIZE_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
  case(M_RESUME_STATIC_CALC)
    instr = instr + 1; i_stack(instr) = I_END_RPSI
    instr = instr + 1; i_stack(instr) = I_SCF
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
    instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
  case(M_START_UNOCC_STATES)
    instr = instr + 1; i_stack(instr) = I_END_UNOCC
    instr = instr + 1; i_stack(instr) = I_UNOCC_RUN
    instr = instr + 1; i_stack(instr) = I_END_RPSI
    instr = instr + 1; i_stack(instr) = I_RANDOMIZE_UNOCC
    instr = instr + 1; i_stack(instr) = I_SETUP_UNOCC
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN
    instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
  case(M_RESUME_UNOCC_STATES)
    instr = instr + 1; i_stack(instr) = I_END_UNOCC
    instr = instr + 1; i_stack(instr) = I_UNOCC_RUN
    call m_load_unocc()
  case(M_START_TD)
    instr = instr + 1; i_stack(instr) = I_END_TD
    instr = instr + 1; i_stack(instr) = I_END_OCC_AN
    instr = instr + 1; i_stack(instr) = I_TD
    instr = instr + 1; i_stack(instr) = I_SETUP_OCC_AN
    instr = instr + 1; i_stack(instr) = I_INIT_ZPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_TD
  case(M_RESUME_TD)
    instr = instr + 1; i_stack(instr) = I_END_TD
    instr = instr + 1; i_stack(instr) = I_END_OCC_AN
    instr = instr + 1; i_stack(instr) = I_TD
    instr = instr + 1; i_stack(instr) = I_SETUP_OCC_AN
    instr = instr + 1; i_stack(instr) = I_LOAD_ZPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_TD    
  case(M_START_STATIC_POL)
    instr = instr + 1; i_stack(instr) = I_END_RPSI
    instr = instr + 1; i_stack(instr) = I_POL_SCF
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
    instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
    instr = instr + 1; i_stack(instr) = I_START_POL
  case(M_RESUME_STATIC_POL)
    instr = instr + 1; i_stack(instr) = I_END_RPSI
    instr = instr + 1; i_stack(instr) = I_POL_SCF
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
    instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
  case(M_GEOM_OPT)
    instr = instr + 1; i_stack(instr) = I_END_RPSI
    instr = instr + 1; i_stack(instr) = I_GEOM_OPT
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
    instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
  case(M_PHONONS)
    instr = instr + 1; i_stack(instr) = I_END_RPSI
    instr = instr + 1; i_stack(instr) = I_PHONONS
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
    instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
  case(M_OPT_CONTROL)
    instr = instr + 1; i_stack(instr) = I_END_TD
    instr = instr + 1; i_stack(instr) = I_OPT_CONTROL
    instr = instr + 1; i_stack(instr) = I_SETUP_TD
  case(M_PULPO_A_FEIRA)
    instr = instr + 1; i_stack(instr) = I_PULPO
  end select

end subroutine define_run_modes

end module run_prog
