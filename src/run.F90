#include "config.h"

module run_prog
use global
use liboct
use io
use units
use states
use system
use hamiltonian
use lcao
use scf
use unocc
use static_pol
use timedep
use pulpo

implicit none

type(system_type) :: sys
type(hamiltonian_type) :: h
type(scf_type) :: scfv
type(unocc_type) :: unoccv
integer :: calc_mode

! run stack
integer, private :: i_stack(100), instr
integer, private, parameter :: &
     M_START_STATIC_CALC   = 1, &
     M_RESUME_STATIC_CALC  = 2, &
     M_START_UNOCC_STATES  = 3, &
     M_RESUME_UNOCC_STATES = 4, &
     M_START_TD            = 5, &
     M_RESUME_TD           = 6, &
     M_START_STATIC_POL    = 7, &
     M_RESUME_STATIC_POL   = 8, &
     M_BO_MD               = 9, &
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

     I_PULPO               = 99

contains

subroutine run()
  type(td_type), pointer :: td
  R_TYPE, pointer :: aux_psi(:,:,:,:)
  real(r8), pointer :: aux_eigen(:,:)
  integer :: i, aux_i1, aux_i2, iunit
  logical :: log

  sub_name = 'run'; call push_sub()

  call run_init()

  instr = 0
  call define_run_modes()

  program: do
    if(instr <= 0) exit program

    select case(i_stack(instr))
    case(I_SETUP_RPSI)
      message(1) = 'Info: Allocating rpsi.'
      call write_info(1)
        
      allocate(sys%st%R_FUNC(psi)(0:sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))

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
      ! we will need some starting density in order to define the hamiltonian
      call lcao_dens(sys, sys%st%nspin, sys%st%rho)

      ! TODO: non-collinear spin
      ! the off-diagonal densities are set to zero
      !if(sys%st%ispin > 2) then
      !  sys%st%rho(:,i+1:sys%st%ispin) = 0._r8
      !end if

    case(I_LOAD_RPSI)
      message(1) = 'Info: Loading rpsi.'
      call write_info(1)
      
      if(R_FUNC(states_load_restart)(trim(sys%sysname)//".restart", &
           sys%m, sys%st)) then
        call R_FUNC(calcdens)(sys%st, sys%m%np, sys%st%rho)
      else
        ! run scf unless it is already in the stack
        if(calc_mode .ne. M_RESUME_STATIC_CALC .and. &
             calc_mode .ne. M_START_STATIC_POL .and. calc_mode .ne. M_RESUME_STATIC_POL) then
          i_stack(instr) = I_SCF;               instr = instr + 1
          i_stack(instr) = I_SETUP_HAMILTONIAN; instr = instr + 1
        end if
        i_stack(instr) = I_RANDOMIZE_RPSI
        cycle program         
      end if

    case(I_SETUP_HAMILTONIAN)
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)

      call hamiltonian_setup(h, sys)                    ! get potentials
      call R_FUNC(hamiltonian_eigenval)(h, sys, 1, sys%st%nst) ! eigenvalues
      call states_fermi(sys%st)                         ! occupations
      call hamiltonian_energy(h, sys, -1)               ! get the total energy

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
      if(log) then
        message(1) = 'Info: Performing LCAO calculation.'
        call write_info(1)
        call lcao_wf(sys, h)
      end if

    case(I_SETUP_UNOCC)
      message(1) = 'Info: Initializing unoccupied states.'
      call write_info(1)

      call unocc_init(unoccv, sys%m, sys%st)
      
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

      call states_generate_random(unoccv%st, sys%m, sys%st%nst+1)

    case(I_LOAD_UNOCC)
      message(1) = 'Info: Loading unoccupied states.'
      call write_info(1)

      if(.not.R_FUNC(states_load_restart)(trim(sys%sysname)//".occ_restart", &
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

      
      call scf_init(scfv, sys)
      call unocc_run(unoccv, scfv, sys, h)
      call scf_end(scfv)

    case(I_SETUP_TD)
      message(1) = 'Info: Setting up TD.'
      call write_info(1)

      ! store old states
      aux_psi => sys%st%R_FUNC(psi)
      nullify(sys%st%R_FUNC(psi))
      aux_i1 = sys%st%st_start
      aux_i2 = sys%st%st_end

      ! initialize structure
      allocate(td)
      ! this allocates zpsi
      call td_init(td, sys, sys%m, sys%st)

    case(I_END_TD)
      message(1) = 'Info: Cleaning up TD.'
      call write_info(1)

      call td_end(td)
      deallocate(td); nullify(td)
      deallocate(sys%st%zpsi); nullify(sys%st%zpsi)
      
    case(I_INIT_ZPSI)
      message(1) = 'Info: Initializing zpsi.'
      call write_info(1)

      do i = sys%st%st_start, sys%st%st_end
        sys%st%zpsi(:,:, i,:) = aux_psi(:,:, i,:)
      end do

      message(1) = 'Info: Deallocating rpsi.'
      call write_info(1)
      deallocate(aux_psi) ! clean up old wf
      nullify(aux_psi)

    case(I_LOAD_ZPSI)
      message(1) = 'Info: Loading zpsi.'
      call write_info(1)

      if(zstates_load_restart(trim(td%filename), &
           sys%m, sys%st, iter=td%iter, v1=td%v_old1, v2=td%v_old2)) then

        !TODO -> ions

        ! define density and hamiltonian
        call zcalcdens(sys%st, sys%m%np, sys%st%rho, reduce=.true.)
        call hamiltonian_setup(h, sys)
        call zhamiltonian_eigenval (h, sys, 1, sys%st%nst)
        call hamiltonian_energy(h, sys, -1, reduce=.true.)        
        
      else
        ! reset variables
        sys%st%st_start = aux_i1
        sys%st%st_end   = aux_i2

        i_stack(instr) = I_INIT_ZPSI
        instr = instr + 1; i_stack(instr) = I_SETUP_TD
        instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
        instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
        instr = instr + 1; i_stack(instr) = I_END_TD
        cycle program
      end if

    case(I_TD)
      message(1) = 'Info: Time-dependent simulation.'
      call write_info(1)

      call td_run(td, unoccv%st, sys, h)

    case(I_SETUP_OCC_AN)
      message(1) = 'Info: Seting up occupational analysis.'
      call write_info(1)

      if(td%occ_analysis) then
        instr = instr - 1
        call m_load_unocc()
        cycle program
      end if

    case(I_END_OCC_AN)
      message(1) = 'Info: Cleaning up occupational analysis.'
      call write_info(1)

      if(td%occ_analysis) then
        i_stack(instr) = I_END_UNOCC
        cycle program
      end if

    case(I_START_POL)
      ! just delete the pol file
      message(1) = 'Info: Starting static polarizability calculation'
      call write_info(1)

      call io_assign(iunit)
      open(iunit, file=trim(sys%sysname)//'.pol_restart', status='unknown')
      write(iunit, '(a)') ' '
      call io_close(iunit)
        
    case(I_POL_SCF)
      message(1) = 'Info: Calculating static polarizability'
      call write_info(1)

      call scf_init(scfv, sys)
      call static_pol_run(scfv, sys, h)
      call scf_end(scfv)

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

  call oct_parse_int(C_string('CalculationMode'), 1, calc_mode)
  if( (calc_mode < 1 .or. calc_mode > 9) .and. (calc_mode .ne. M_PULPO_A_FEIRA)) then
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
    message(11)= '                   =99 <= prints out the "Pulpo a Feira" recipe'
    call write_fatal(11)
  end if

  if(calc_mode .ne. M_PULPO_A_FEIRA) then
    call units_init()
    call system_init(sys)
    call hamiltonian_init(h, sys)
    if(h%classic_pot) call generate_classic_pot(h, sys)
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
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
    instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
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
  case(M_PULPO_A_FEIRA)
    instr = instr + 1; i_stack(instr) = I_PULPO
  end select

end subroutine define_run_modes

end module run_prog
