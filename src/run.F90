#include "config.h"

module run_prog
use global
use fdf
use units
use states
use system
use hamiltonian
use lcao
use scf
use timedep
use pulpo

implicit none

type(system_type) :: sys
type(hamiltonian_type) :: h
integer :: calc_mode

! run stack
integer, private :: i_stack(100), instr
integer, private, parameter :: &
     M_START_STATIC_CALC             = 1, &
     M_RESUME_STATIC_CALC            = 2, &
     M_CALCULATE_UNOCC_STATES        = 3, &
     M_RESUME_UNOCC_STATES           = 4, &
     M_START_TD                      = 5, &
     M_RESUME_TD                     = 6, &
     M_START_STATIC_POL              = 7, &
     M_RESUME_STATIC_POL             = 8, &
     M_BO_MD                         = 9, &
     M_PULPO_A_FEIRA                 = 99

integer, private, parameter :: &
     I_SETUP_RPSI         =  1,  &
     I_END_RPSI           =  2,  &
     I_RANDOMIZE_RPSI     =  3,  &
     I_LOAD_RPSI          =  4,  &
     I_SETUP_HAMILTONIAN  =  5,  &
     I_SCF                =  6,  &
     I_LCAO               =  7,  &

     I_SETUP_OCC_RPSI     = 11,  &
     I_END_OCC_RPSI       = 12,  &
     I_RANDOMIZE_OCC_RPSI = 13,  &
     I_LOAD_OCC_RPSI      = 14,  &
     I_OCC_SCF            = 15,  &
    
     I_SETUP_TD           = 21,  &
     I_END_TD             = 22,  &
     I_INIT_ZPSI          = 23,  &
     I_LOAD_ZPSI          = 24,  &
     I_TD                 = 25,  &
     I_SETUP_OCC_AN       = 26,  &
     I_END_OCC_AN         = 27,  &

     I_PULPO              = 99

contains

subroutine run()
  type(td_type), pointer :: td
  R_TYPE, pointer :: aux_psi(:,:,:,:)
  integer :: i, aux_i1, aux_i2

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

    case(I_LCAO)
      if(fdf_boolean("LCAOStart", .true.)) then
        message(1) = 'Info: Performing LCAO calculation.'
        call write_info(1)
        call lcao_wf(sys, h)
      end if

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
      call scf_run(sys, h)

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

      call td_run(td, sys, h)
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

  calc_mode = fdf_integer('CalculationMode', 1)
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
  case(M_PULPO_A_FEIRA)
    instr = instr + 1; i_stack(instr) = I_PULPO
  end select

end subroutine define_run_modes

end module run_prog
