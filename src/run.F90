#include "config.h"

module run_prog
use system
use hamiltonian

implicit none

type(system_type) :: sys
type(hamiltonian_type) :: h
integer :: calc_mode

! run stack
integer, private :: i_stack(100), instr
integer, private, parameter :: &
     I_SETUP_RPSI         =  1,  &
     I_END_RPSI           =  2,  &
     I_RANDOMIZE_RPSI     =  3,  &
     I_LOAD_RPSI          =  4,  &
     I_SETUP_HAMILTONIAN  =  5,  &
     I_SCF                =  6,  &
     I_LCAO               =  7,  &
     I_PULPO              = 99

integer, private, parameter :: M_START_STATIC_CALC             = 1, &
                               M_RESUME_STATIC_CALC            = 2, &
                               M_CALCULATE_UNOCC_STATES        = 3, &
                               M_RESUME_UNOCC_STATES           = 4, &
                               M_START_TD                      = 5, &
                               M_RESUME_TD                     = 6, &
                               M_START_STATIC_POL              = 7, &
                               M_RESUME_STATIC_POL             = 8, &
                               M_BO_MD                         = 9, &
                               M_PULPO_A_FEIRA                 = 99

contains

subroutine run()
  use global
  use units
  implicit none

  sub_name = 'run'; call push_sub()

  call run_init()

  instr = 0
  call define_run_modes()

  program: do
    if(instr <= 0) exit program

    select case(i_stack(instr))
    case(I_SETUP_RPSI)
      message(1) = 'Info: Allocating rpsi'
      call write_info(1)
        
      allocate(sys%st%R_FUNC(psi)(0:sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))

    case(I_END_RPSI)
      message(1) = 'Info: Deallocating rpsi'
      call write_info(1)

      if(associated(sys%st%R_FUNC(psi))) then
        deallocate(sys%st%R_FUNC(psi)); nullify(sys%st%R_FUNC(psi))
      end if

    case(I_RANDOMIZE_RPSI)
      message(1) = 'Info: Random generating starting wavefunctions'
      call write_info(1)

      call states_generate_random(sys%st, sys%m)

    case(I_PULPO)
      message(1) = ''
      message(2) = 'PULPO A FEIRA:'
      message(3) = ''
      message(4) = 'Ingredientes: Para 4 personas - 2 kg. de pulpo  - 2 kg. de patatas '
      message(5) = '- 100 grs. de pimentón picante - 100 grs. de sal gorda - aceite'
      message(6) = ''
      message(7) = 'Preparación: Se lava el pulpo en agua fría, se pone una olla de cobre'
      message(8) = 'con agua al fuego y cuando rompa a hervir se coge el pulpo, se mete y '
      message(9) = 'se saca del agua tres veces dejando que en cada intervalo vuelva a'
      message(10) = 'hervir el agua. Se deja cocer el pulpo durante unos 20 minutos '
      message(11) = 'retirándolo del fuego y dejándolo reposar durante 5 minutos. A'
      message(12) = 'continuación, se quita del agua y se corta en trozos finos con '
      message(13) = 'unas tijeras. Para servirlo se pone en unos platos de madera '
      message(14)= 'condimentándolo por este orden: sal, pimentón, aceite y se añaden '
      message(15)= 'unos cachelos (Patatas).'
      call write_info(15)

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
  case(M_PULPO_A_FEIRA)
    instr = instr + 1; i_stack(instr) = I_PULPO
  end select

end subroutine define_run_modes

end module run_prog
