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
     I_LCAO               =  7

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
  call units_init()
  call system_init(sys)
  call hamiltonian_init(h, sys)
  call generate_external_pot(h, sys)

  calc_mode = fdf_integer('CalculationMode', 1)
  if(calc_mode < 1 .or. calc_mode > 9) then
    write(message(1), '(a,i2,a)') "Input: '", calc_mode, "' is not a valid CalculationMode"
    message(2) = '  Calculation Mode = 1 <= start static calculation'
    message(3) = '                   = 2 <= resume static calculation'
    message(5) = '                   = 3 <= calculate unocuppied states'
    message(6) = '                   = 4 <= resume unocuppied states calculation'
    message(7) = '                   = 5 <= start td'
    message(8) = '                   = 6 <= resume td'
    message(9) = '                   = 7 <= start static polarizability'
    message(9) = '                   = 8 <= resume static polarizability'
    message(10)= '                   = 9 <= perform Born-Oppenheimer MD'
    call write_fatal(11)
  end if
  
end subroutine run_init

subroutine run_end()
  call system_end(sys)
  call hamiltonian_end(h)
end subroutine run_end

subroutine define_run_modes()
  select case(calc_mode)
  case(1)
    instr = instr + 1; i_stack(instr) = I_END_RPSI
    instr = instr + 1; i_stack(instr) = I_SCF
    instr = instr + 1; i_stack(instr) = I_LCAO
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
    instr = instr + 1; i_stack(instr) = I_RANDOMIZE_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI
  case(2)
    instr = instr + 1; i_stack(instr) = I_END_RPSI
    instr = instr + 1; i_stack(instr) = I_SCF
    instr = instr + 1; i_stack(instr) = I_SETUP_HAMILTONIAN    
    instr = instr + 1; i_stack(instr) = I_LOAD_RPSI
    instr = instr + 1; i_stack(instr) = I_SETUP_RPSI  
  end select

end subroutine define_run_modes

end module run_prog
