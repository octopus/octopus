!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module initst_oct_m
  use density_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ions_oct_m
  use messages_oct_m
  use mesh_function_oct_m
  use opt_control_state_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use states_elec_oct_m
  use states_elec_restart_oct_m
  use string_oct_m
  use electrons_oct_m
  use td_oct_m
  use v_ks_oct_m
  use varinfo_oct_m
  use types_oct_m
  use xc_oct_m

  implicit none

  private
  public :: initial_state_init

  integer, parameter ::  &
    oct_is_groundstate      = 1,      &
    oct_is_excited          = 2,      &
    oct_is_gstransformation = 3,      &
    oct_is_userdefined      = 4         


contains


  ! ---------------------------------------------------------
  subroutine initial_state_init(sys, qcstate)
    type(electrons_t),                 intent(inout) :: sys
    type(opt_control_state_t), target, intent(inout) :: qcstate

    integer           :: ik, ib, idim, inst, inik, id, is, ip, ierr, &
                         no_states, istype, freeze_orbitals
    type(block_t)     :: blk
    FLOAT             :: xx(1:sys%space%dim), rr, psi_re, psi_im
    type(restart_t) :: restart
    CMPLX, allocatable :: zpsi(:, :)

    type(states_elec_t), pointer :: psi

    PUSH_SUB(initial_state_init)

    call opt_control_state_init(qcstate, sys%st, sys%ions)
    psi => opt_control_point_qs(qcstate)
    call states_elec_deallocate_wfns(psi)
    call states_elec_allocate_wfns(psi, sys%gr%mesh, TYPE_CMPLX)

    !%Variable OCTInitialState
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default oct_is_groundstate
    !%Description
    !% Describes the initial state of the quantum system.
    !% Possible arguments are:
    !%Option oct_is_groundstate 1
    !% Start in the ground state.
    !%Option oct_is_excited 2
    !% Currently not in use.
    !%Option oct_is_gstransformation 3
    !% Start in a transformation of the ground-state orbitals, as defined in the
    !% block <tt>OCTInitialTransformStates</tt>.
    !%Option oct_is_userdefined 4
    !% Start in a user-defined state.
    !%End
    call parse_variable(sys%namespace, 'OCTInitialState', oct_is_groundstate, istype)
    if(.not.varinfo_valid_option('OCTInitialState', istype)) call messages_input_error(sys%namespace, 'OCTInitialState')    

    select case(istype)
    case(oct_is_groundstate) 
      message(1) =  'Info: Using ground state for initial state.'
      call messages_info(1)
      call restart_init(restart, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
      if(ierr == 0) then
        call states_elec_load(restart, sys%namespace, sys%space, psi, sys%gr%mesh, sys%kpoints, ierr)
      end if
      if (ierr /= 0) then
        message(1) = "Unable to read wavefunctions."
        call messages_fatal(1)
      end if
      call restart_end(restart)

    case(oct_is_excited)  
      message(1) = 'Using an excited state as the starting state for an '
      message(2) = 'optimal-control run is not possible yet.'
      message(3) = 'Try using "OCTInitialState = oct_is_transformation" instead.'
      call messages_fatal(3)

    case(oct_is_gstransformation)   
      message(1) =  'Info: Using superposition of states for initial state.'
      call messages_info(1)


      !%Variable OCTInitialTransformStates
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTInitialState = oct_is_gstransformation</tt>, you must specify an
      !% <tt>OCTInitialTransformStates</tt> block, in order to specify which linear
      !% combination of the states present in <tt>restart/gs</tt> is used to
      !% create the initial state.
      !% 
      !% The syntax is the same as the <tt>TransformStates</tt> block.
      !%End

      if(.not. parse_is_defined(sys%namespace, "OCTInitialTransformStates")) then
        message(1) = 'If "OCTInitialState = oct_is_gstransformation", then you must'
        message(2) = 'supply an "OCTInitialTransformStates" block to define the transformation.'
        call messages_fatal(2)
      end if

      call restart_init(restart, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
      if(ierr /= 0) then
        message(1) = "Could not read states for OCTInitialTransformStates."
        call messages_fatal(1)
      end if
      
      call states_elec_transform(psi, sys%namespace, sys%space, restart, sys%gr%mesh, sys%kpoints, prefix = "OCTInitial")
      call restart_end(restart)

    case(oct_is_userdefined) 
      message(1) =  'Info: Building user-defined initial state.'
      call messages_info(1)
      
      !%Variable OCTInitialUserdefined
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% Define an initial state. Syntax follows the one of the <tt>UserDefinedStates</tt> block.
      !% Example:
      !%
      !% <tt>%OCTInitialUserdefined
      !% <br>&nbsp;&nbsp; 1 | 1 | 1 |  "exp(-r^2)*exp(-i*0.2*x)"
      !% <br>%</tt>
      !%  
      !%End
      if(parse_block(sys%namespace, 'OCTInitialUserdefined', blk) == 0) then

        SAFE_ALLOCATE(zpsi(1:sys%gr%mesh%np, 1:psi%d%dim))
        
        no_states = parse_block_n(blk)
        do ib = 1, no_states
          call parse_block_integer(blk, ib - 1, 0, idim)
          call parse_block_integer(blk, ib - 1, 1, inst)
          call parse_block_integer(blk, ib - 1, 2, inik)

          ! read formula strings and convert to C strings
          do id = 1, psi%d%dim
            do is = 1, psi%nst
              do ik = 1, psi%d%nik   
                
                ! does the block entry match and is this node responsible?
                if(.not. (id  ==  idim .and. is  ==  inst .and. ik  ==  inik    &
                  .and. psi%st_start  <=  is .and. psi%st_end >= is) ) cycle
                
                ! parse formula string
                call parse_block_string(                            &
                  blk, ib - 1, 3, psi%user_def_states(id, is, ik))
                ! convert to C string
                call conv_to_C_string(psi%user_def_states(id, is, ik))
                
                do ip = 1, sys%gr%mesh%np
                  xx = sys%gr%mesh%x(ip, :)
                  rr = sqrt(sum(xx**2))
                  
                  ! parse user-defined expressions
                  call parse_expression(psi_re, psi_im, &
                    sys%space%dim, xx, rr, M_ZERO, psi%user_def_states(id, is, ik))
                  ! fill state
                  zpsi(ip, id) = TOCMPLX(psi_re, psi_im)
                end do
                call states_elec_set_state(psi, sys%gr%mesh, id, is, ik, zpsi(:, id))
              end do
            end do
          end do
        end do
        call parse_block_end(blk)
        do ik = 1, psi%d%nik
          do is = psi%st_start, psi%st_end
            call states_elec_get_state(psi, sys%gr%mesh, is, ik, zpsi)
            call zmf_normalize(sys%gr%mesh, psi%d%dim, zpsi)
            call states_elec_set_state(psi, sys%gr%mesh, is, ik, zpsi)
          end do
        end do
        SAFE_DEALLOCATE_A(zpsi)
      else
        call messages_variable_is_block(sys%namespace, 'OCTInitialUserdefined')
      end if
      
    case default
      write(message(1),'(a)') "No valid initial state defined."
      write(message(2),'(a)') "Choosing the ground state."
      call messages_info(2)
    end select

    ! Check whether we want to freeze some of the deeper orbitals.
    call parse_variable(sys%namespace, 'TDFreezeOrbitals', 0, freeze_orbitals)
    if(freeze_orbitals > 0) then
      ! In this case, we first freeze the orbitals, then calculate the Hxc potential.
      call states_elec_freeze_orbitals(psi, sys%namespace, sys%gr, sys%mc, sys%kpoints, &
                   freeze_orbitals, family_is_mgga(sys%ks%xc_family))
      write(message(1),'(a,i4,a,i4,a)') 'Info: The lowest', freeze_orbitals, &
        ' orbitals have been frozen.', psi%nst, ' will be propagated.'
      call messages_info(1)
      call density_calc(psi, sys%gr, psi%rho)
      call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions, calc_eigenval = .true.)
    elseif(freeze_orbitals < 0) then
      ! This means SAE approximation. We calculate the Hxc first, then freeze all
      ! orbitals minus one.
      write(message(1),'(a)') 'Info: The single-active-electron approximation will be used.'
      call messages_info(1)
      call density_calc(psi, sys%gr, psi%rho)
      call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions, calc_eigenval = .true.)
      call states_elec_freeze_orbitals(psi, sys%namespace, sys%gr, sys%mc, sys%kpoints, &
                   psi%nst - 1, family_is_mgga(sys%ks%xc_family))
      call v_ks_freeze_hxc(sys%ks)
      call density_calc(psi, sys%gr, psi%rho)
    else
      ! Normal run.
      call density_calc(psi, sys%gr, psi%rho)
      call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions, calc_eigenval = .true.)
    end if
    
    POP_SUB(initial_state_init)
  end subroutine initial_state_init

end module initst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
