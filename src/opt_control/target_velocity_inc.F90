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


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_velocity(gr, namespace, ions, tg, oct, td, ep)
    type(grid_t),      intent(in)    :: gr
    type(namespace_t), intent(in)    :: namespace
    type(ions_t),      intent(in)    :: ions
    type(target_t),    intent(inout) :: tg
    type(oct_t),       intent(in)    :: oct
    type(td_t),        intent(in)    :: td
    type(epot_t),      intent(in)    :: ep

    integer             :: iatom, ist, jst, jj
    FLOAT, allocatable  :: vl(:), vl_grad(:,:)
    type(block_t)       :: blk
    character(len=1024) :: expression

    PUSH_SUB(target_init_velocity)

    !%Variable OCTVelocityTarget
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% If <tt>OCTTargetOperator = oct_tg_velocity</tt>, then one must supply the 
    !% target to optimize in terms of the ionic velocities. This is done by 
    !% supplying a string through the block <tt>OCTVelocityTarget</tt>.
    !% Each velocity component is supplied by <tt>"v[n_atom,vec_comp]"</tt>,
    !% where <tt>n_atom</tt> is the atom number, corresponding to the 
    !% <tt>Coordinates</tt> block, and <tt>vec_comp</tt> is the corresponding
    !% vector component of the velocity. The target string can be
    !% supplied by using several lines in this block.
    !% As an example, the following target can be used to maximize the
    !% velocity difference between atom 1 and 2 (in a 3D system):
    !%
    !% <tt>%OCTVelocityTarget
    !% <br> "(v[1,1]-v[2,1])^2 + (v[1,2]-v[2,2])^2 + "
    !% <br> "(v[1,3]-v[2,3])^2"
    !% <br>%</tt>
    !%
    !%End

    !%Variable OCTVelocityDerivatives
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% If <tt>OCTTargetOperator = oct_tg_velocity</tt>, and
    !% <tt>OCTScheme = oct_cg</tt> or <tt>OCTScheme = oct_bfgs</tt>
    !% then you must supply the target in terms of the ionic velocities AND
    !% the derivatives of the target with respect to the ionic velocity components.
    !% The derivatives are supplied via strings through the block
    !% <tt>OCTVelocityDerivatives</tt>.
    !% Each velocity component is supplied by <tt>"v[n_atom,vec_comp]"</tt>,
    !% while <tt>n_atom</tt> is the atom number, corresponding to the 
    !% <tt>Coordinates</tt> block, and <tt>vec_comp</tt> is the corresponding
    !% vector component of the velocity. The first line of the 
    !% <tt>OCTVelocityDerivatives</tt> block contains the derivatives
    !% with respect to <tt>v[1,*]</tt>, the second with respect to <tt>v[2,*]</tt> and so
    !% on. The first column contains all derivatives with respect <tt>v[*,1]</tt>,
    !% the second with respect to <tt>v[*,2]</tt> and the third w.r.t. <tt>v[*,3]</tt>.
    !% As an example, we show the <tt>OCTVelocityDerivatives</tt> block
    !% corresponding to the target shown in the <tt>OCTVelocityTarget</tt>
    !% help section:
    !%
    !% <tt>%OCTVelocityDerivatives
    !% <br> " 2*(v[1,1]-v[2,1])" | " 2*(v[1,2]-v[2,2])" | " 2*(v[1,3]-v[2,3])"
    !% <br> "-2*(v[1,1]-v[2,1])" | "-2*(v[1,2]-v[2,2])" | "-2*(v[1,3]-v[2,3])"
    !% <br>%</tt>
    !%
    !%End
       
    if(parse_block(namespace, 'OCTVelocityTarget', blk)==0) then
      tg%vel_input_string = " "
      do jj=0, parse_block_n(blk)-1
        call parse_block_string(blk, jj, 0, expression)
        tg%vel_input_string = trim(tg%vel_input_string) // trim(expression)
      end do
      call parse_block_end(blk)
    else
      message(1) = 'If OCTTargetOperator = oct_tg_velocity, then you must give the shape'
      message(2) = 'of this target in the block "OCTVelocityTarget".'
      call messages_fatal(2)
    end if
       
    tg%move_ions = ion_dynamics_ions_move(td%ions_dyn)
    if(tg%move_ions) then
      message(1) = 'If OCTTargetOperator = oct_tg_velocity, then you must not allow the ions'
      message(2) = 'to move. If you want to move the ions, then you can get the same functionality'
      message(3) = 'with OCTTargetOperator = oct_tg_classical.'
      call messages_fatal(3)
    end if
       
    if(oct%algorithm  ==  OPTION__OCTSCHEME__OCT_CG .or. oct%algorithm == OPTION__OCTSCHEME__OCT_BFGS) then
      if(parse_block(namespace, 'OCTVelocityDerivatives', blk)==0) then
        SAFE_ALLOCATE(tg%vel_der_array(1:ions%natoms,1:gr%sb%dim))
        do ist=0, ions%natoms-1
          do jst=0, gr%sb%dim-1
            call parse_block_string(blk, ist, jst, tg%vel_der_array(ist+1, jst+1))
          end do
        end do
        call parse_block_end(blk)
      else
        message(1) = 'If OCTTargetOperator = oct_tg_velocity, and'
        message(2) = 'OCTScheme = oct_cg, or OCTScheme = oct_bfgs then you must define the'
        message(3) = 'blocks "OCTVelocityTarget" AND "OCTVelocityDerivatives"'
        call messages_fatal(3)
      end if
    end if

      SAFE_ALLOCATE(tg%grad_local_pot(1:ions%natoms, 1:gr%mesh%np, 1:gr%sb%dim))
      SAFE_ALLOCATE(vl(1:gr%mesh%np_part))
      SAFE_ALLOCATE(vl_grad(1:gr%mesh%np, 1:gr%sb%dim))
      SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))

      ! calculate gradient of each species potential
      do iatom = 1, ions%natoms
        vl(:) = M_ZERO
        vl_grad(:,:) = M_ZERO
        call epot_local_potential(ep, namespace, ions%space, ions%latt, gr%mesh, ions%atom(iatom)%species, &
          ions%pos(:, iatom), iatom, vl)
        call dderivatives_grad(gr%der, vl, vl_grad)
        do ist = 1, gr%mesh%np
          do jst=1, gr%sb%dim
            tg%grad_local_pot(iatom, ist, jst) = vl_grad(ist, jst)
          end do
        end do
      end do
      SAFE_DEALLOCATE_A(vl)
      SAFE_DEALLOCATE_A(vl_grad)

      ! Note that the calculation of the gradient of the potential
      ! is wrong at the borders of the box, since it assumes zero boundary
      ! conditions. The best way to solve this problems is to define the 
      ! target making use of the definition of the forces based on the gradient
      ! of the density, rather than on the gradient of the potential.

    tg%dt = td%dt
    SAFE_ALLOCATE(tg%td_fitness(0:td%max_iter))
    tg%td_fitness = M_ZERO

    POP_SUB(target_init_velocity)
  end subroutine target_init_velocity


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_velocity(tg, oct)
    type(target_t), intent(inout) :: tg
    type(oct_t),    intent(in)    :: oct

    PUSH_SUB(target_end_velocity)

    if(oct%algorithm  ==  OPTION__OCTSCHEME__OCT_CG .or. oct%algorithm  ==  OPTION__OCTSCHEME__OCT_BFGS) then
      SAFE_DEALLOCATE_A(tg%vel_der_array)
      SAFE_DEALLOCATE_A(tg%grad_local_pot)
      SAFE_DEALLOCATE_A(tg%rho)
     end if
     SAFE_DEALLOCATE_A(tg%td_fitness)

     POP_SUB(target_end_velocity)
  end subroutine target_end_velocity


  ! ----------------------------------------------------------------------
  subroutine target_output_velocity(tg, namespace, space, gr, dir, ions, hm, outp)
    type(target_t),      intent(in) :: tg
    type(namespace_t),   intent(in) :: namespace
    type(space_t),       intent(in) :: space
    type(grid_t),        intent(in) :: gr
    character(len=*),    intent(in) :: dir
    type(ions_t),        intent(in) :: ions
    type(hamiltonian_elec_t), intent(in) :: hm
    type(output_t),      intent(in) :: outp

    PUSH_SUB(target_output_velocity)
    
    call io_mkdir(trim(dir), namespace)
    call output_states(outp, namespace, space, trim(dir), tg%st, gr, ions, hm, -1)

    POP_SUB(target_output_velocity)
  end subroutine target_output_velocity
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_velocity(tg, ions) result(j1)
    type(target_t),   intent(in) :: tg
    type(ions_t),     intent(in) :: ions

    integer :: i
    FLOAT :: f_re, dummy(3)
    FLOAT, allocatable :: x(:, :)
    character(len=4096) :: inp_string
    PUSH_SUB(target_j1_velocity)

    SAFE_ALLOCATE(x(1:ions%natoms, 1:ions%space%dim))
    do i = 1, ions%natoms
      x(i, :) = ions%vel(:, i)
    end do

    f_re = M_ZERO
    dummy(:) = M_ZERO
    inp_string = tg%vel_input_string
    call parse_array(inp_string, x, 'v')
    call conv_to_C_string(inp_string)
    call parse_expression(f_re, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), inp_string)
    j1 = f_re

    SAFE_DEALLOCATE_A(x)
    POP_SUB(target_j1_velocity)
  end function target_j1_velocity


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_velocity(gr, tg, chi_out, ions)
    type(grid_t),        intent(in)    :: gr
    type(target_t),      intent(inout) :: tg
    type(states_elec_t), intent(inout) :: chi_out
    type(ions_t),        intent(in)    :: ions

    integer :: ip, ist, jst, ik, ib
    character(len=1024) :: temp_string
    FLOAT :: df_dv, dummy(3)
    FLOAT, allocatable :: x(:, :)
    PUSH_SUB(target_chi_velocity)

    !we have a time-dependent target --> Chi(T)=0
    do ik = chi_out%d%kpt%start, chi_out%d%kpt%end
      do ib = chi_out%group%block_start, chi_out%group%block_end
        call batch_set_zero(chi_out%group%psib(ib, ik))
      end do
    end do

    SAFE_ALLOCATE(x(1:ions%natoms, 1:ions%space%dim))
    do ip = 1, ions%natoms
      x(ip, :) = ions%vel(:, ip)
    end do

    !calculate dF/dn, which is the time-independent part of the inhomogenous term for the propagation of Chi
    df_dv = M_ZERO
    dummy(:) = M_ZERO
    tg%rho(:) = M_ZERO
    do ist = 1, ions%natoms
      do jst=1, gr%sb%dim
        temp_string = tg%vel_der_array(ist, jst)
        call parse_array(temp_string, x, 'v')
        call conv_to_C_string(temp_string)
        call parse_expression(df_dv, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), temp_string)
        tg%rho(:) = tg%rho(:) + df_dv*tg%grad_local_pot(ist,:,jst)/ions%mass(ist)
      end do
    end do

    SAFE_DEALLOCATE_A(x)
    POP_SUB(target_chi_velocity)
  end subroutine target_chi_velocity


  ! ---------------------------------------------------------
  !> 
  !!
  subroutine target_tdcalc_velocity(tg, hm, gr, ions, psi, time, max_time)
    type(target_t),      intent(inout) :: tg
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(grid_t),        intent(in)    :: gr
    type(ions_t),        intent(inout) :: ions
    type(states_elec_t), intent(in)    :: psi
    integer,             intent(in)    :: time
    integer,             intent(in)    :: max_time

    CMPLX, allocatable :: opsi(:, :), zpsi(:, :)
    integer :: iatom, ik, ist, idim
    FLOAT :: dt
    PUSH_SUB(target_tdcalc_velocity)

    tg%td_fitness(time) = M_ZERO

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np_part, 1:1))
    SAFE_ALLOCATE(opsi(1:gr%mesh%np_part, 1:1))
    opsi = M_z0
    ! WARNING This does not work for spinors.
    do iatom = 1, ions%natoms
      ions%tot_force(:, iatom) = hm%ep%fii(1:gr%sb%dim, iatom)
      do ik = 1, psi%d%nik
        do ist = 1, psi%nst
          do idim = 1, gr%sb%dim
            call states_elec_get_state(psi, gr%mesh, ist, ik, zpsi)
            opsi(1:gr%mesh%np, 1) = tg%grad_local_pot(iatom, 1:gr%mesh%np, idim)*zpsi(1:gr%mesh%np, 1)
            ions%tot_force(idim, iatom) = ions%tot_force(idim, iatom) &
              + TOFLOAT(psi%occ(ist, ik)*zmf_dotp(gr%mesh, psi%d%dim, opsi, zpsi))
          end do
        end do
      end do
    end do
    SAFE_DEALLOCATE_A(opsi)
    SAFE_DEALLOCATE_A(zpsi)
    
    dt = tg%dt
    if( (time  ==  0) .or. (time  ==  max_time) ) dt = tg%dt * M_HALF
    do iatom = 1, ions%natoms
      ions%vel(:, iatom) = ions%vel(:, iatom) + ions%tot_force(:, iatom) * dt / ions%mass(iatom)
    end do

    POP_SUB(target_tdcalc_velocity)
  end subroutine target_tdcalc_velocity
  ! ----------------------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
