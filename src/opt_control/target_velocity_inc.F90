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
!! $Id: target_velocity_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_velocity(gr, geo, tg, oct, td, ep)
    type(grid_t),     intent(in)    :: gr
    type(geometry_t), intent(in)    :: geo
    type(target_t),   intent(inout) :: tg
    type(oct_t),      intent(in)    :: oct
    type(td_t),       intent(in)    :: td
    type(epot_t),     intent(inout) :: ep

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
    !% while "n_atom" is the respective atom number, corresponding to the 
    !% <tt>Coordinates</tt> block and "vec_comp" is the corresponding
    !% vector component of the velocity. The target string can be
    !% supplied by using several lines in the OCTTargetOperator block.
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
    !% <tt>OCTScheme = oct_algorithm_cg</tt> then you must supply 
    !% the target in terms of the ionic velocities AND the derivatives
    !% of the target with respect to the ionic velocity components.
    !% The derivatives are supplied via strings through the block
    !% <tt>OCTVelocityDerivatives</tt>.
    !% Each velocity component is supplied by <tt>"v[n_atom,vec_comp]"</tt>,
    !% while "n_atom" is the atom number, corresponding to the 
    !% <tt>Coordinates</tt> block and "vec_comp" is the corresponding
    !% vector component of the velocity. The first line of the 
    !% <tt>OCTVelocityDerivatives</tt> block contains the derivatives
    !% with respect to "v[1,*]", the second with respect to "v[2,*]" and so
    !% on. The first column contains all derivatives with respect "v[*,1]",
    !% the second with respect to "v[*,2]" and the third w.r.t. "v[*,3]".
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
       
    if(parse_block(datasets_check('OCTVelocityTarget'),blk)==0) then
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
       
    tg%move_ions = ion_dynamics_ions_move(td%ions)
       
    if(oct%algorithm  ==  oct_algorithm_cg) then
      if(parse_block(datasets_check('OCTVelocityDerivatives'),blk)==0) then
        SAFE_ALLOCATE(tg%vel_der_array(1:geo%natoms,1:gr%sb%dim))
        do ist=0, geo%natoms-1
          do jst=0, gr%sb%dim-1
            call parse_block_string(blk, ist, jst, tg%vel_der_array(ist+1, jst+1))
          end do
        end do
        call parse_block_end(blk)
      else
        message(1) = 'If OCTTargetOperator = oct_tg_velocity, and'
        message(2) = 'OCTScheme = oct_algorithm_cg, then you must define the'
        message(3) = 'blocks "OCTVelocityTarget" AND "OCTVelocityDerivatives"'
        call messages_fatal(3)
      end if
    end if
          
      SAFE_ALLOCATE(tg%grad_local_pot(1:geo%natoms, 1:gr%mesh%np, 1:gr%sb%dim))
      SAFE_ALLOCATE(vl(1:gr%mesh%np_part))
      SAFE_ALLOCATE(vl_grad(1:gr%mesh%np, 1:gr%sb%dim))
      SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))

      ! calculate gradient of each species potential
      do iatom=1, geo%natoms
        vl(:) = M_ZERO
        vl_grad(:,:) = M_ZERO
        call epot_local_potential(ep, gr%der, gr%dgrid, geo, iatom, vl)
        call dderivatives_grad(gr%der, vl, vl_grad)
        forall(ist=1:gr%mesh%np, jst=1:gr%sb%dim)
          tg%grad_local_pot(iatom, ist, jst) = vl_grad(ist, jst)
        end forall
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
    type(target_t),   intent(inout) :: tg
    type(oct_t), intent(in)       :: oct
    PUSH_SUB(target_end_velocity)
    if(oct%algorithm  ==  oct_algorithm_cg) then
      SAFE_DEALLOCATE_P(tg%vel_der_array)
      SAFE_DEALLOCATE_P(tg%grad_local_pot)
      SAFE_DEALLOCATE_P(tg%rho)
     end if
     SAFE_DEALLOCATE_P(tg%td_fitness)
    POP_SUB(target_end_velocity)
  end subroutine target_end_velocity


  ! ----------------------------------------------------------------------
  subroutine target_output_velocity(tg, gr, dir, geo, outp)
    type(target_t), intent(inout) :: tg
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp

    PUSH_SUB(target_output_velocity)
    
    call loct_mkdir(trim(dir))
    call output_states(tg%st, gr, geo, trim(dir), outp)

    POP_SUB(target_output_velocity)
  end subroutine target_output_velocity
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_velocity(tg, geo) result(j1)
    type(target_t),   intent(inout) :: tg
    type(geometry_t), intent(in)    :: geo

    integer :: i
    FLOAT :: f_re, dummy(3)
    FLOAT, allocatable :: x(:, :)
    character(len=4096) :: inp_string
    PUSH_SUB(target_j1_velocity)

    SAFE_ALLOCATE(x(1:geo%natoms, 1:geo%space%dim))
    forall(i=1: geo%natoms) x(i, 1:geo%space%dim) = geo%atom(i)%v(1:geo%space%dim)

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
  subroutine target_chi_velocity(gr, tg, chi_out, geo)
    type(grid_t),     intent(inout) :: gr
    type(target_t),   intent(inout) :: tg
    type(states_t),   intent(inout) :: chi_out
    type(geometry_t), intent(in)    :: geo

    integer :: ip, idim, ist, jst, ik
    character(len=1024) :: temp_string
    FLOAT :: df_dv, dummy(3)
    FLOAT, allocatable :: x(:, :)
    PUSH_SUB(target_chi_velocity)

    !we have a time-dependent target --> Chi(T)=0
    forall(ip=1:gr%mesh%np, idim=1:chi_out%d%dim, ist=chi_out%st_start:chi_out%st_end, ik=1:chi_out%d%nik)
       chi_out%zpsi(ip, idim, ist, ik) = M_z0
    end forall

    SAFE_ALLOCATE(x(1:geo%natoms, 1:geo%space%dim))
    forall(ip=1: geo%natoms) x(ip, 1:geo%space%dim) = geo%atom(ip)%v(1:geo%space%dim)
      
    !calculate dF/dn, which is the time-independent part of the inhomogenous term for the propagation of Chi
    df_dv = M_ZERO
    dummy(:) = M_ZERO
    tg%rho(:) = M_ZERO
    do ist=1, geo%natoms
      do jst=1, gr%sb%dim
        temp_string = tg%vel_der_array(ist, jst)
        call parse_array(temp_string, x, 'v')
        call conv_to_C_string(temp_string)
        call parse_expression(df_dv, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), temp_string)
        tg%rho(:) = tg%rho(:) + df_dv*tg%grad_local_pot(ist,:,jst)/species_weight(geo%atom(ist)%spec)
      end do
    end do

    SAFE_DEALLOCATE_A(x)
    POP_SUB(target_chi_velocity)
  end subroutine target_chi_velocity


  ! ---------------------------------------------------------
  !> 
  !!
  subroutine target_tdcalc_velocity(tg, hm, gr, geo, psi, time, max_time)
    type(target_t),      intent(inout) :: tg
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(states_t),      intent(inout) :: psi
    integer,             intent(in)    :: time
    integer,             intent(in)    :: max_time

    CMPLX, allocatable :: opsi(:, :)
    integer :: iatom, ik, ist, idim
    FLOAT :: dt
    PUSH_SUB(target_tdcalc_velocity)

    tg%td_fitness(time) = M_ZERO

    ! If the ions move, the target is computed in the propagation routine.
    if(.not.target_move_ions(tg)) then

      SAFE_ALLOCATE(opsi(1:gr%mesh%np_part, 1:1))
      opsi = M_z0
      ! WARNING This does not work for spinors.
      do iatom = 1, geo%natoms
        geo%atom(iatom)%f(1:gr%sb%dim) = hm%ep%fii(1:gr%sb%dim, iatom)
        do ik = 1, psi%d%nik
          do ist = 1, psi%nst
            do idim = 1, gr%sb%dim
              opsi(1:gr%mesh%np, 1) = tg%grad_local_pot(iatom, 1:gr%mesh%np, idim) * psi%zpsi(1:gr%mesh%np, 1, ist, ik)
              geo%atom(iatom)%f(idim) = geo%atom(iatom)%f(idim) + real(psi%occ(ist, ik) * &
                zmf_dotp(gr%mesh, psi%d%dim, opsi, psi%zpsi(:, :, ist, ik)), REAL_PRECISION)
            end do
          end do
        end do
      end do
      SAFE_DEALLOCATE_A(opsi)

    end if

    dt = tg%dt
    if( (time  ==  0) .or. (time  ==  max_time) ) dt = tg%dt * M_HALF
    do iatom = 1, geo%natoms
       geo%atom(iatom)%v(1:MAX_DIM) = geo%atom(iatom)%v(1:MAX_DIM) + &
         geo%atom(iatom)%f(1:MAX_DIM) * dt / species_weight(geo%atom(iatom)%spec)
    end do

    POP_SUB(target_tdcalc_velocity)
  end subroutine target_tdcalc_velocity
  ! ----------------------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
