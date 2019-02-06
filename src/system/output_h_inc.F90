!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch.
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

  ! ---------------------------------------------------------
  subroutine output_hamiltonian(hm, st, der, dir, outp, geo, gr, grp)
    type(hamiltonian_t),       intent(in)    :: hm
    type(states_t),            intent(inout) :: st
    type(derivatives_t),       intent(inout) :: der
    character(len=*),          intent(in)    :: dir
    type(output_t),            intent(in)    :: outp
    type(geometry_t),          intent(in)    :: geo
    type(grid_t),              intent(in)    :: gr
    type(mpi_grp_t), optional, intent(in)    :: grp !< the group that shares the same data, must contain the domains group

    integer :: is, err, idir, ispin, ik, ib
    character(len=MAX_PATH_LEN) :: fname
    FLOAT,         dimension(:),   pointer :: xpot
    FLOAT, allocatable :: v0(:,:), nxc(:), potential(:)
    FLOAT, allocatable :: current_kpt(:, :)
    FLOAT, allocatable :: density_kpt(:), density_tmp(:,:)
    type(density_calc_t) :: dens_calc
    FLOAT, allocatable :: gradvh(:, :), heat_current(:, :, :)

    PUSH_SUB(output_hamiltonian)
   

    if(bitand(outp%what, OPTION__OUTPUT__POTENTIAL) /= 0) then
      SAFE_ALLOCATE(v0(1:der%mesh%np, 1:hm%d%dim))
      v0(1:der%mesh%np, 1) = hm%ep%vpsl(1:der%mesh%np)
      call dio_function_output(outp%how, dir, "v0", der%mesh, v0(:, 1), units_out%energy, err, geo = geo, grp = grp)
      SAFE_DEALLOCATE_A(v0)

      if(associated(hm%ep%v_static)) then
        call dio_function_output(outp%how, dir, "vext", der%mesh, hm%ep%v_static, units_out%energy, err, geo = geo, grp = grp)
      end if

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        call dio_function_output(outp%how, dir, 'vh', der%mesh, hm%vhartree, units_out%energy, err, geo = geo, grp = grp)
        if(bitand(outp%what, OPTION__OUTPUT__POTENTIAL_GRADIENT) /= 0) then
          SAFE_ALLOCATE(gradvh(1:der%mesh%np, 1:der%mesh%sb%dim))
          call dderivatives_grad(der, hm%vhartree(1:der%mesh%np_part), gradvh(1:der%mesh%np, 1:der%mesh%sb%dim))
          call io_function_output_vector(outp%how, dir, 'grad_vh', der%mesh, gradvh(:, :), der%mesh%sb%dim, units_out%force, err,&
                   geo = geo, grp = grp, vector_dim_labels = (/'x', 'y', 'z'/))

          SAFE_ALLOCATE(v0(1:der%mesh%np_part, 1))
          v0(1:der%mesh%np, 1) = hm%ep%vpsl(1:der%mesh%np)
          call dderivatives_grad(der, v0(1:der%mesh%np_part, 1), gradvh(1:der%mesh%np, 1:der%mesh%sb%dim))
          call io_function_output_vector(outp%how, dir, 'grad_v0', der%mesh, gradvh(:, :), der%mesh%sb%dim, units_out%force, err,&
                   geo = geo, grp = grp, vector_dim_labels = (/'x', 'y', 'z'/))
          SAFE_DEALLOCATE_A(v0)
          SAFE_DEALLOCATE_A(gradvh)
        end if
        
        SAFE_ALLOCATE(potential(1:der%mesh%np))
        do is = 1, min(hm%d%ispin, 2)
          if(hm%d%ispin == 1) then
            write(fname, '(a)') 'vxc'
          else
            write(fname, '(a,i1)') 'vxc-sp', is
          end if
          call dio_function_output(outp%how, dir, fname, der%mesh, hm%vxc(:, is), units_out%energy, err, geo = geo, grp = grp)
          
          ! finally the full KS potential (without non-local PP contributions)
          potential = hm%ep%vpsl + hm%vhxc(:, is)
          if(hm%d%ispin == 1) then
            write(fname, '(a)') 'vks'
          else
            write(fname, '(a,i1)') 'vks-sp', is
          end if
          call dio_function_output(outp%how, dir, fname, der%mesh, potential, units_out%energy, err, geo = geo, grp = grp)
        end do
        SAFE_DEALLOCATE_A(potential)
      end if
    end if
 
    POP_SUB(output_hamiltonian)
  end subroutine output_hamiltonian


  ! ---------------------------------------------------------
  subroutine output_scalar_pot(outp, gr, geo, hm, dir, time)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(hamiltonian_t),  intent(inout) :: hm
    type(output_t),       intent(in)    :: outp
    character(len=*),     intent(in)    :: dir
    FLOAT, optional,      intent(in)    :: time

    integer :: is, err
    character(len=80) :: fname
    FLOAT, allocatable :: scalar_pot(:)

    PUSH_SUB(output_scalar_pot)

    if(bitand(outp%what, OPTION__OUTPUT__EXTERNAL_TD_POTENTIAL) /= 0) then
      SAFE_ALLOCATE(scalar_pot(1:gr%mesh%np))
      do is = 1, hm%ep%no_lasers
        write(fname, '(a,i1)') 'scalar_pot-', is
        scalar_pot = M_ZERO
        call laser_potential(hm%ep%lasers(is), gr%mesh, scalar_pot, time=time)
        call dio_function_output(outp%how, dir, fname, gr%mesh, scalar_pot, units_out%energy, err, geo = geo)
      end do
      SAFE_DEALLOCATE_A(scalar_pot)
    end if

    POP_SUB(output_scalar_pot)
  end subroutine output_scalar_pot


  ! ---------------------------------------------------------
  subroutine output_kick(outp, mesh, geo, kick, dir)
    type(mesh_t),     intent(in) :: mesh
    type(geometry_t), intent(in) :: geo
    type(kick_t),     intent(in) :: kick
    type(output_t),   intent(in) :: outp
    character(len=*), intent(in) :: dir

    integer :: err
    CMPLX, allocatable :: kick_function(:)
    
    PUSH_SUB(output_kick)

    if(bitand(outp%what, OPTION__OUTPUT__DELTA_PERTURBATION) /= 0) then
      SAFE_ALLOCATE(kick_function(1:mesh%np))
      call kick_function_get(mesh, kick, kick_function)
      call zio_function_output(outp%how, dir, "kick_function", mesh, kick_function(:), &
        units_out%energy, err, geo = geo)
      SAFE_DEALLOCATE_A(kick_function)
    end if
    
    POP_SUB(output_kick)
  end subroutine output_kick


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
