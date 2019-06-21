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

    integer :: is, err, idir, ik, ib
    character(len=MAX_PATH_LEN) :: fname
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

      if(hm%ep%classical_pot > 0) then
        call dio_function_output(outp%how, dir, "vc", der%mesh, hm%ep%Vclassical, units_out%energy, err, geo = geo, grp = grp)
      end if

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
          if (hm%ep%classical_pot > 0) then
            call dio_function_output(outp%how, dir, fname, der%mesh, &
              potential + hm%ep%Vclassical, units_out%energy, err, geo = geo, grp = grp)
          else
            call dio_function_output(outp%how, dir, fname, der%mesh, &
                potential, units_out%energy, err, geo = geo, grp = grp)
          end if
        end do
        SAFE_DEALLOCATE_A(potential)
      end if

      !PCM potentials
      if(hm%theory_level == KOHN_SHAM_DFT .and. hm%pcm%run_pcm ) then
        if (hm%pcm%solute .and. hm%pcm%localf) then
          call dio_function_output(outp%how, dir, 'vpcm', der%mesh, hm%pcm%v_e_rs + hm%pcm%v_n_rs + hm%pcm%v_ext_rs , & 
                                   units_out%energy, err, geo = geo, grp = grp)
          call dio_function_output(outp%how, dir, 'vpcm_sol', der%mesh, hm%pcm%v_e_rs + hm%pcm%v_n_rs , & 
                                   units_out%energy, err, geo = geo, grp = grp)
          call dio_function_output(outp%how, dir, 'vpcm_e', der%mesh, hm%pcm%v_e_rs , & 
                                   units_out%energy, err, geo = geo, grp = grp)
          call dio_function_output(outp%how, dir, 'vpcm_n', der%mesh, hm%pcm%v_n_rs , & 
                                   units_out%energy, err, geo = geo, grp = grp)
          call dio_function_output(outp%how, dir, 'vpcm_ext', der%mesh, hm%pcm%v_ext_rs , & 
                                   units_out%energy, err, geo = geo, grp = grp)
        else if (hm%pcm%solute .and. .not.hm%pcm%localf) then
          call dio_function_output(outp%how, dir, 'vpcm_sol', der%mesh, hm%pcm%v_e_rs + hm%pcm%v_n_rs , & 
                                   units_out%energy, err, geo = geo, grp = grp)
          call dio_function_output(outp%how, dir, 'vpcm_e', der%mesh, hm%pcm%v_e_rs , & 
                                   units_out%energy, err, geo = geo, grp = grp)
          call dio_function_output(outp%how, dir, 'vpcm_n', der%mesh, hm%pcm%v_n_rs , & 
                                   units_out%energy, err, geo = geo, grp = grp)
        else if (.not.hm%pcm%solute .and. hm%pcm%localf) then
          call dio_function_output(outp%how, dir, 'vpcm_ext', der%mesh, hm%pcm%v_ext_rs , & 
                                   units_out%energy, err, geo = geo, grp = grp)
        end if
      end if

      if(hm%self_induced_magnetic) then
        ! unit of magnetic field is same as of electric field, and same as force (since e = 1)
        select case(der%mesh%sb%dim)
        case(3)
          do idir = 1, der%mesh%sb%dim
            call dio_function_output(outp%how, dir, 'Bind_'//index2axis(idir), der%mesh, hm%b_ind(:, idir), &
              units_out%force, err, geo = geo, grp = grp)
          end do
        case(2)
          call dio_function_output(outp%how, dir, 'Bind_z', der%mesh, hm%b_ind(:, 1), units_out%force, err, geo = geo, grp = grp)
        end select
      end if
    end if


    if(bitand(outp%what, OPTION__OUTPUT__XC_DENSITY) /= 0 .and. hm%theory_level /= INDEPENDENT_PARTICLES) then
      SAFE_ALLOCATE(v0(1:der%mesh%np_part, 1))
      SAFE_ALLOCATE(nxc(1:der%mesh%np))

      do is = 1, min(hm%d%ispin, 2)
        if(hm%d%ispin == 1) then
          write(fname, '(a)') 'nxc'
        else
          write(fname, '(a,i1)') 'nxc-sp', is
        end if
                
        v0(1:der%mesh%np, 1) = hm%vxc(1:der%mesh%np, is)

        call dderivatives_lapl(der, v0(:, 1), nxc)

        call dio_function_output(outp%how, dir, fname, der%mesh, nxc, units_out%energy, err, geo = geo, grp = grp)
        
      end do

      SAFE_DEALLOCATE_A(v0)
      SAFE_DEALLOCATE_A(nxc)
    end if

    if(bitand(outp%what, OPTION__OUTPUT__CURRENT) /= 0) then
      if(states_are_complex(st)) then
        do is = 1, hm%d%nspin

          if(st%d%nspin == 1) then
            write(fname, '(2a)') 'current'
          else
            write(fname, '(a,i1)') 'current-sp', is
          end if
          
          call io_function_output_vector(outp%how, dir, fname, der%mesh, &
            st%current(:, :, is), der%mesh%sb%dim, (unit_one/units_out%time)*units_out%length**(1 - der%mesh%sb%dim), err, &
            geo = geo, grp = st%dom_st_kpt_mpi_grp, vector_dim_labels = (/'x', 'y', 'z'/))

        end do
      else
        message(1) = 'No current density output for real states since it is identically zero.'
        call messages_warning(1)
      end if
    end if
   
    if(bitand(outp%whatBZ, OPTION__OUTPUT_KPT__CURRENT_KPT) /= 0) then
      if(states_are_complex(st)) then
      
        SAFE_ALLOCATE(current_kpt(st%d%kpt%start:st%d%kpt%end, der%mesh%sb%dim)) 
        do ik = st%d%kpt%start,st%d%kpt%end
          do idir = 1,der%mesh%sb%dim
            current_kpt(ik,idir) = dmf_integrate(der%mesh, st%current_kpt(:, idir, ik))
          end do
        end do
        write(fname, '(2a)') 'current_kpt'
        call io_function_output_vector_BZ(outp%how, dir, fname, der%mesh, st%d%kpt, &
            current_kpt(:, :), der%mesh%sb%dim, (unit_one/units_out%time)*units_out%length**(1 - der%mesh%sb%dim), err, &
            grp = st%st_kpt_mpi_grp, vector_dim_labels = (/'x', 'y', 'z'/))
        SAFE_DEALLOCATE_A(current_kpt)
      else
        message(1) = 'No current density output for real states since it is identically zero.'
        call messages_warning(1)
      end if
    end if

    if(bitand(outp%what, OPTION__OUTPUT__HEAT_CURRENT) /= 0) then

      SAFE_ALLOCATE(heat_current(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:st%d%nspin))
      
      call current_heat_calculate(der, hm, geo, st, heat_current)
      
      do is = 1, hm%d%nspin
        if(st%d%nspin == 1) then
          write(fname, '(2a)') 'heat_current'
        else
          write(fname, '(a,i1)') 'heat_current-sp', is
        end if
        
        call io_function_output_vector(outp%how, dir, fname, der%mesh, &
          st%current(:, :, is), der%mesh%sb%dim, (unit_one/units_out%time)*units_out%length**(1 - der%mesh%sb%dim), err, &
          geo = geo, grp = st%dom_st_kpt_mpi_grp, vector_dim_labels = (/'x', 'y', 'z'/))
        
        SAFE_DEALLOCATE_A(heat_current)
      end do
    end if
    
    if(bitand(outp%whatBZ, OPTION__OUTPUT_KPT__DENSITY_KPT) /= 0) then
      SAFE_ALLOCATE(density_kpt(1:st%d%nik))
      density_kpt(1:st%d%nik) = M_ZERO

      SAFE_ALLOCATE(density_tmp(1:gr%fine%mesh%np, st%d%nspin))

      !These two conditions should be copied from the density_calc_end routine
      !We cannot call this routine as we must not symmetrize of reduce on kpoints
      ASSERT(.not.states_are_packed(st))
      ASSERT(.not.gr%have_fine_mesh)

      do ik = st%d%kpt%start,st%d%kpt%end
        call density_calc_init(dens_calc, st, gr, density_tmp)
        do ib = st%group%block_start, st%group%block_end
          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))
        end do
 
        density_kpt(ik) = M_ZERO
        do is = 1, st%d%nspin
          density_kpt(ik) = density_kpt(ik) + dmf_integrate(der%mesh, density_tmp(:,is))
        end do
      end do

   #if defined(HAVE_MPI)        
      if(st%parallel_in_states .or. st%d%kpt%parallel) then
        call comm_allreduce(st%st_kpt_mpi_grp%comm, density_kpt)
      end if
   #endif  

      call io_function_output_global_BZ(outp%how, dir, "density_kpt", gr%mesh, density_kpt, unit_one, err)
      SAFE_DEALLOCATE_A(density_tmp)
      SAFE_DEALLOCATE_A(density_kpt)
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
