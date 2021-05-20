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
  subroutine output_hamiltonian(outp, namespace, space, dir, hm, st, der, ions, gr, iter, grp)
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(hamiltonian_elec_t),  intent(in)    :: hm
    type(states_elec_t),       intent(inout) :: st
    type(derivatives_t),       intent(in)    :: der
    type(ions_t),              intent(in)    :: ions
    type(grid_t),              intent(in)    :: gr
    integer,                   intent(in)    :: iter
    type(mpi_grp_t), optional, intent(in)    :: grp !< the group that shares the same data, must contain the domains group

    integer :: is, err, idir, ik, ib
    character(len=MAX_PATH_LEN) :: fname
    FLOAT, allocatable :: vh(:), v0(:,:), nxc(:), potential(:)
    FLOAT, allocatable :: current_kpt(:, :)
    FLOAT, allocatable :: density_kpt(:, :), density_tmp(:,:)
    type(density_calc_t) :: dens_calc
    FLOAT, allocatable :: gradvh(:, :), heat_current(:, :, :)

    PUSH_SUB(output_hamiltonian)
   

    if (outp%what_now(OPTION__OUTPUT__POTENTIAL, iter)) then
      SAFE_ALLOCATE(v0(1:der%mesh%np, 1:hm%d%dim))
      v0(1:der%mesh%np, 1) = hm%ep%vpsl(1:der%mesh%np)
      call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, "v0", namespace, &
        space, der%mesh, v0(:, 1), units_out%energy, err, ions = ions, grp = grp)
      SAFE_DEALLOCATE_A(v0)

      if (hm%ep%classical_pot > 0) then
        call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, "vc", namespace, &
          space, der%mesh, hm%ep%Vclassical, units_out%energy, err, ions = ions, grp = grp)
      end if

      if (allocated(hm%v_static)) then
        call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, "vext", namespace, &
          space, der%mesh, hm%v_static, units_out%energy, err, ions = ions, grp = grp)
      end if

      if (hm%theory_level /= INDEPENDENT_PARTICLES) then
        call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vh', namespace, &
          space, der%mesh, hm%vhartree, units_out%energy, err, ions = ions, grp = grp)
        if (outp%what(OPTION__OUTPUT__POTENTIAL_GRADIENT)) then
          SAFE_ALLOCATE(vh(1:der%mesh%np_part))
          SAFE_ALLOCATE(gradvh(1:der%mesh%np, 1:space%dim))
          vh(1:der%mesh%np) = hm%vhartree(1:der%mesh%np)
          call dderivatives_grad(der, vh, gradvh(1:der%mesh%np, 1:space%dim))
          call io_function_output_vector(outp%how(OPTION__OUTPUT__POTENTIAL_GRADIENT), dir, 'grad_vh', namespace, &
            space, der%mesh, gradvh(:, :), units_out%force, err, ions = ions, grp = grp)

          SAFE_ALLOCATE(v0(1:der%mesh%np_part, 1))
          v0(1:der%mesh%np, 1) = hm%ep%vpsl(1:der%mesh%np)
          call dderivatives_grad(der, v0(1:der%mesh%np_part, 1), gradvh(1:der%mesh%np, 1:space%dim))
          call io_function_output_vector(outp%how(OPTION__OUTPUT__POTENTIAL_GRADIENT), dir, 'grad_v0', namespace, &
            space, der%mesh, gradvh(:, :), units_out%force, err, ions = ions, grp = grp)
          SAFE_DEALLOCATE_A(v0)
          SAFE_DEALLOCATE_A(vh)
          SAFE_DEALLOCATE_A(gradvh)
        end if
        
        SAFE_ALLOCATE(potential(1:der%mesh%np))
        do is = 1, hm%d%nspin
          if (hm%d%nspin == 1) then
            write(fname, '(a)') 'vxc'
          else
            write(fname, '(a,i1)') 'vxc-sp', is
          end if
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, fname, namespace, space, &
            der%mesh, hm%vxc(:, is), units_out%energy, err, ions = ions, grp = grp)
          
          ! finally the full KS potential (without non-local PP contributions)
          potential = hm%ep%vpsl + hm%vhxc(:, is)
          if(hm%d%nspin == 1) then
            write(fname, '(a)') 'vks'
          else
            write(fname, '(a,i1)') 'vks-sp', is
          end if
          if (hm%ep%classical_pot > 0) then
            call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, fname, namespace, space, &
              der%mesh, potential + hm%ep%Vclassical, units_out%energy, err, ions = ions, grp = grp)
          else
            call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, fname, namespace, space, &
              der%mesh, potential, units_out%energy, err, ions = ions, grp = grp)
          end if
        end do
        SAFE_DEALLOCATE_A(potential)
      end if

      !PCM potentials
      if (hm%theory_level == KOHN_SHAM_DFT .and. hm%pcm%run_pcm) then
        if (hm%pcm%solute .and. hm%pcm%localf) then
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vpcm', namespace, space, &
            der%mesh, hm%pcm%v_e_rs + hm%pcm%v_n_rs + hm%pcm%v_ext_rs , & 
            units_out%energy, err, ions = ions, grp = grp)
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vpcm_sol', namespace, space, &
            der%mesh, hm%pcm%v_e_rs + hm%pcm%v_n_rs , & 
            units_out%energy, err, ions = ions, grp = grp)
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vpcm_e', namespace, space, &
            der%mesh, hm%pcm%v_e_rs , & 
            units_out%energy, err, ions = ions, grp = grp)
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vpcm_n', namespace, space, &
            der%mesh, hm%pcm%v_n_rs , & 
            units_out%energy, err, ions = ions, grp = grp)
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vpcm_ext', namespace, space, &
            der%mesh, hm%pcm%v_ext_rs , & 
            units_out%energy, err, ions = ions, grp = grp)
        else if (hm%pcm%solute .and. .not.hm%pcm%localf) then
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vpcm_sol', namespace, space, &
            der%mesh, hm%pcm%v_e_rs + hm%pcm%v_n_rs , & 
            units_out%energy, err, ions = ions, grp = grp)
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vpcm_e', namespace, space, &
            der%mesh, hm%pcm%v_e_rs , & 
            units_out%energy, err, ions = ions, grp = grp)
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vpcm_n', namespace, space, &
            der%mesh, hm%pcm%v_n_rs , & 
            units_out%energy, err, ions = ions, grp = grp)
        else if (.not.hm%pcm%solute .and. hm%pcm%localf) then
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'vpcm_ext', namespace, space, &
            der%mesh, hm%pcm%v_ext_rs , & 
            units_out%energy, err, ions = ions, grp = grp)
        end if
      end if

      if(hm%self_induced_magnetic) then
        ! unit of magnetic field is same as of electric field, and same as force (since e = 1)
        select case(space%dim)
        case(3)
          do idir = 1, space%dim
            call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'Bind_'//index2axis(idir), namespace, &
              space, der%mesh, hm%b_ind(:, idir), &
              units_out%force, err, ions = ions, grp = grp)
          end do
        case(2)
          call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), dir, 'Bind_z', namespace, &
            space, der%mesh, hm%b_ind(:, 1), units_out%force, err, ions = ions, grp = grp)
        end select
      end if
    end if


    if (outp%what_now(OPTION__OUTPUT__XC_DENSITY, iter) .and. hm%theory_level /= INDEPENDENT_PARTICLES) then
      SAFE_ALLOCATE(v0(1:der%mesh%np_part, 1))
      SAFE_ALLOCATE(nxc(1:der%mesh%np))

      do is = 1, hm%d%nspin
        if(hm%d%nspin == 1) then
          write(fname, '(a)') 'nxc'
        else
          write(fname, '(a,i1)') 'nxc-sp', is
        end if
                
        v0(1:der%mesh%np, 1) = hm%vxc(1:der%mesh%np, is)

        call dderivatives_lapl(der, v0(:, 1), nxc)

        call dio_function_output(outp%how(OPTION__OUTPUT__XC_DENSITY), dir, fname, namespace, &
          space, der%mesh, nxc, units_out%energy, err, ions = ions, grp = grp)
        
      end do

      SAFE_DEALLOCATE_A(v0)
      SAFE_DEALLOCATE_A(nxc)
    end if

    if (outp%what_now(OPTION__OUTPUT__CURRENT, iter)) then
      
      if(states_are_complex(st)) then
        ASSERT(allocated(st%current))

        do is = 1, hm%d%nspin

          if(st%d%nspin == 1) then
            write(fname, '(2a)') 'current'
          else
            write(fname, '(a,i1)') 'current-sp', is
          end if
          
          call io_function_output_vector(outp%how(OPTION__OUTPUT__CURRENT), dir, fname, namespace, space, der%mesh, &
            st%current(:, :, is), (unit_one/units_out%time)*units_out%length**(1 - space%dim), err, &
            ions = ions, grp = st%dom_st_kpt_mpi_grp)

        end do
      else
        message(1) = 'No current density output for real states since it is identically zero.'
        call messages_warning(1, namespace=namespace)
      endif
    end if
   
    if (outp%what_now(OPTION__OUTPUT__CURRENT_KPT, iter)) then

      if(states_are_complex(st)) then
      
        SAFE_ALLOCATE(current_kpt(st%d%kpt%start:st%d%kpt%end, 1:space%dim)) 
        do ik = st%d%kpt%start,st%d%kpt%end
          do idir = 1, space%dim
            current_kpt(ik,idir) = dmf_integrate(der%mesh, st%current_kpt(:, idir, ik), reduce = .false.)
          end do
        end do
        if(der%mesh%parallel_in_domains) then
          call der%mesh%allreduce(current_kpt, dim = (/st%d%kpt%end-st%d%kpt%start+1, space%dim/))
        end if

        write(fname, '(2a)') 'current_kpt'
        call io_function_output_vector_BZ(outp%how(OPTION__OUTPUT__CURRENT_KPT), dir, fname, namespace, space, &
          st%d%kpt, hm%kpoints, current_kpt(:, :), (unit_one/units_out%time)*units_out%length**(1 - space%dim), err, &
          grp = st%st_kpt_mpi_grp)
        SAFE_DEALLOCATE_A(current_kpt)
      else
        message(1) = 'No current density output for real states since it is identically zero.'
        call messages_warning(1, namespace=namespace)
      endif
    end if

    if (outp%what_now(OPTION__OUTPUT__HEAT_CURRENT, iter)) then
      
      if(states_are_complex(st)) then

        SAFE_ALLOCATE(heat_current(1:der%mesh%np_part, 1:space%dim, 1:st%d%nspin))
      
        call current_heat_calculate(space, der, hm, st, heat_current)
      
        do is = 1, hm%d%nspin
          if(st%d%nspin == 1) then
            write(fname, '(2a)') 'heat_current'
          else
            write(fname, '(a,i1)') 'heat_current-sp', is
          end if
        
          call io_function_output_vector(outp%how(OPTION__OUTPUT__HEAT_CURRENT), dir, fname, namespace, space, der%mesh, &
            st%current(:, :, is), (unit_one/units_out%time)*units_out%length**(1 - space%dim), err, &
            ions = ions, grp = st%dom_st_kpt_mpi_grp)
        
          SAFE_DEALLOCATE_A(heat_current)
        end do
      else
        message(1) = 'No current density output for real states since it is identically zero.'
        call messages_warning(1)
      endif
    end if
    
    if (outp%what_now(OPTION__OUTPUT__DENSITY_KPT, iter)) then
      SAFE_ALLOCATE(density_kpt(1:st%d%nik, 1:st%d%nspin))
      density_kpt(1:st%d%nik, 1:st%d%nspin) = M_ZERO

      SAFE_ALLOCATE(density_tmp(1:gr%mesh%np, st%d%nspin))

      !Compute the k-resolved density and integrate it over the mesh
      do ik = st%d%kpt%start,st%d%kpt%end
        call density_calc_init(dens_calc, st, gr, density_tmp)
        do ib = st%group%block_start, st%group%block_end
          call density_calc_accumulate(dens_calc, st%group%psib(ib, ik))
        end do
        call density_calc_end(dens_calc, allreduce=.false., symmetrize=.false.)
 
        do is = 1, st%d%nspin
          density_kpt(ik, is) = density_kpt(ik, is) + dmf_integrate(der%mesh, density_tmp(:,is), reduce = .false.)
        end do
      end do
    
      if(st%parallel_in_states .or. st%d%kpt%parallel) then
        call comm_allreduce(st%dom_st_kpt_mpi_grp, density_kpt)
      end if

      do is = 1, st%d%nspin
        if(st%d%nspin == 1) then
          write(fname, '(2a)') 'density_kpt'
        else
          write(fname, '(a,i1)') 'density_kpt-sp', is
        end if
        call io_function_output_global_BZ(outp%how(OPTION__OUTPUT__DENSITY_KPT), dir, fname, namespace, &
          hm%kpoints, density_kpt(:, is), unit_one, err)
      end do
      SAFE_DEALLOCATE_A(density_tmp)
      SAFE_DEALLOCATE_A(density_kpt)
    end if
 
    POP_SUB(output_hamiltonian)
  end subroutine output_hamiltonian


  ! ---------------------------------------------------------
  subroutine output_scalar_pot(outp, namespace, space, dir, gr, ions, hm, time)
    type(output_t),           intent(in)    :: outp
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=*),         intent(in)    :: dir
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(in)    :: ions
    type(hamiltonian_elec_t), intent(inout) :: hm
    FLOAT, optional,          intent(in)    :: time

    integer :: is, err
    character(len=80) :: fname
    FLOAT, allocatable :: scalar_pot(:)

    PUSH_SUB(output_scalar_pot)

    if (outp%what(OPTION__OUTPUT__EXTERNAL_TD_POTENTIAL)) then
      SAFE_ALLOCATE(scalar_pot(1:gr%mesh%np))
      do is = 1, hm%ext_lasers%no_lasers
        write(fname, '(a,i1)') 'scalar_pot-', is
        scalar_pot = M_ZERO
        call laser_potential(hm%ext_lasers%lasers(is), gr%mesh, scalar_pot, time=time)
        call dio_function_output(outp%how(OPTION__OUTPUT__EXTERNAL_TD_POTENTIAL), dir, fname, namespace, &
          space, gr%mesh, scalar_pot, units_out%energy, err, ions = ions)
      end do
      SAFE_DEALLOCATE_A(scalar_pot)
    end if

    POP_SUB(output_scalar_pot)
  end subroutine output_scalar_pot


  ! ---------------------------------------------------------
  subroutine output_kick(outp, namespace, space, dir, mesh, ions, kick)
    type(output_t),    intent(in) :: outp
    type(namespace_t), intent(in) :: namespace
    type(space_t),     intent(in) :: space
    character(len=*),  intent(in) :: dir
    type(mesh_t),      intent(in) :: mesh
    type(ions_t),      intent(in) :: ions
    type(kick_t),      intent(in) :: kick

    integer :: err
    CMPLX, allocatable :: kick_function(:)
    
    PUSH_SUB(output_kick)

    if (outp%what(OPTION__OUTPUT__DELTA_PERTURBATION)) then
      SAFE_ALLOCATE(kick_function(1:mesh%np))
      call kick_function_get(space, mesh, kick, kick_function, 1)
      call zio_function_output(outp%how(OPTION__OUTPUT__DELTA_PERTURBATION), dir, "kick_function", namespace, &
        space, mesh, kick_function(:), units_out%energy, err, ions = ions)
      SAFE_DEALLOCATE_A(kick_function)
    end if
    
    POP_SUB(output_kick)
  end subroutine output_kick

  ! ---------------------------------------------------------
  subroutine output_xc_torque(outp, namespace, dir, mesh, hm, st, ions, space)
    type(output_t),           intent(in) :: outp
    type(namespace_t),        intent(in) :: namespace
    character(len=*),         intent(in) :: dir
    type(mesh_t),             intent(in) :: mesh
    type(hamiltonian_elec_t), intent(in) :: hm
    type(states_elec_t),      intent(in) :: st
    type(ions_t),             intent(in) :: ions
    type(space_t),            intent(in) :: space
 

    FLOAT, allocatable :: torque(:,:)
    type(unit_t) :: fn_unit
    integer :: err

    PUSH_SUB(output_xc_torque)

    if (outp%what(OPTION__OUTPUT__XC_TORQUE)) then
      SAFE_ALLOCATE(torque(1:mesh%np, 1:3))

      call calc_xc_torque(mesh, hm%vxc, st, torque)

      fn_unit = units_out%length**(1 - 2*space%dim)
      call io_function_output_vector(outp%how(OPTION__OUTPUT__XC_TORQUE), dir, 'xc_torque', namespace, space, mesh, &
        torque, fn_unit, err, ions = ions, grp = st%dom_st_kpt_mpi_grp)

      SAFE_DEALLOCATE_A(torque)
    end if

    POP_SUB(output_xc_torque)
  end subroutine output_xc_torque


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
