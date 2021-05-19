!! Copyright (C) 2019 F. Bonafe, R. Jestaedt, H. Appel
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
  subroutine output_mxll_init(outp, namespace, space)
    type(output_t),            intent(out) :: outp
    type(namespace_t),         intent(in)  :: namespace
    type(space_t),             intent(in)  :: space

    PUSH_SUB(output_mxll_init)

    !%Variable MaxwellOutput
    !%Type block
    !%Default none
    !%Section Output
    !%Description
    !% Specifies what to print. The output files are written at the end of the run into the output directory for the
    !% Maxwell run.
    !% Time-dependent simulations print only per iteration, including always the last. The frequency of output per iteration
    !% is set by <tt>OutputInterval</tt> and the directory is set by <tt>OutputIterDir</tt>.
    !% Each option must be in a separate row. Optionally individual output formats and output intervals can be defined
    !% for each row or they can be read separately from <tt>OutputFormat</tt> and <tt>MaxwellOutputInterval</tt> variables
    !% in the input file.
    !%
    !% Example:
    !% <br><br><tt>%MaxwellOutput
    !% <br>&nbsp;&nbsp;electric_field
    !% <br>&nbsp;&nbsp;magnetic_field
    !% <br>%<br></tt>
    !% This block supports all the formats of the <tt>Output</tt> block.
    !% See <tt>Output</tt>.
    !%Option electric_field 1
    !% Output of the electric field
    !%Option magnetic_field 2
    !% Output of the magnetic field
    !%Option trans_electric_field 3
    !% Output of the transversal electric field
    !%Option trans_magnetic_field 4
    !% Output of the transversal magnetic field
    !%Option long_electric_field 5
    !% Output of the longitudinal electric field
    !%Option long_magnetic_field 6
    !% Output of the longitudinal magnetic field
    !%Option div_electric_field 7
    !% Output of the divergence of the electric field
    !%Option div_magnetic_field 8
    !% Output of the divergence of the magnetic field
    !%Option maxwell_vector_potential 9
    !% Output of the vector potential
    !%Option poynting_vector 10
    !% Output of the Maxwell Poynting vector
    !%Option maxwell_energy_density 11
    !% Output of the electromagnetic density
    !%Option maxwell_current 12
    !% Output of the Maxwell current i.e. matter current plus external current.
    !%Option external_current 13
    !% Output of the external Maxwell current
    !%Option electric_dipole_potential 14
    !% Output of the electric dipole potential
    !%Option electric_quadrupole_potential 15
    !% Output of the electric quadrupole potential
    !%Option charge_density 16
    !% Output of the charge density calculated by the divergence of the electric field.
    !%Option charge_density_diff 17
    !% Output of the charge density difference, one calculated via matter calculation, the other
    !% one calculated with Div(E) = rho/epsilon.
    !%Option maxwell_pol_density 18
    !% Output of the polarization density on the Maxwell grid.
    !%Option medium_variables_electric 19
    !% Output of the electric displacement field, the polarization density and the electric susceptibility.
    !%Option medium_variables_magnetic 20
    !% Output of the H-magnetic field, the magnetization and the magnetic susceptibility.
    !%Option test_output 21
    !% Output of a test function for debugging
    !%End

    !%Variable MaxwellOutputInterval
    !%Type integer
    !%Default 50
    !%Section Output
    !%Description
    !% The output requested by variable <tt>MaxwellOutput</tt> is written
    !% to the directory <tt>MaxwellOutputIterDir</tt>
    !% when the iteration number is a multiple of the <tt>MaxwellOutputInterval</tt> variable.
    !% Subdirectories are named Y.X, where Y is <tt>td</tt>, <tt>scf</tt>, or <tt>unocc</tt>, and
    !% X is the iteration number. To use the working directory, specify <tt>"."</tt>
    !% (Output of restart files is instead controlled by <tt>MaxwellRestartWriteInterval</tt>.)
    !% Must be >= 0. If it is 0, then no output is written. 
    !% This variable can also be defined inside the <tt>MaxwellOutput</tt> block.
    !% See <tt>MaxwellOutput</tt>.
    !%End

    outp%what = .false.
    call io_function_read_what_how_when(namespace, space, outp%what, outp%how, outp%output_interval, &
      'MaxwellOutput', 'OutputFormat', 'MaxwellOutputInterval')



    !%Variable MaxwellOutputIterDir
    !%Default "output_iter"
    !%Type string
    !%Section Output
    !%Description
    !% The name of the directory where <tt>Octopus</tt> stores information
    !% such as the density, forces, etc. requested by variable <tt>MaxwellOutput</tt>
    !% in the format specified by <tt>OutputHow</tt>.
    !% This information is written while iterating <tt>CalculationMode = maxwell</tt>
    !% according to <tt>OutputInterval</tt>, and has nothing to do with the restart information.
    !%End
    call parse_variable(namespace, 'MaxwellOutputIterDir', "output_iter", outp%iter_dir)
    if (any(outp%what) .and. maxval(outp%output_interval) > 0) then
      call io_mkdir(outp%iter_dir, namespace)
    end if
    call add_last_slash(outp%iter_dir)

    !%Variable MaxwellRestartWriteInterval
    !%Type integer
    !%Default 50
    !%Section Execution::IO
    !%Description
    !% Restart data is written when the iteration number is a multiple of the
    !% <tt>MaxwellRestartWriteInterval</tt> variable. (Other output is controlled by <tt>MaxwellOutputInterval</tt>.)
    !%End
    call parse_variable(namespace, 'MaxwellRestartWriteInterval', 50, outp%restart_write_interval)
    if (outp%restart_write_interval <= 0) then
      message(1) = "MaxwellRestartWriteInterval must be > 0."
      call messages_fatal(1, namespace=namespace)
    end if

    if (outp%what(OPTION__MAXWELLOUTPUT__ELECTRIC_FIELD)) then
      outp%wfs_list = trim("1-3")
    end if

    if (outp%what(OPTION__MAXWELLOUTPUT__MAGNETIC_FIELD)) then
      outp%wfs_list = trim("1-3")
    end if

    if (outp%what(OPTION__MAXWELLOUTPUT__TRANS_ELECTRIC_FIELD)) then
      outp%wfs_list = trim("1-3")
    end if

    if (outp%what(OPTION__MAXWELLOUTPUT__TRANS_MAGNETIC_FIELD)) then
      outp%wfs_list = trim("1-3")
    end if

    if (outp%what(OPTION__MAXWELLOUTPUT__MAXWELL_CURRENT)) then
      outp%wfs_list = trim("1-3")
    end if

    POP_SUB(output_mxll_init)
  end subroutine output_mxll_init


  ! ---------------------------------------------------------
  subroutine output_mxll(outp, namespace, space, gr_mxll, st_mxll, hm_mxll, time, ions, dir, gr_elec, st_elec, hm_elec)
    type(output_t),                     intent(in)    :: outp
    type(namespace_t),                  intent(in)    :: namespace
    type(space_t),                      intent(in)    :: space
    type(states_mxll_t),                intent(inout) :: st_mxll
    type(hamiltonian_mxll_t),           intent(in)    :: hm_mxll
    type(grid_t),                       intent(inout) :: gr_mxll
    FLOAT,                              intent(in)    :: time
    type(ions_t),                       intent(inout) :: ions
    character(len=*),                   intent(in)    :: dir
    type(grid_t),             optional, intent(inout) :: gr_elec
    type(states_elec_t),      optional, intent(inout) :: st_elec
    type(hamiltonian_elec_t), optional, intent(in)    :: hm_elec

    PUSH_SUB(output_mxll)

    if (any(outp%what)) then
      message(1) = "Info: Writing output to " // trim(dir)
      call messages_info(1)
      call io_mkdir(dir, namespace)
    endif

    call output_states_mxll(outp, namespace, space, dir, st_mxll, gr_mxll%mesh, hm_mxll, ions)
    call output_energy_density_mxll(outp, namespace, space, dir, hm_mxll, gr_mxll%mesh, ions)
    call output_poynting_vector(outp, namespace, space, dir, st_mxll, gr_mxll, hm_mxll, ions)
    call output_transverse_rs_state(outp, namespace, space, dir, st_mxll, gr_mxll%mesh, ions)
    call output_longitudinal_rs_state(outp, namespace, space, dir, st_mxll, gr_mxll%mesh, ions)
    call output_divergence_rs_state(outp, namespace, space, dir, st_mxll, gr_mxll, ions)
    call output_vector_potential(outp, namespace, space, dir, st_mxll, gr_mxll%mesh, hm_mxll, ions)
    call output_external_current_density(outp, namespace, space, dir, st_mxll, gr_mxll%mesh, hm_mxll, ions, time)
    call output_charge_density_mxll(outp, namespace, space, dir, st_mxll, gr_mxll, hm_mxll, ions)

    if (present(hm_elec) .and. present(gr_elec) .and. present(st_elec)) then
      call output_coupling_potentials(outp, namespace, dir, hm_elec, gr_elec, ions)
      call output_current_density(outp, namespace, dir, st_mxll, gr_mxll, hm_mxll, st_elec, gr_elec, hm_elec, ions, time)
      call output_medium_variables_electric(outp, namespace, space, dir, st_mxll, gr_mxll%mesh, hm_mxll, ions, st_elec, &
        gr_elec, hm_elec)
      call output_medium_variables_magnetic(outp, namespace, space, dir, st_mxll, gr_mxll%mesh, hm_mxll, ions, st_elec, &
        gr_elec, hm_elec)
    end if

    POP_SUB(output_mxll)
  end subroutine output_mxll


  ! ---------------------------------------------------------
  subroutine output_states_mxll(outp, namespace, space, dir, st, mesh, hm, ions)
    type(output_t),           intent(in)    :: outp
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=*),         intent(in)    :: dir
    type(states_mxll_t),      intent(inout) :: st
    type(mesh_t),             intent(in)    :: mesh
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(ions_t),             intent(in)    :: ions

    integer :: idim, ierr
    character(len=MAX_PATH_LEN) :: fname
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp(:,:)

    PUSH_SUB(output_states_mxll)

    ! Electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__ELECTRIC_FIELD)) then
      fn_unit = units_out%energy/units_out%length
      SAFE_ALLOCATE(dtmp(1:mesh%np, 1:st%dim))
      call get_electric_field_state(st%rs_state, mesh, dtmp, st%ep, mesh%np)
      do idim = 1, st%dim
        write(fname, '(2a)') 'e_field-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__ELECTRIC_FIELD), dir, fname, namespace, space, &
          mesh, dtmp(:, idim), fn_unit, ierr, ions = ions)
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    ! Magnetic field
    if (outp%what(OPTION__MAXWELLOUTPUT__MAGNETIC_FIELD)) then
      fn_unit = unit_one/units_out%length**2
      SAFE_ALLOCATE(dtmp(1:mesh%np, 1:st%dim))
      call get_magnetic_field_state(st%rs_state, mesh, st%rs_sign, dtmp, st%mu, mesh%np)
      do idim = 1, st%dim
        write(fname, '(2a)') 'b_field-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__MAGNETIC_FIELD), dir, fname, namespace, space, &
          mesh, dtmp(:, idim), fn_unit, ierr, ions = ions)
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    POP_SUB(output_states_mxll)
  end subroutine output_states_mxll


  !----------------------------------------------------------
  subroutine output_energy_density_mxll(outp, namespace, space, dir, hm, mesh, ions)   !< have to set unit output correctly
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(mesh_t),              intent(in)    :: mesh
    type(ions_t),              intent(in)    :: ions
    type(output_t),            intent(in)    :: outp

    integer :: ierr

    PUSH_SUB(output_energy_density_mxll)

    ! Maxell energy density
    if (outp%what(OPTION__MAXWELLOUTPUT__MAXWELL_ENERGY_DENSITY)) then
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__MAXWELL_ENERGY_DENSITY), dir, "maxwell_energy_density", &
        namespace, space, mesh, hm%energy%energy_density(:), units_out%energy/units_out%length**3, ierr, ions = ions)
    end if

    POP_SUB(output_energy_density_mxll)
  end subroutine output_energy_density_mxll


  !----------------------------------------------------------
  subroutine output_poynting_vector(outp, namespace, space, dir, st, gr, hm, ions)
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(inout) :: st
    type(grid_t),              intent(inout) :: gr
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    type(ions_t),              intent(in)    :: ions

    integer :: idim, ierr
    character(len=MAX_PATH_LEN) :: fname
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp(:,:)

    PUSH_SUB(output_poynting_vector)

    ! Maxwell Poynting vector
    if (outp%what(OPTION__MAXWELLOUTPUT__POYNTING_VECTOR)) then
      fn_unit = units_out%energy/(units_out%time*units_out%length**2)
      SAFE_ALLOCATE(dtmp(1:gr%mesh%np, 1:st%dim))
      call get_poynting_vector(gr, st, st%rs_state, st%rs_sign, dtmp, st%ep, st%mu)
      do idim = 1, st%dim
        write(fname, '(2a)') 'poynting_vector-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__POYNTING_VECTOR), dir, fname, namespace, space, &
          gr%mesh, dtmp(:, idim), fn_unit, ierr, ions = ions)
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    POP_SUB(output_poynting_vector)
  end subroutine output_poynting_vector


  !----------------------------------------------------------
  subroutine output_transverse_rs_state(outp, namespace, space, dir, st, mesh, ions)    !< have to set unit output correctly
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(in)    :: st
    type(mesh_t),              intent(in)    :: mesh
    type(ions_t),              intent(in)    :: ions

    character(len=MAX_PATH_LEN) :: fname
    integer :: idim, ierr
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp(:,:)

    PUSH_SUB(output_transverse_rs_state)

    ! transverse component of the electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__TRANS_ELECTRIC_FIELD)) then
      fn_unit = units_out%energy/units_out%length
      SAFE_ALLOCATE(dtmp(1:mesh%np, 1:st%dim))
      call get_electric_field_state(st%rs_state_trans, mesh, dtmp, st%ep(1:mesh%np), mesh%np)
      do idim=1, st%dim
        write(fname, '(2a)') 'e_field_trans-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__TRANS_ELECTRIC_FIELD), dir, fname, namespace, space, &
          mesh, dtmp(:, idim), fn_unit, ierr, ions = ions)
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    ! transverse component of the magnetic field
    if (outp%what(OPTION__MAXWELLOUTPUT__TRANS_MAGNETIC_FIELD)) then
      fn_unit = unit_one/units_out%length**2
      SAFE_ALLOCATE(dtmp(1:mesh%np, 1:st%dim))
      call get_magnetic_field_state(st%rs_state_trans, mesh, st%rs_sign, dtmp, st%ep(1:mesh%np), mesh%np)
      do idim=1, st%dim
        write(fname, '(2a)') 'b_field_trans-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__TRANS_MAGNETIC_FIELD), dir, fname, namespace, space, &
          mesh, dtmp(:, idim), fn_unit, ierr, ions = ions)
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    POP_SUB(output_transverse_rs_state)
  end subroutine output_transverse_rs_state


  !----------------------------------------------------------
  subroutine output_longitudinal_rs_state(outp, namespace, space, dir, st, mesh, ions)    !< have to set unit output correctly
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(in)    :: st
    type(mesh_t),              intent(in)    :: mesh
    type(ions_t),              intent(in)    :: ions

    character(len=MAX_PATH_LEN) :: fname
    integer :: idim, ierr
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp(:,:)

    PUSH_SUB(output_longitudinal_rs_state)

    ! longitudinal component of the electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__LONG_ELECTRIC_FIELD)) then
      fn_unit = units_out%energy/units_out%length
      SAFE_ALLOCATE(dtmp(1:mesh%np, 1:st%dim))
      call get_electric_field_state(st%rs_state_long, mesh, dtmp, st%ep(1:mesh%np), mesh%np)
      do idim = 1, st%dim
        write(fname, '(2a)') 'e_field_long-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__LONG_ELECTRIC_FIELD), dir, fname, namespace, space, &
          mesh, dtmp(:, idim), fn_unit, ierr, ions = ions)
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    ! longitudinal component of the magnetic field
    if (outp%what(OPTION__MAXWELLOUTPUT__LONG_MAGNETIC_FIELD)) then
      fn_unit = unit_one/units_out%length**2
      SAFE_ALLOCATE(dtmp(1:mesh%np, 1:st%dim))
      call get_magnetic_field_state(st%rs_state_long, mesh, st%rs_sign, dtmp, st%mu(1:mesh%np), mesh%np)
      do idim = 1, st%dim
        write(fname, '(2a)') 'b_field_long-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__LONG_MAGNETIC_FIELD), dir, fname, namespace, space, &
          mesh, dtmp(:, idim), fn_unit, ierr, ions = ions)
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    POP_SUB(output_longitudinal_rs_state)
  end subroutine output_longitudinal_rs_state


  !----------------------------------------------------------
  subroutine output_divergence_rs_state(outp, namespace, space, dir, st, gr, ions)    !< have to set unit output correctly
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(in)    :: st
    type(grid_t),              intent(in)    :: gr
    type(ions_t),              intent(in)    :: ions

    integer :: ierr
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp_1(:,:), dtmp_2(:)

    PUSH_SUB(output_divergence_rs_state)

    ! divergence of the electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__DIV_ELECTRIC_FIELD)) then
      fn_unit = units_out%length**(-space%dim)
      SAFE_ALLOCATE(dtmp_1(1:gr%mesh%np_part, 1:st%dim))
      SAFE_ALLOCATE(dtmp_2(1:gr%mesh%np))
      dtmp_1 = M_ZERO
      call get_electric_field_state(st%rs_state, gr%mesh, dtmp_1(1:gr%mesh%np,:), st%mu(1:gr%mesh%np), gr%mesh%np)
      call get_divergence_field(gr, dtmp_1, dtmp_2, .false.)
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__DIV_ELECTRIC_FIELD), dir, "e_field_div", namespace, space, &
        gr%mesh, dtmp_2, fn_unit, ierr, ions = ions)
      SAFE_DEALLOCATE_A(dtmp_1)
      SAFE_DEALLOCATE_A(dtmp_2)
    end if

    ! divergence of the magnetic field
    if (outp%what(OPTION__MAXWELLOUTPUT__DIV_MAGNETIC_FIELD)) then
      ! unit does not matter, should be zero
      fn_unit = unit_one
      SAFE_ALLOCATE(dtmp_1(1:gr%mesh%np_part, 1:st%dim))
      SAFE_ALLOCATE(dtmp_2(1:gr%mesh%np))
      dtmp_1 = M_ZERO
      call get_magnetic_field_state(st%rs_state, gr%mesh, st%rs_sign, dtmp_1(1:gr%mesh%np,:), st%mu(1:gr%mesh%np), gr%mesh%np)
      call get_divergence_field(gr, dtmp_1, dtmp_2, .false.)
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__DIV_MAGNETIC_FIELD), dir, "b_field_div", namespace, space, &
        gr%mesh, dtmp_2, fn_unit, ierr, ions = ions)
      SAFE_DEALLOCATE_A(dtmp_1)
      SAFE_DEALLOCATE_A(dtmp_2)
    end if

    POP_SUB(output_divergence_rs_state)
  end subroutine output_divergence_rs_state


  !----------------------------------------------------------
  subroutine output_vector_potential(outp, namespace, space, dir, st, mesh, hm, ions)     !< have to set unit output correctly
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(in)    :: st
    type(mesh_t),              intent(in)    :: mesh
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    type(ions_t),              intent(in)    :: ions

    character(len=MAX_PATH_LEN) :: fname
    integer :: idim, ierr
    type(unit_t) :: fn_unit

    PUSH_SUB(output_vector_potential)

    ! Maxwell vector potential
    if (outp%what(OPTION__MAXWELLOUTPUT__MAXWELL_VECTOR_POTENTIAL)) then
      fn_unit = unit_one
      do idim = 1, space%dim
        write(fname, '(2a)') 'vector_potential-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__MAXWELL_VECTOR_POTENTIAL),&
          dir, fname, namespace, space, mesh, hm%vector_potential(1:mesh%np, idim), fn_unit, &
          ierr, ions = ions)
      end do
    end if

    POP_SUB(output_vector_potential)
  end subroutine output_vector_potential


  !------------------------------------------------------------
  subroutine output_coupling_potentials(outp, namespace, dir, hm, gr, ions)    !< have to set unit output correctly
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    character(len=*),          intent(in)    :: dir
    type(hamiltonian_elec_t),  intent(in)    :: hm
    type(grid_t),              intent(in)    :: gr
    type(ions_t),              intent(in)    :: ions

    message(1) = "Maxwell-matter coupling potentials not implemented yet."
    call messages_fatal(1, namespace=namespace)

  end subroutine output_coupling_potentials


  !----------------------------------------------------------
  subroutine output_current_density(outp, namespace, dir, st_mxll, gr_mxll, hm_mxll, st_elec, gr_elec, hm_elec, ions, time)
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(inout) :: st_mxll
    type(grid_t),              intent(inout) :: gr_mxll
    type(hamiltonian_mxll_t),  intent(in)    :: hm_mxll
    type(states_elec_t),       intent(inout) :: st_elec
    type(grid_t),              intent(inout) :: gr_elec
    type(hamiltonian_elec_t),  intent(in)    :: hm_elec
    type(ions_t),              intent(in)    :: ions
    FLOAT,                     intent(in)    :: time

    message(1) = "Current density can not yet be calculated in a Maxwell propagation."
    call messages_fatal(1, namespace=namespace)

  end subroutine output_current_density


  !----------------------------------------------------------
  subroutine output_external_current_density(outp, namespace, space, dir, st, mesh, hm, ions, time)
    type(output_t),            intent(in)    :: outp
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    character(len=*),          intent(in)    :: dir
    type(states_mxll_t),       intent(inout) :: st
    type(mesh_t),              intent(in)    :: mesh
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    type(ions_t),              intent(in)    :: ions
    FLOAT,                     intent(in)    :: time

    character(len=MAX_PATH_LEN) :: fname
    integer                     :: ierr, idim
    CMPLX, allocatable          :: ztmp(:,:)
    FLOAT, allocatable          :: dtmp(:,:)
    type(unit_t)                :: fn_unit

    PUSH_SUB(output_external_current_density)

     ! Maxwell current density
    if (outp%what(OPTION__MAXWELLOUTPUT__EXTERNAL_CURRENT)) then
      if (hm%current_density_ext_flag) then
        fn_unit = (unit_one/units_out%time)/(units_out%length**2)      !< test both if its the same
        SAFE_ALLOCATE(ztmp(1:mesh%np, 1:st%dim))
        SAFE_ALLOCATE(dtmp(1:mesh%np, 1:st%dim))
        call get_rs_density_ext(st, mesh, time, ztmp)
        call get_current_state(ztmp, dtmp, mesh, st%ep, mesh%np)
        do idim = 1, st%dim
          write(fname, '(2a)') 'external_current-', index2axis(idim)
          call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__EXTERNAL_CURRENT), dir, fname, namespace, space, &
            mesh, dtmp(:, idim), fn_unit, ierr, ions = ions)
        end do
        SAFE_DEALLOCATE_A(ztmp)
        SAFE_DEALLOCATE_A(dtmp)
      end if
    end if

    POP_SUB(output_external_current_density)
  end subroutine output_external_current_density


  !----------------------------------------------------------
  subroutine output_charge_density_mxll(outp, namespace, space, dir, st, gr, hm, ions)    !< have to set unit output correctly
    type(output_t),                intent(in)    :: outp
    type(namespace_t),             intent(in)    :: namespace
    type(space_t),                 intent(in)    :: space
    character(len=*),              intent(in)    :: dir
    type(states_mxll_t),           intent(in)    :: st
    type(grid_t),                  intent(in)    :: gr
    type(hamiltonian_mxll_t),      intent(in)    :: hm
    type(ions_t),                  intent(in)    :: ions

    integer :: ierr
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp_1(:,:), dtmp_2(:)

    PUSH_SUB(output_charge_density_mxll)

    ! charge density calculated by the divergence of the electric field
    if (outp%what(OPTION__MAXWELLOUTPUT__CHARGE_DENSITY)) then
      fn_unit = units_out%length**(-space%dim)
      SAFE_ALLOCATE(dtmp_1(1:gr%mesh%np_part,1:st%dim))
      SAFE_ALLOCATE(dtmp_2(1:gr%mesh%np))
      call get_electric_field_state(st%rs_state, gr%mesh, dtmp_1, st%ep, gr%mesh%np)
      call get_divergence_field(gr, dtmp_1, dtmp_2, .true.)
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__CHARGE_DENSITY), dir, "charge_density", namespace, space, &
        gr%mesh, dtmp_2(:), fn_unit, ierr, ions = ions)
      SAFE_DEALLOCATE_A(dtmp_1)
      SAFE_DEALLOCATE_A(dtmp_2)
    end if

    POP_SUB(output_charge_density_mxll)
  end subroutine output_charge_density_mxll


  !----------------------------------------------------------
  subroutine output_medium_variables_electric(outp, namespace, space, dir, st, mesh, hm, ions, st_elec, gr_elec, hm_elec)    !< have to set unit output correctly
    type(output_t),                     intent(in)    :: outp
    type(namespace_t),                  intent(in)    :: namespace
    type(space_t),                      intent(in)    :: space
    character(len=*),                   intent(in)    :: dir
    type(states_mxll_t),                intent(in)    :: st
    type(mesh_t),                       intent(in)    :: mesh
    type(hamiltonian_mxll_t),           intent(in)    :: hm
    type(ions_t),                       intent(in)    :: ions
    type(grid_t),             optional, intent(in)    :: gr_elec
    type(states_elec_t),      optional, intent(in)    :: st_elec
    type(hamiltonian_elec_t), optional, intent(in)    :: hm_elec

    integer                     :: idim, ip, ierr
    character(len=MAX_PATH_LEN) :: fname
    type(unit_t)                :: fn_unit
    FLOAT, allocatable          :: polarization_ma_gr(:,:), polarization_mx_gr(:,:), E_field(:,:), D_field(:,:), e_susceptibility(:)

    PUSH_SUB(output_medium_variables_electric)

    if (outp%what(OPTION__MAXWELLOUTPUT__MEDIUM_VARIABLES_ELECTRIC)) then
      fn_unit = units_out%length**(1 - space%dim)
      SAFE_ALLOCATE(polarization_mx_gr(1:mesh%np, 1:st%dim))
      SAFE_ALLOCATE(E_field(1:mesh%np, 1:st%dim))
      SAFE_ALLOCATE(D_field(1:mesh%np, 1:st%dim))
      SAFE_ALLOCATE(e_susceptibility(1:mesh%np))
      call get_electric_field_state(st%rs_state, mesh, E_field, st%ep, mesh%np)
      polarization_ma_gr = M_ZERO
      polarization_mx_gr = M_ZERO
      do idim = 1, st%dim
        do ip = 1, mesh%np
          D_field(ip, idim) = E_field(ip, idim)*P_ep - polarization_mx_gr(ip, idim)
        end do
        write(fname, '(2a)') 'd_field-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__MEDIUM_VARIABLES_ELECTRIC), dir, fname, namespace, space, &
          mesh, D_field(:, idim), fn_unit, ierr, ions = ions)
      end do
      do ip = 1, mesh%np
        e_susceptibility(ip) = sqrt(sum(D_field(ip, :)**2))/(P_ep*sqrt(sum(E_field(ip, :)**2))) - M_ONE
      end do
      write(fname, '(1a)') 'electric_susceptibility'
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__MEDIUM_VARIABLES_ELECTRIC), dir, fname, namespace, space, &
        mesh, e_susceptibility(:), fn_unit, ierr, ions = ions)
      SAFE_DEALLOCATE_A(polarization_mx_gr)
      SAFE_DEALLOCATE_A(E_field)
      SAFE_DEALLOCATE_A(D_field)
      SAFE_DEALLOCATE_A(e_susceptibility)
    end if

    POP_SUB(output_medium_variables_electric)
  end subroutine output_medium_variables_electric


  !----------------------------------------------------------
  subroutine output_medium_variables_magnetic(outp, namespace, space, dir, st, mesh, hm, ions, st_elec, gr_elec, hm_elec)    !< have to set unit output correctly
    type(output_t),                     intent(in)    :: outp
    type(namespace_t),                  intent(in)    :: namespace
    type(space_t),                      intent(in)    :: space
    character(len=*),                   intent(in)    :: dir
    type(states_mxll_t),                intent(in)    :: st
    type(mesh_t),                       intent(in)    :: mesh
    type(hamiltonian_mxll_t),           intent(in)    :: hm
    type(ions_t),                       intent(in)    :: ions
    type(grid_t),             optional, intent(inout) :: gr_elec
    type(states_elec_t),      optional, intent(inout) :: st_elec
    type(hamiltonian_elec_t), optional, intent(in)    :: hm_elec

    integer                     :: idim, ip, ierr
    character(len=MAX_PATH_LEN) :: fname
    type(unit_t)                :: fn_unit
    FLOAT, allocatable          :: magnetization_ma_gr(:,:), magnetization_mx_gr(:,:), B_field(:,:), H_field(:,:), b_susceptibility(:)

    PUSH_SUB(output_medium_variables_magnetic)

    if (outp%what(OPTION__MAXWELLOUTPUT__MEDIUM_VARIABLES_MAGNETIC)) then
      fn_unit = units_out%length**(1 - space%dim)
      SAFE_ALLOCATE(magnetization_mx_gr(1:mesh%np, 1:st%dim))
      SAFE_ALLOCATE(B_field(1:mesh%np, 1:st%dim))
      SAFE_ALLOCATE(H_field(1:mesh%np, 1:st%dim))
      SAFE_ALLOCATE(b_susceptibility(1:mesh%np))
      call get_magnetic_field_state(st%rs_state, mesh, st%rs_sign, B_field, st%mu, mesh%np)
      magnetization_ma_gr = M_ZERO
      magnetization_mx_gr = M_ZERO
      do idim = 1, st%dim
        do ip = 1, mesh%np
          H_field(ip, idim) = B_field(ip,idim)/P_mu - magnetization_mx_gr(ip, idim)
        end do
        write(fname, '(2a)') 'h_field-', index2axis(idim)
        call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__MEDIUM_VARIABLES_MAGNETIC), dir, fname, namespace, space, &
          mesh, H_field(:, idim), fn_unit, ierr, ions = ions)
      end do
      do ip = 1, mesh%np
        b_susceptibility(ip) = sqrt(sum(H_field(ip, :)**2))/(P_ep*sqrt(sum(B_field(ip,:)**2))) - M_ONE
      end do
      write(fname, '(1a)') 'magnetic_susceptibility'
      call dio_function_output(outp%how(OPTION__MAXWELLOUTPUT__MEDIUM_VARIABLES_MAGNETIC), dir, fname, namespace, space, &
        mesh, b_susceptibility(:), fn_unit, ierr, ions = ions)
      SAFE_DEALLOCATE_A(magnetization_mx_gr)
      SAFE_DEALLOCATE_A(B_field)
      SAFE_DEALLOCATE_A(H_field)
      SAFE_DEALLOCATE_A(b_susceptibility)
    end if

    POP_SUB(output_medium_variables_magnetic)
  end subroutine output_medium_variables_magnetic

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
