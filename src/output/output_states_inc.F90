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

! ---------------------------------------------------------
subroutine output_states(outp, namespace, space, dir, st, gr, ions, hm, iter)
  type(output_t),           intent(in) :: outp
  type(namespace_t),        intent(in) :: namespace
  type(space_t),            intent(in) :: space
  character(len=*),         intent(in) :: dir
  type(states_elec_t),      intent(in) :: st
  type(grid_t),             intent(in) :: gr
  type(ions_t),             intent(in) :: ions
  type(hamiltonian_elec_t), intent(in) :: hm
  integer,                  intent(in) :: iter

  integer :: ik, ist, idim, idir, is, ierr, ip
  character(len=MAX_PATH_LEN) :: fname
  type(unit_t) :: fn_unit
  FLOAT, allocatable :: dtmp(:), elf(:,:), polarization(:, :)
  CMPLX, allocatable :: ztmp(:)
  type(dos_t) :: dos

  PUSH_SUB(output_states)

  if (outp%what_now(OPTION__OUTPUT__DENSITY, iter)) then
    fn_unit = units_out%length**(-space%dim)
    do is = 1, st%d%nspin
      if (st%d%nspin == 1) then
        write(fname, '(a)') 'density'
      else
        write(fname, '(a,i1)') 'density-sp', is
      end if
      call dio_function_output(outp%how(OPTION__OUTPUT__DENSITY), dir, fname, namespace, space, gr%mesh, &
        st%rho(:, is), fn_unit, ierr, ions = ions, grp = st%dom_st_kpt_mpi_grp)
    end do
  end if

  if (outp%what_now(OPTION__OUTPUT__POL_DENSITY, iter)) then
    fn_unit = units_out%length**(1 - space%dim)
    SAFE_ALLOCATE(polarization(1:gr%mesh%np, 1:space%dim))

    do is = 1, st%d%nspin
      do idir = 1, space%dim
        do ip = 1, gr%mesh%np
          polarization(ip, idir) = st%rho(ip, is)*gr%mesh%x(ip, idir)
        end do
      end do

      if (st%d%nspin == 1) then
        write(fname, '(a)') 'dipole_density'
      else
        write(fname, '(a,i1)') 'dipole_density-sp', is
      end if
      call io_function_output_vector(outp%how(OPTION__OUTPUT__POL_DENSITY),&
        dir, fname, namespace, space, gr%mesh, polarization, fn_unit, ierr, &
        ions = ions, grp = st%dom_st_kpt_mpi_grp)
    end do

    SAFE_DEALLOCATE_A(polarization)
  end if

  if (outp%what_now(OPTION__OUTPUT__WFS, iter)) then
    fn_unit = sqrt(units_out%length**(-space%dim))

    if (states_are_real(st)) then
      SAFE_ALLOCATE(dtmp(1:gr%mesh%np))
    else
      SAFE_ALLOCATE(ztmp(1:gr%mesh%np))
    end if

    do ist = st%st_start, st%st_end
      if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = st%d%kpt%start, st%d%kpt%end
          do idim = 1, st%d%dim
            if(st%d%nik > 1) then
              if(st%d%dim > 1) then
                write(fname, '(a,i3.3,a,i4.4,a,i1)') 'wf-k', ik, '-st', ist, '-sp', idim
              else
                write(fname, '(a,i3.3,a,i4.4)')      'wf-k', ik, '-st', ist
              end if
            else
              if(st%d%dim > 1) then
                write(fname, '(a,i4.4,a,i1)')        'wf-st', ist, '-sp', idim
              else
                write(fname, '(a,i4.4)')             'wf-st', ist
              end if
            end if

            if(states_are_real(st)) then
              call states_elec_get_state(st, gr%mesh, idim, ist, ik, dtmp)
              call dio_function_output(outp%how(OPTION__OUTPUT__WFS), dir, fname, namespace, space, gr%mesh, dtmp, &
                fn_unit, ierr, ions = ions)
            else
              call states_elec_get_state(st, gr%mesh, idim, ist, ik, ztmp)
              call zio_function_output(outp%how(OPTION__OUTPUT__WFS), dir, fname, namespace, space, gr%mesh, ztmp, &
                fn_unit, ierr, ions = ions)
            end if
          end do
        end do
      end if
    end do

    SAFE_DEALLOCATE_A(dtmp)
    SAFE_DEALLOCATE_A(ztmp)
  end if

  if (outp%what_now(OPTION__OUTPUT__WFS_SQMOD, iter)) then
    fn_unit = units_out%length**(-space%dim)
    SAFE_ALLOCATE(dtmp(1:gr%mesh%np_part))
    if (states_are_complex(st)) then
      SAFE_ALLOCATE(ztmp(1:gr%mesh%np))
    end if
    do ist = st%st_start, st%st_end
      if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = st%d%kpt%start, st%d%kpt%end
          do idim = 1, st%d%dim
            if(st%d%nik > 1) then
              if(st%d%dim > 1) then
                write(fname, '(a,i3.3,a,i4.4,a,i1)') 'sqm-wf-k', ik, '-st', ist, '-sp', idim
              else
                write(fname, '(a,i3.3,a,i4.4)')      'sqm-wf-k', ik, '-st', ist
              end if
            else
              if(st%d%dim > 1) then
                write(fname, '(a,i4.4,a,i1)')        'sqm-wf-st', ist, '-sp', idim
              else
                write(fname, '(a,i4.4)')             'sqm-wf-st', ist
              end if
            end if

            if (states_are_real(st)) then
              call states_elec_get_state(st, gr%mesh, idim, ist, ik, dtmp)
              dtmp(1:gr%mesh%np) = abs(dtmp(1:gr%mesh%np))**2
            else
              call states_elec_get_state(st, gr%mesh, idim, ist, ik, ztmp)
              dtmp(1:gr%mesh%np) = abs(ztmp(1:gr%mesh%np))**2
            end if
            call dio_function_output(outp%how(OPTION__OUTPUT__WFS_SQMOD), dir, fname, namespace, space, gr%mesh, &
              dtmp, fn_unit, ierr, ions = ions)
          end do
        end do
      end if
    end do
    SAFE_DEALLOCATE_A(dtmp)
    SAFE_DEALLOCATE_A(ztmp)
  end if

  if (outp%what_now(OPTION__OUTPUT__KINETIC_ENERGY_DENSITY, iter)) then
    fn_unit = units_out%energy * units_out%length**(-space%dim)
    SAFE_ALLOCATE(elf(1:gr%mesh%np, 1:st%d%nspin))
    call states_elec_calc_quantities(gr%der, st, hm%kpoints, .false., kinetic_energy_density = elf)
    select case(st%d%ispin)
    case(UNPOLARIZED)
      write(fname, '(a)') 'tau'
      call dio_function_output(outp%how(OPTION__OUTPUT__KINETIC_ENERGY_DENSITY), dir, trim(fname), namespace, space, &
        gr%mesh, elf(:,1), fn_unit, ierr, ions = ions, grp = st%dom_st_kpt_mpi_grp)
    case(SPIN_POLARIZED, SPINORS)
      do is = 1, st%d%nspin
        write(fname, '(a,i1)') 'tau-sp', is
        call dio_function_output(outp%how(OPTION__OUTPUT__KINETIC_ENERGY_DENSITY), dir, trim(fname), namespace, space, &
          gr%mesh, elf(:, is), fn_unit, ierr, ions = ions, grp = st%dom_st_kpt_mpi_grp)
      end do
    end select
    SAFE_DEALLOCATE_A(elf)
  end if

  if (outp%what_now(OPTION__OUTPUT__DOS, iter)) then
    call dos_init(dos, namespace, st, hm%kpoints)
    call dos_write_dos(dos, trim(dir), st, gr%sb, ions, gr%mesh, hm, namespace)
  end if

  if (outp%what_now(OPTION__OUTPUT__TPA, iter)) then
    call states_elec_write_tpa(trim(dir), namespace, gr, st)
  end if

  if (outp%what_now(OPTION__OUTPUT__MMB_DEN, iter) .or. outp%what_now(OPTION__OUTPUT__MMB_WFS, iter)) then
    if (states_are_real(st)) then
      call doutput_modelmb(outp, namespace, space, trim(dir), gr, st, ions)
    else
      call zoutput_modelmb(outp, namespace, space, trim(dir), gr, st, ions)
    end if
  end if

  POP_SUB(output_states)

end subroutine output_states


! ---------------------------------------------------------
subroutine output_current_flow(outp, namespace, dir, gr, st, kpoints)
  type(output_t),       intent(in) :: outp
  type(namespace_t),    intent(in) :: namespace
  character(len=*),     intent(in) :: dir
  type(grid_t),         intent(in) :: gr
  type(states_elec_t),  intent(in) :: st
  type(kpoints_t),      intent(in) :: kpoints

  integer :: iunit, ip, idir, rankmin
  FLOAT   :: flow, dmin
  FLOAT, allocatable :: j(:, :, :)

  PUSH_SUB(output_current_flow)

  if(mpi_grp_is_root(mpi_world)) then

    call io_mkdir(dir, namespace)
    iunit = io_open(trim(dir)//'/'//'current-flow', namespace, action='write')

    select case(gr%sb%dim)
    case(3)
      write(iunit,'(a)')        '# Plane:'
      write(iunit,'(3a,3f9.5)') '# origin [', trim(units_abbrev(units_out%length)), '] = ', &
        (units_from_atomic(units_out%length, outp%plane%origin(idir)), idir = 1, 3)
      write(iunit,'(a,3f9.5)')  '# u = ', outp%plane%u(1), outp%plane%u(2), outp%plane%u(3)
      write(iunit,'(a,3f9.5)')  '# v = ', outp%plane%v(1), outp%plane%v(2), outp%plane%v(3)
      write(iunit,'(a,3f9.5)')  '# n = ', outp%plane%n(1), outp%plane%n(2), outp%plane%n(3)
      write(iunit,'(a, f9.5)')  '# spacing = ', units_from_atomic(units_out%length, outp%plane%spacing)
      write(iunit,'(a,2i4)')    '# nu, mu = ', outp%plane%nu, outp%plane%mu
      write(iunit,'(a,2i4)')    '# nv, mv = ', outp%plane%nv, outp%plane%mv

    case(2)
      write(iunit,'(a)')        '# Line:'
      write(iunit,'(3a,2f9.5)') '# origin [',  trim(units_abbrev(units_out%length)), '] = ', &
        (units_from_atomic(units_out%length, outp%line%origin(idir)), idir = 1, 2)
      write(iunit,'(a,2f9.5)')  '# u = ', outp%line%u(1), outp%line%u(2)
      write(iunit,'(a,2f9.5)')  '# n = ', outp%line%n(1), outp%line%n(2)
      write(iunit,'(a, f9.5)')  '# spacing = ', units_from_atomic(units_out%length, outp%line%spacing)
      write(iunit,'(a,2i4)')    '# nu, mu = ', outp%line%nu, outp%line%mu

    case(1)
      write(iunit,'(a)')        '# Point:'
      write(iunit,'(3a, f9.5)') '# origin [',  trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, outp%line%origin(1))

    end select
  end if

  if(states_are_complex(st)) then
    SAFE_ALLOCATE(j(1:gr%mesh%np, 1:gr%sb%dim, 1:st%d%nspin))
    call states_elec_calc_quantities(gr%der, st, kpoints, .false., paramagnetic_current = j)

    do idir = 1, gr%sb%dim
      do ip = 1, gr%mesh%np
        j(ip, idir, 1) = sum(j(ip, idir, 1:st%d%nspin))
      end do
    end do

    select case(gr%sb%dim)
    case(3); flow = mf_surface_integral(gr%mesh, j(:, :, 1), outp%plane)
    case(2); flow = mf_line_integral(gr%mesh, j(:, :, 1), outp%line)
    case(1); flow = j(mesh_nearest_point(gr%mesh, outp%line%origin(1), dmin, rankmin), 1, 1)
    end select

    SAFE_DEALLOCATE_A(j)
  else
    flow = M_ZERO
  end if

  if(mpi_grp_is_root(mpi_world)) then
    write(iunit,'(3a,e20.12)') '# Flow [', trim(units_abbrev(unit_one/units_out%time)), '] = ', &
      units_from_atomic(unit_one/units_out%time, flow)
    call io_close(iunit)
  end if

  POP_SUB(output_current_flow)
end subroutine output_current_flow

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
