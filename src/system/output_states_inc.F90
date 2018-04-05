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
subroutine output_states(st, gr, geo, dir, outp)
  type(states_t),         intent(inout) :: st
  type(grid_t),           intent(inout) :: gr
  type(geometry_t),       intent(in)    :: geo
  character(len=*),       intent(in)    :: dir
  type(output_t),         intent(in)    :: outp

  integer :: ik, ist, idim, idir, is, ierr, ip, nspin
  character(len=MAX_PATH_LEN) :: fname
  type(unit_t) :: fn_unit
  type(base_density_iterator_t)        :: iter
  type(base_density_t),        pointer :: subsys_density
  type(base_density_t),        pointer :: base_density
  character(len=BASE_DENSITY_NAME_LEN) :: name
  FLOAT,       dimension(:,:), pointer :: density
  FLOAT, allocatable :: dtmp(:), elf(:,:), polarization(:, :)
  CMPLX, allocatable :: ztmp(:)
  type(dos_t) :: dos

  PUSH_SUB(output_states)

  nullify(subsys_density, base_density, density)
  if(iand(outp%what, OPTION__OUTPUT__DENSITY) /= 0) then
    fn_unit = units_out%length**(-gr%mesh%sb%dim)
    if(associated(st%subsys_st))then
      call base_states_get(st%subsys_st, subsys_density)
      ASSERT(associated(subsys_density))
      call base_density_get(subsys_density, density)
      ASSERT(associated(density))
    else
      density => st%rho
    end if
    do is = 1, st%d%nspin
      if(st%d%nspin == 1) then
        write(fname, '(a)') 'density'
      else
        write(fname, '(a,i1)') 'density-sp', is
      end if
      if(associated(st%zrho%Im)) then !cmplxscl
        SAFE_ALLOCATE(ztmp(1:gr%fine%mesh%np))
        forall (ip = 1:gr%fine%mesh%np) ztmp(ip) = st%zrho%Re(ip, is) + M_zI * st%zrho%Im(ip, is)
        call zio_function_output(outp%how, dir, fname, gr%fine%mesh, &
          ztmp(:), fn_unit, ierr, geo = geo, grp = st%dom_st_kpt_mpi_grp)
        SAFE_DEALLOCATE_A(ztmp)
      else
        call dio_function_output(outp%how, dir, fname, gr%fine%mesh, &
          density(:, is), fn_unit, ierr, geo = geo, grp = st%dom_st_kpt_mpi_grp)
      end if
    end do
    if(associated(subsys_density))then
      call base_density_init(iter, subsys_density)
      do
        nullify(base_density, density)
        call base_density_next(iter, name, base_density, ierr)
        if(ierr/=BASE_DENSITY_OK)exit
        ASSERT(associated(base_density))
        call base_density_get(base_density, density)
        ASSERT(associated(density))
        call base_density_get(base_density, nspin=nspin)
        ASSERT(nspin>0)
        do is = 1, nspin
          if(nspin>1) then
            write(fname, "(a,i1,'-',a)") "density-sp", is, trim(adjustl(name))
          else
            write(fname, "(a,'-',a)") "density", trim(adjustl(name))
          end if
          call dio_function_output(outp%how, dir, fname, gr%fine%mesh, &
            density(:,is), fn_unit, ierr, geo=geo, grp=st%dom_st_kpt_mpi_grp)
        end do
      end do 
      call base_density_end(iter)
      nullify(subsys_density, base_density, density)
    end if    
  end if

  if(iand(outp%what, OPTION__OUTPUT__POL_DENSITY) /= 0) then
    fn_unit = units_out%length**(1-gr%mesh%sb%dim)
    SAFE_ALLOCATE(polarization(1:gr%fine%mesh%np, 1:gr%sb%dim))

    do is = 1, st%d%nspin
      forall(ip = 1:gr%fine%mesh%np, idir = 1:gr%sb%dim) polarization(ip, idir) = st%rho(ip, is)*gr%fine%mesh%x(ip, idir)

      if(st%d%nspin == 1) then
        write(fname, '(a)') 'dipole_density'
      else
        write(fname, '(a,i1)') 'dipole_density-sp', is
      end if
      call io_function_output_vector(outp%how, dir, fname, gr%fine%mesh, polarization, gr%sb%dim, fn_unit, ierr, &
        geo = geo, grp = st%dom_st_kpt_mpi_grp, vector_dim_labels = (/'x', 'y', 'z'/))
    end do
    
    SAFE_DEALLOCATE_A(polarization)
  end if

  if(iand(outp%what, OPTION__OUTPUT__WFS) /= 0) then
    fn_unit = sqrt(units_out%length**(-gr%mesh%sb%dim))

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

            if (states_are_real(st)) then
              call states_get_state(st, gr%mesh, idim, ist, ik, dtmp)
              call dio_function_output(outp%how, dir, fname, gr%mesh, dtmp, &
                fn_unit, ierr, geo = geo)
            else
              call states_get_state(st, gr%mesh, idim, ist, ik, ztmp)
              call zio_function_output(outp%how, dir, fname, gr%mesh, ztmp, &
                fn_unit, ierr, geo = geo)
              if(st%have_left_states) then
                fname = 'L'// trim(fname)
                call states_get_state(st, gr%mesh, idim, ist, ik, ztmp, left = .true.)
                call zio_function_output(outp%how, dir, fname, gr%mesh, ztmp, &
                  fn_unit, ierr, geo = geo)                    
              end if
            end if
          end do
        end do
      end if
    end do

    SAFE_DEALLOCATE_A(dtmp)
    SAFE_DEALLOCATE_A(ztmp)
  end if

  if(iand(outp%what, OPTION__OUTPUT__WFS_SQMOD) /= 0) then
    fn_unit = units_out%length**(-gr%mesh%sb%dim)
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
              call states_get_state(st, gr%mesh, idim, ist, ik, dtmp)
              dtmp(1:gr%mesh%np) = abs(dtmp(1:gr%mesh%np))**2
            else
              call states_get_state(st, gr%mesh, idim, ist, ik, ztmp)
              dtmp(1:gr%mesh%np) = abs(ztmp(1:gr%mesh%np))**2
            end if
            call dio_function_output (outp%how, dir, fname, gr%mesh, &
              dtmp, fn_unit, ierr, geo = geo)
          end do
        end do
      end if
    end do
    SAFE_DEALLOCATE_A(dtmp)
    SAFE_DEALLOCATE_A(ztmp)
  end if

  if(iand(outp%what, OPTION__OUTPUT__KINETIC_ENERGY_DENSITY) /= 0) then
    fn_unit = units_out%energy * units_out%length**(-gr%mesh%sb%dim)
    SAFE_ALLOCATE(elf(1:gr%mesh%np, 1:st%d%nspin))
    call states_calc_quantities(gr%der, st, kinetic_energy_density = elf)
    select case(st%d%ispin)
    case(UNPOLARIZED)
      write(fname, '(a)') 'tau'
      call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
        elf(:,1), unit_one, ierr, geo = geo, grp = st%dom_st_kpt_mpi_grp)
    case(SPIN_POLARIZED, SPINORS)
      do is = 1, 2
        write(fname, '(a,i1)') 'tau-sp', is
        call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
          elf(:, is), unit_one, ierr, geo = geo, grp = st%dom_st_kpt_mpi_grp)
      end do
    end select
    SAFE_DEALLOCATE_A(elf)
  end if

  if(iand(outp%what, OPTION__OUTPUT__DOS) /= 0) then
    call dos_init(dos, st)
    call dos_write_dos (dos, trim(dir), st, gr%sb, geo, gr%mesh)
    call dos_end(dos)
  end if

  if(iand(outp%what, OPTION__OUTPUT__TPA) /= 0) then
    call states_write_tpa (trim(dir), gr, st)
  end if

  if(iand(outp%what, OPTION__OUTPUT__MMB_DEN) /= 0 .or. iand(outp%what, OPTION__OUTPUT__MMB_WFS) /= 0) then
    if (states_are_real(st)) then
      call doutput_modelmb (trim(dir), gr, st, geo, outp)
    else
      call zoutput_modelmb (trim(dir), gr, st, geo, outp)
    end if
  end if

  POP_SUB(output_states)

end subroutine output_states


! ---------------------------------------------------------
subroutine output_current_flow(gr, st, dir, outp)
  type(grid_t),         intent(inout) :: gr
  type(states_t),       intent(inout) :: st
  character(len=*),     intent(in)    :: dir
  type(output_t),       intent(in)    :: outp

  integer :: iunit, ip, idir, rankmin
  FLOAT   :: flow, dmin
  FLOAT, allocatable :: j(:, :, :)

  PUSH_SUB(output_current_flow)

  if(iand(outp%what, OPTION__OUTPUT__J_FLOW) == 0) then
    POP_SUB(output_current_flow)
    return
  end if

  if(mpi_grp_is_root(mpi_world)) then

    call io_mkdir(dir)
    iunit = io_open(trim(dir)//'/'//'current-flow', action='write')

    select case(gr%mesh%sb%dim)
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
    SAFE_ALLOCATE(j(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%nspin))
    call states_calc_quantities(gr%der, st, paramagnetic_current = j)

    do idir = 1, gr%mesh%sb%dim
      do ip = 1, gr%mesh%np
        j(ip, idir, 1) = sum(j(ip, idir, 1:st%d%nspin))
      end do
    end do

    select case(gr%mesh%sb%dim)
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
