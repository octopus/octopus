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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: states.F90 2515 2006-10-24 17:13:30Z acastro $

  ! ---------------------------------------------------------
  subroutine output_states(st, gr, geo, dir, outp)

    type(states_t),         intent(inout) :: st
    type(grid_t),           intent(inout) :: gr
    type(geometry_t),       intent(in)    :: geo
    character(len=*),       intent(in)    :: dir
    type(output_t),         intent(in)    :: outp

    integer :: ik, ist, idim, idir, is, ierr, ip
    character(len=80) :: fname
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: dtmp(:), elf(:,:)
    FLOAT, allocatable :: current(:, :, :)

    PUSH_SUB(output_states)

    if(iand(outp%what, C_OUTPUT_DENSITY) .ne. 0) then
      fn_unit = units_out%length**(-gr%mesh%sb%dim)
      do is = 1, st%d%nspin
        if(st%d%nspin == 1) then
          write(fname, '(a)') 'density'
        else
          write(fname, '(a,i1)') 'density-sp', is
        endif
        call dio_function_output(outp%how, dir, fname, gr%fine%mesh, &
          st%rho(:, is), fn_unit, ierr, is_tmp = .false., geo = geo, grp = st%mpi_grp)
      end do
    end if

    if(iand(outp%what, C_OUTPUT_POL_DENSITY) .ne. 0) then
      fn_unit = units_out%length**(1-gr%mesh%sb%dim)
      SAFE_ALLOCATE(dtmp(1:gr%fine%mesh%np))
      do idir = 1, gr%sb%dim
        do is = 1, st%d%nspin
          forall (ip = 1:gr%fine%mesh%np) dtmp(ip) = st%rho(ip, is) * gr%fine%mesh%x(ip, idir)
          if(st%d%nspin == 1) then
            write(fname, '(2a)') 'dipole_density-', index2axis(idir)
          else
            write(fname, '(a,i1,2a)') 'dipole_density-sp', is, '-', index2axis(idir)
          endif
          call dio_function_output(outp%how, dir, fname, gr%fine%mesh, &
            dtmp(:), fn_unit, ierr, is_tmp = .false., geo = geo, grp = st%mpi_grp)
        end do
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    if(iand(outp%what, C_OUTPUT_CURRENT) .ne. 0) then
      if(states_are_complex(st)) then
        fn_unit = units_out%time * units_out%length**(-gr%mesh%sb%dim)
        ! calculate current first
        SAFE_ALLOCATE(current(1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:st%d%nspin))
        call states_calc_quantities(gr%der, st, paramagnetic_current = current)
        do is = 1, st%d%nspin
          do idir = 1, gr%mesh%sb%dim
            if(st%d%nspin == 1) then
              write(fname, '(2a)') 'current-', index2axis(idir)
            else
              write(fname, '(a,i1,2a)') 'current-sp', is, '-', index2axis(idir)
            endif
            call dio_function_output(outp%how, dir, fname, gr%mesh, &
              current(:, idir, is), fn_unit, ierr, is_tmp = .false., geo = geo, grp = st%mpi_grp)
          end do
        end do
        SAFE_DEALLOCATE_A(current)
      else
        message(1) = 'No current density output for real states since it is identically zero.'
        call messages_warning(1)
      endif
    end if

    if(iand(outp%what, C_OUTPUT_WFS).ne.0) then
      fn_unit = sqrt(units_out%length**(-gr%mesh%sb%dim))
      do ist = st%st_start, st%st_end
        if(loct_isinstringlist(ist, outp%wfs_list)) then
          do ik = st%d%kpt%start, st%d%kpt%end
            do idim = 1, st%d%dim
              if(st%d%nik > 1) then
                if(st%d%dim > 1) then
                  write(fname, '(a,i3.3,a,i4.4,a,i1)') 'wf-k', ik, '-st', ist, '-sp', idim
                else
                  write(fname, '(a,i3.3,a,i4.4)')      'wf-k', ik, '-st', ist
                endif
              else
                if(st%d%dim > 1) then
                  write(fname, '(a,i4.4,a,i1)')        'wf-st', ist, '-sp', idim
                else
                  write(fname, '(a,i4.4)')             'wf-st', ist
                endif
              endif
                
              if (states_are_real(st)) then
                call dio_function_output(outp%how, dir, fname, gr%mesh, &
                     st%dpsi(1:, idim, ist, ik), fn_unit, ierr, is_tmp = .false., geo = geo)
              else
                call zio_function_output(outp%how, dir, fname, gr%mesh, &
                     st%zpsi(1:, idim, ist, ik), fn_unit, ierr, is_tmp = .false., geo = geo)
              end if
            end do
          end do
        end if
      end do
    end if

    if(iand(outp%what, C_OUTPUT_WFS_SQMOD).ne.0) then
      fn_unit = units_out%length**(-gr%mesh%sb%dim)
      SAFE_ALLOCATE(dtmp(1:gr%mesh%np_part))
      do ist = st%st_start, st%st_end
        if(loct_isinstringlist(ist, outp%wfs_list)) then
          do ik = st%d%kpt%start, st%d%kpt%end
            do idim = 1, st%d%dim
              if(st%d%nik > 1) then
                if(st%d%dim > 1) then
                  write(fname, '(a,i3.3,a,i4.4,a,i1)') 'sqm-wf-k', ik, '-st', ist, '-sp', idim
                else
                  write(fname, '(a,i3.3,a,i4.4)')      'sqm-wf-k', ik, '-st', ist
                endif
              else
                if(st%d%dim > 1) then
                  write(fname, '(a,i4.4,a,i1)')        'sqm-wf-st', ist, '-sp', idim
                else
                  write(fname, '(a,i4.4)')             'sqm-wf-st', ist
                endif
              endif

              if (states_are_real(st)) then
                dtmp = abs(st%dpsi(:, idim, ist, ik))**2
              else
                dtmp = abs(st%zpsi(:, idim, ist, ik))**2
              end if
              call dio_function_output (outp%how, dir, fname, gr%mesh, &
                dtmp, fn_unit, ierr, is_tmp = .false., geo = geo)
            end do
          end do
        end if
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    if(iand(outp%what, C_OUTPUT_KED).ne.0) then
      fn_unit = units_out%energy * units_out%length**(-gr%mesh%sb%dim)
      SAFE_ALLOCATE(elf(1:gr%mesh%np, 1:st%d%nspin))
      call states_calc_quantities(gr%der, st, kinetic_energy_density = elf)
      select case(st%d%ispin)
        case(UNPOLARIZED)
          write(fname, '(a)') 'tau'
          call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
            elf(:,1), unit_one, ierr, is_tmp = .false., geo = geo)
        case(SPIN_POLARIZED, SPINORS)
          do is = 1, 2
            write(fname, '(a,i1)') 'tau-sp', is
            call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
              elf(:, is), unit_one, ierr, is_tmp = .false., geo = geo, grp = st%mpi_grp)
          end do
      end select
      SAFE_DEALLOCATE_A(elf)
    end if

    if(iand(outp%what, C_OUTPUT_DOS).ne.0) then
      call states_write_dos (trim(dir), st)
    end if

    if(iand(outp%what, C_OUTPUT_TPA).ne.0) then
      call states_write_tpa (trim(dir), gr, st)
    end if

    if(iand(outp%what, C_OUTPUT_MODELMB).ne.0) then
      call output_modelmb (trim(dir), gr, st, geo, outp)
    end if

    POP_SUB(output_states)

  end subroutine output_states


  ! ---------------------------------------------------------
  !
  !  routine for output of model many-body quantities.
  !
  subroutine output_modelmb (dir, gr, st, geo, outp)
    type(states_t),         intent(inout) :: st
    type(grid_t),           intent(inout) :: gr
    character(len=*),       intent(in)    :: dir
    type(geometry_t),       intent(in)    :: geo
    type(output_t),         intent(in)    :: outp

    ! local vars
    integer :: mm, iunit, itype
    integer :: iyoung
    logical :: symmetries_satisfied, impose_exch_symmetry
    CMPLX, allocatable :: wf(:)
    character(len=80) :: dirname
    character(len=80) :: filename
    type(modelmb_denmat_t) :: denmat
    type(modelmb_density_t) :: den

    PUSH_SUB(output_modelmb)

    impose_exch_symmetry = .true.

    ! make sure directory exists
    call io_mkdir(trim(dir))
    ! all model mb stuff should be in this directory
    dirname = trim(dir)//'modelmb'
    call io_mkdir(trim(dirname))

    SAFE_ALLOCATE(wf(1:gr%mesh%np_part_global))

    call modelmb_density_matrix_nullify(denmat)
    if(iand(outp%what, C_OUTPUT_DENSITY_MATRIX).ne.0) then
      call modelmb_density_matrix_init(dirname, st, denmat)
    end if
 
    call modelmb_density_nullify(den)
    if(iand(outp%what, C_OUTPUT_DENSITY).ne.0) then
      call modelmb_density_init (dirname, st, den)
    end if
 
    ! open file for Young diagrams and projection info
    write (filename,'(a,a)') trim(dirname), '/youngprojections'
    iunit = io_open(trim(filename), action='write')

    ! just treat particle type 1 for the moment
    itype = 1
    call young_write_allspins (iunit, st%modelmbparticles%nparticles_per_type(itype))

    ! write header
    write (iunit, '(a)') '  state      eigenvalue   ptype    Young#    nspindown    projection'

    iyoung = 1
    do mm = 1, st%nst
      ! FIXME make this into some preprocessed X() stuff, along with dens and dens_mat
      if(states_are_real(st)) then
        wf = cmplx(st%dpsi(1:gr%mesh%np_part_global, 1, mm, 1), M_ZERO)
      else
        wf = st%zpsi(1:gr%mesh%np_part_global, 1, mm, 1)
      end if

      if (impose_exch_symmetry) then
        if (mm > 1) then
          ! if eigenval is degenerate increment iyoung
          if (abs(st%eigenval(mm,1) - st%eigenval(mm-1,1)) < 1.e-5) then
            iyoung = iyoung + 1
          else
            iyoung = 1
          end if
        end if
 
        call modelmb_sym_state(st%eigenval(mm,1), iyoung, iunit, gr, mm, geo, &
             st%modelmbparticles, wf, symmetries_satisfied)
      end if

      if(iand(outp%what, C_OUTPUT_DENSITY_MATRIX).ne.0 .and. symmetries_satisfied) then
        call modelmb_density_matrix_write(gr, st, wf, mm, denmat)
      end if

      if(      iand(outp%what, C_OUTPUT_DENSITY).ne.0 .and. &
         .not. iand(outp%what, C_OUTPUT_DENSITY_MATRIX).ne.0 .and. &
         symmetries_satisfied) then
        call modelmb_density_write(gr, st, wf, mm, den)
      end if

      if(iand(outp%what, C_OUTPUT_WFS).ne.0 .and. symmetries_satisfied) then
        call zio_function_out_text(trim(dirname), gr%mesh, mm, wf)
      end if

    end do

    call io_close(iunit)

    SAFE_DEALLOCATE_A(wf)

    if(iand(outp%what, C_OUTPUT_DENSITY_MATRIX).ne.0) then
      call modelmb_density_matrix_end (denmat)
    end if

    if(iand(outp%what, C_OUTPUT_DENSITY).ne.0) then
      call modelmb_density_end (den)
    end if
 
    POP_SUB(output_modelmb)

  end subroutine output_modelmb


  ! ---------------------------------------------------------
  subroutine output_current_flow(gr, st, dir, outp, geo)
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    character(len=*),     intent(in)    :: dir
    type(output_t),       intent(in)    :: outp
    type(geometry_t),     intent(in)    :: geo

    integer :: iunit, ip, idir, rankmin
    FLOAT   :: flow, dmin
    FLOAT, allocatable :: j(:, :, :)

    PUSH_SUB(output_current_flow)

    if(iand(outp%what, C_OUTPUT_J_FLOW) == 0) then
      POP_SUB(output_current_flow)
      return
    end if

    if(mpi_grp_is_root(mpi_world)) then

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
             (units_from_atomic(units_out%length, outp%plane%origin(idir)), idir = 1, 2)
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
      SAFE_ALLOCATE(j(1:gr%mesh%np, 1:MAX_DIM, 0:st%d%nspin))
      call states_calc_quantities(gr%der, st, paramagnetic_current = j(:, :, 1:))

      do idir = 1, gr%mesh%sb%dim
        do ip = 1, gr%mesh%np
          j(ip, idir, 0) = sum(j(ip, idir, 1:st%d%nspin))
        end do
      end do

      select case(gr%mesh%sb%dim)
      case(3); flow = mf_surface_integral(gr%mesh, j(:, :, 0), outp%plane)
      case(2); flow = mf_line_integral(gr%mesh, j(:, :, 0), outp%line)
      case(1); flow = j(mesh_nearest_point(gr%mesh, outp%line%origin(1), dmin, rankmin), 1, 0)
      end select

      SAFE_DEALLOCATE_A(j)
    else
      flow = M_ZERO
    end if

    if(mpi_grp_is_root(mpi_world)) then
      write(iunit,'(3a,e20.12)') '# Flow [', trim(units_abbrev(unit_one / units_out%time)), '] = ', &
           units_from_atomic(unit_one / units_out%time, flow)
      call io_close(iunit)
    end if

    POP_SUB(output_current_flow)
  end subroutine output_current_flow

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
