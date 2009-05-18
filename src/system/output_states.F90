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
  subroutine h_sys_output_states(st, gr, geo, dir, outp)

    type(states_t),         intent(inout) :: st
    type(grid_t),           intent(inout) :: gr
    type(geometry_t),       intent(in)    :: geo
    character(len=*),       intent(in)    :: dir
    type(h_sys_output_t),   intent(in)    :: outp

    integer :: ik, ist, idim, is, id, ierr
    character(len=80) :: fname, dirname
    FLOAT :: u
    FLOAT, allocatable :: dtmp(:), elf(:,:)
    FLOAT, allocatable :: current(:, :, :)

    call push_sub('output_states.h_sys_output_states')

    u = M_ONE/units_out%length%factor**gr%mesh%sb%dim

    if(iand(outp%what, output_density).ne.0) then
      do is = 1, st%d%nspin
        write(fname, '(a,i1)') 'density-', is
        call doutput_function(outp%how, dir, fname, gr%mesh, gr%sb, &
          st%rho(:, is), u, ierr, is_tmp = .false., geo = geo, grp = st%mpi_grp)
      end do
    end if

    if(iand(outp%what, output_pol_density).ne.0) then
      SAFE_ALLOCATE(dtmp(1:gr%mesh%np))
      do idim = 1, gr%mesh%sb%dim
        do is = 1, st%d%nspin
          dtmp(1:gr%mesh%np)=st%rho(1:gr%mesh%np,is)*gr%mesh%x(1:gr%mesh%np,idim)
          write(fname, '(a,i1,a,i1)') 'dipole_density-', is, '-',idim
          call doutput_function(outp%how, dir, fname, gr%mesh, gr%sb, &
            dtmp(:), u, ierr, is_tmp = .false., geo = geo, grp = st%mpi_grp)
        end do
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    if( (iand(outp%what, output_current).ne.0) .and. (st%wfs_type == M_CMPLX) ) then
      ! calculate current first
      SAFE_ALLOCATE(current(1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:st%d%nspin))
      call states_calc_tau_jp_gn(gr, st, jp = current)
      do is = 1, st%d%nspin
        do id = 1, gr%mesh%sb%dim
          write(fname, '(a,i1,a,a)') 'current-', is, '-', index2axis(id)
          call doutput_function(outp%how, dir, fname, gr%mesh, gr%sb, &
            current(:, id, is), u, ierr, is_tmp = .false., geo = geo, grp = st%mpi_grp)
        end do
      end do
      SAFE_DEALLOCATE_A(current)
    end if

    if(iand(outp%what, output_wfs).ne.0) then
      do ist = st%st_start, st%st_end
        if(loct_isinstringlist(ist, outp%wfs_list)) then
          do ik = st%d%kpt%start, st%d%kpt%end
            do idim = 1, st%d%dim
              write(fname, '(a,i3.3,a,i4.4,a,i1)') 'wf-', ik, '-', ist, '-', idim
              if (st%wfs_type == M_REAL) then
                call doutput_function(outp%how, dir, fname, gr%mesh, gr%sb, &
                     st%dpsi(1:, idim, ist, ik), sqrt(u), ierr, is_tmp = .false., geo = geo)
              else
                call zoutput_function(outp%how, dir, fname, gr%mesh, gr%sb, &
                     st%zpsi(1:, idim, ist, ik), sqrt(u), ierr, is_tmp = .false., geo = geo)
              end if
            end do
          end do
        end if
      end do
    end if

    if(iand(outp%what, output_wfs_sqmod).ne.0) then
      SAFE_ALLOCATE(dtmp(1:gr%mesh%np_part))
      do ist = st%st_start, st%st_end
        if(loct_isinstringlist(ist, outp%wfs_list)) then
          do ik = st%d%kpt%start, st%d%kpt%end
            do idim = 1, st%d%dim
              write(fname, '(a,i3.3,a,i4.4,a,i1)') 'sqm-wf-', ik, '-', ist, '-', idim
              if (st%wfs_type == M_REAL) then
                dtmp = abs(st%dpsi(:, idim, ist, ik))**2
              else
                dtmp = abs(st%zpsi(:, idim, ist, ik))**2
              end if
              call doutput_function (outp%how, dir, fname, gr%mesh, gr%sb, &
                dtmp, u, ierr, is_tmp = .false., geo = geo)
            end do
          end do
        end if
      end do
      SAFE_DEALLOCATE_A(dtmp)
    end if

    if(iand(outp%what, output_ked).ne.0) then
      SAFE_ALLOCATE(elf(1:gr%mesh%np, 1:st%d%nspin))
      call states_calc_tau_jp_gn(gr, st, tau=elf)
      select case(st%d%ispin)
        case(UNPOLARIZED)
          write(fname, '(a)') 'tau'
          call doutput_function(outp%how, dir, trim(fname), gr%mesh, gr%sb, &
            elf(:,1), M_ONE, ierr, is_tmp = .false., geo = geo)
        case(SPIN_POLARIZED, SPINORS)
          do is = 1, 2
            write(fname, '(a,a,i1)') 'tau', '-', is
            call doutput_function(outp%how, dir, trim(fname), gr%mesh, gr%sb, &
              elf(:, is), M_ONE, ierr, is_tmp = .false., geo = geo, grp = st%mpi_grp)
          end do
      end select
      SAFE_DEALLOCATE_A(elf)
    end if

    if(iand(outp%what, output_dos).ne.0) then
      call states_write_dos (trim(dir), st)
    end if

    if(iand(outp%what, output_tpa).ne.0) then
      call states_write_tpa (trim(dir), gr, st)
    end if

    if(iand(outp%what, output_modelmb).ne.0) then
      call h_sys_output_modelmb (trim(dir), gr, st, geo, outp)
    end if

    call pop_sub()

  end subroutine h_sys_output_states


  !
  !  routine for output of model many-body quantities.
  !
  subroutine h_sys_output_modelmb (dir, gr, st, geo, outp)

    type(states_t),         intent(inout) :: st
    type(grid_t),           intent(inout) :: gr
    character(len=*),       intent(in)    :: dir
    type(geometry_t),       intent(in)    :: geo
    type(h_sys_output_t),   intent(in)    :: outp

    ! local vars
    integer :: mm
    logical :: symmetries_satisfied, impose_exch_symmetry
    CMPLX, allocatable :: wf(:)
    character(len=80) :: dirname
    type(modelmb_denmat_t) :: denmat
    type(modelmb_density_t) :: den

    call push_sub('system.h_sys_output_modelmb')

    impose_exch_symmetry = .true.

    ! make sure directory exists
    call loct_mkdir(trim(dir))
    ! all model mb stuff should be in this directory
    dirname = trim(dir)//'/modelmb'
    call loct_mkdir(trim(dirname))

    SAFE_ALLOCATE(wf(1:gr%mesh%np_part_global))

    call modelmb_density_matrix_nullify(denmat)
    if(iand(outp%what, output_density_matrix).ne.0) then
      call modelmb_density_matrix_init(dirname, st, denmat)
    end if
 
    call modelmb_density_nullify(den)
    if(iand(outp%what, output_density).ne.0) then
      call modelmb_density_init (dirname, st, den)
    end if
 
    do mm = 1, st%nst
      ! NOTE!!!! do not make this into some preprocessed X() stuff until I am dead and
      ! buried. Thanks - mjv
      if(states_are_real(st)) then
        wf = cmplx(st%dpsi(1:gr%mesh%np_part_global, 1, mm, 1),M_ZERO)
      else
        wf = st%zpsi(1:gr%mesh%np_part_global, 1, mm, 1)
      end if

      if (impose_exch_symmetry) then
        call modelmb_sym_state(dirname, gr, mm, geo, st%modelmbparticles, wf, symmetries_satisfied)
      end if

      if(iand(outp%what, output_density_matrix).ne.0 .and. symmetries_satisfied) then
        call modelmb_density_matrix_write(gr, st, wf, mm, denmat)
      end if

      if(      iand(outp%what, output_density).ne.0 .and. &
         .not. iand(outp%what, output_density_matrix).ne.0 .and. &
         symmetries_satisfied) then
        call modelmb_density_write(gr, st, wf, mm, den)
      end if

      if(iand(outp%what, output_wfs).ne.0 .and. symmetries_satisfied) then
        call zio_function_out_text(trim(dirname), gr%mesh, mm, wf)
      end if

    end do

    SAFE_DEALLOCATE_A(wf)

    if(iand(outp%what, output_density_matrix).ne.0) then
      call modelmb_density_matrix_end (denmat)
    end if

    if(iand(outp%what, output_density).ne.0) then
      call modelmb_density_end (den)
    end if
 
    call pop_sub()

  end subroutine h_sys_output_modelmb


  ! ---------------------------------------------------------
  ! Prints out the multipole matrix elements between KS states.
  ! It prints the states to the file opened in iunit.
  ! It prints the (l,m) multipole moment, for
  ! the Kohn-Sham states in the irreducible subspace ik.
  ! ---------------------------------------------------------
  subroutine h_sys_write_multipole_matrix(st, gr, l, m, ik, iunit, geo)
    type(states_t), intent(in) :: st
    type(grid_t), intent(in) :: gr
    integer, intent(in) :: l, m, ik, iunit
    type(geometry_t), intent(in)    :: geo

    integer :: ii, jj, i
    FLOAT, allocatable :: multipole(:, :)
    CMPLX :: multip_element
    FLOAT :: r, x(MAX_DIM), ylm

    call push_sub('output_states.h_sys_write_multipole_matrix')

    write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | r**l * Y_{lm}(theta,phi) | Phi_j>' 
    write(iunit, fmt = '(a,i2,a,i2)') '# l =', l, '; m =', m
    write(iunit, fmt = '(a,i4)')      '# ik =', ik
    if(l>1) then
      write(iunit, fmt = '(a,i1)') '# Units = ['//trim(units_out%length%abbrev)//']^',l
    else
      write(iunit, fmt = '(a)')    '# Units = ['//trim(units_out%length%abbrev)//']'
    end if

    SAFE_ALLOCATE(multipole(1:gr%mesh%np_part, 1:st%d%dim))
    multipole = M_ZERO
    do ii = 1, st%d%dim
      do i = 1, gr%mesh%np
        call mesh_r(gr%mesh, i, r, x = x)
        call loct_ylm(1, x(1), x(2), x(3), l, m, ylm)
        multipole(i, ii) = r**l * ylm
      end do
    end do

    do ii = 1, st%nst
      do jj = 1, st%nst
        if (st%wfs_type == M_REAL) then 
          write(iunit,fmt = '(f20.10)', advance = 'no') dmf_dotp(gr%mesh, st%d%dim, &
            st%dpsi(:, :, ii, 1), &
            st%dpsi(:, :, jj, 1) * multipole(:, :)) / units_out%length%factor**l

        else
          multip_element = zmf_dotp(gr%mesh, st%d%dim, &
            st%zpsi(:, :, ii, 1), &
            st%zpsi(:, :, jj, 1) * multipole(:, :)) / units_out%length%factor**l

          write(iunit,fmt = '(f20.12,a1,f20.12,3x)', advance = 'no') &
            real(multip_element, REAL_PRECISION), ',', aimag(multip_element)
        end if
        if(jj==st%nst) write(iunit, '(a)') 
      end do
    end do

    SAFE_DEALLOCATE_A(multipole)
    call pop_sub()
  end subroutine h_sys_write_multipole_matrix


  ! ---------------------------------------------------------
  subroutine h_sys_output_current_flow(gr, st, dir, outp, geo)
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    character(len=*),     intent(in)    :: dir
    type(h_sys_output_t), intent(in)    :: outp
    type(geometry_t),     intent(in)    :: geo

    integer :: iunit, i, k, rankmin
    FLOAT   :: flow, dmin
    FLOAT, allocatable :: j(:, :, :)

    call push_sub('output_states.h_sys_output_current_flow')

    if(iand(outp%what, output_j_flow) == 0) then
      call pop_sub(); return
    end if

    if(mpi_grp_is_root(mpi_world)) then

      iunit = io_open(trim(dir)//'/'//'current-flow', action='write')

      select case(gr%mesh%sb%dim)
      case(3)
        write(iunit,'(a)')       '# Plane:'
        write(iunit,'(a,3f9.5)') '# u = ', outp%plane%u(1), outp%plane%u(2), outp%plane%u(3)
        write(iunit,'(a,3f9.5)') '# v = ', outp%plane%v(1), outp%plane%v(2), outp%plane%v(3)
        write(iunit,'(a,3f9.5)') '# n = ', outp%plane%n(1), outp%plane%n(2), outp%plane%n(3)
        write(iunit,'(a, f9.5)') '# spacing = ', outp%plane%spacing
        write(iunit,'(a,2i4)')   '# nu, mu = ', outp%plane%nu, outp%plane%mu
        write(iunit,'(a,2i4)')   '# nv, mv = ', outp%plane%nv, outp%plane%mv

      case(2)
        write(iunit,'(a)')       '# Line:'
        write(iunit,'(a,2f9.5)') '# u = ', outp%line%u(1), outp%line%u(2)
        write(iunit,'(a,2f9.5)') '# n = ', outp%line%n(1), outp%line%n(2)
        write(iunit,'(a, f9.5)') '# spacing = ', outp%line%spacing
        write(iunit,'(a,2i4)')   '# nu, mu = ', outp%line%nu, outp%line%mu

      case(1)
        write(iunit,'(a)')       '# Point:'
        write(iunit,'(a, f9.5)') '# origin = ', outp%line%origin(1)

      end select
    end if

    if(states_are_complex(st)) then
      SAFE_ALLOCATE(j(1:gr%mesh%np, 1:MAX_DIM, 0:st%d%nspin))
      call states_calc_tau_jp_gn(gr, st, jp = j(:, :, 1:))

      do k = 1, gr%mesh%sb%dim
        do i = 1, gr%mesh%np
          j(i, k, 0) = sum(j(i, k, 1:st%d%nspin))
        end do
      end do

      select case(gr%mesh%sb%dim)
      case(3); flow = mf_surface_integral (gr%mesh, j(:, :, 0), outp%plane)
      case(2); flow = mf_line_integral (gr%mesh, j(:, :, 0), outp%line)
      case(1); flow = j(mesh_nearest_point(gr%mesh, outp%line%origin(1), dmin, rankmin), 1, 0)
      end select

      SAFE_DEALLOCATE_A(j)
    else
      flow = M_ZERO
    end if

    if(mpi_grp_is_root(mpi_world)) write(iunit,'(a,e20.12)') '# Flow = ', flow

    call io_close(iunit)
    call pop_sub()
  end subroutine h_sys_output_current_flow

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
