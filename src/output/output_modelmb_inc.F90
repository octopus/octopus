!! Copyright (C) 2009 N. Helbig and M. Verstraete
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
!
!> routine for output of model many-body quantities.
!
subroutine X(output_modelmb) (outp, namespace, dir, gr, st, geo)
  type(output_t),         intent(in)    :: outp
  type(namespace_t),      intent(in)    :: namespace
  character(len=*),       intent(in)    :: dir
  type(grid_t),           intent(in)    :: gr
  type(states_elec_t),    intent(in)    :: st
  type(geometry_t),       intent(in)    :: geo

  integer :: mm
  integer :: ierr, iunit
  integer :: ncombo
  integer :: itype
  integer, allocatable :: ndiagrams(:)
  integer, allocatable :: young_used(:)
  logical :: symmetries_satisfied
  R_TYPE, allocatable :: wf(:)
  character(len=80) :: dirname
  character(len=80) :: filename
  character(len=500) :: youngstring, tmpstring
  type(modelmb_denmat_t) :: denmat
  type(unit_t)  :: fn_unit

  PUSH_SUB(X(output_modelmb))

  ! make sure directory exists
  call io_mkdir(trim(dir), namespace)
  ! all model mb stuff should be in this directory
  dirname = trim(dir)//'/modelmb'
  call io_mkdir(trim(dirname), namespace)

  ! open file for Young diagrams and projection info
  write (filename,'(a,a)') trim(dirname), '/youngprojections'
  iunit = io_open(trim(filename), namespace, action='write')

  ! treat all particle types
  SAFE_ALLOCATE(ndiagrams(1:st%modelmbparticles%ntype_of_particle))
  ndiagrams = 1
  do itype = 1, st%modelmbparticles%ntype_of_particle
    write (iunit, '(a, i6)') '  Young diagrams for particle type ', itype
    call young_write_allspins (iunit, st%modelmbparticles%nparticles_per_type(itype))
    call young_ndiagrams (st%modelmbparticles%nparticles_per_type(itype), ndiagrams(itype))
  end do

  ncombo = product(ndiagrams)
  write (iunit, '(a, i6)') ' # of possible combinations of Young diagrams for all types = ', &
    ncombo

  SAFE_ALLOCATE(young_used(1:ncombo))
  young_used = 0

  ! write header
  write (iunit, '(a)') '  state        eigenvalue         projection   nspindown Young# for each type'

  SAFE_ALLOCATE(wf(1:gr%mesh%np))

  if(bitand(outp%what, OPTION__OUTPUT__MMB_DEN) /= 0) then
    call modelmb_density_matrix_init(dirname, namespace, st, denmat)
  end if

  do mm = 1, st%nst
!TODO : check if there is another interface for get_states to avoid trivial slice of wf
    call states_elec_get_state(st, gr%mesh, 1, mm, 1, wf)

    youngstring = ""
    if (all(st%mmb_nspindown(:,mm) >= 0)) then
      do itype = 1, st%modelmbparticles%ntype_of_particle
        write (tmpstring, '(3x,I4,1x,I4)') st%mmb_nspindown(itype,mm), st%mmb_iyoung(itype,mm)
        youngstring = trim(youngstring) // trim(tmpstring)
      end do
    else 
      youngstring  = " state does not have an associated Young diagram"
    end if
    write (iunit, '(a,I5,3x,E16.6,5x,E14.6,2x,a)') &
      "  ", mm, st%eigenval(mm,1), st%mmb_proj(mm), trim(youngstring)

    symmetries_satisfied = .true.
    if (st%mmb_proj(mm) < CNST(1.e-6)) then
      symmetries_satisfied = .false.
    end if

    if(bitand(outp%what, OPTION__OUTPUT__MMB_DEN) /= 0 .and. symmetries_satisfied) then
      call X(modelmb_density_matrix_write)(gr, st, wf, mm, denmat, namespace)
    end if

    if(bitand(outp%what, OPTION__OUTPUT__MMB_WFS) /= 0 .and. symmetries_satisfied) then
      fn_unit = units_out%length**(-gr%sb%dim)
      write(filename, '(a,i4.4)') 'wf-st', mm
      call X(io_function_output)(outp%how, trim(dirname), trim(filename), namespace, &
        gr%mesh, wf, fn_unit, ierr, geo = geo)
    end if

  end do

  call io_close(iunit)

  SAFE_DEALLOCATE_A(wf)

  if(bitand(outp%what, OPTION__OUTPUT__MMB_DEN) /= 0) then
    call modelmb_density_matrix_end (denmat)
  end if

  POP_SUB(X(output_modelmb))

end subroutine X(output_modelmb)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
