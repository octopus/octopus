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
subroutine output_states(st, gr, geo, hm, dir, outp)
  type(states_t),         intent(inout) :: st
  type(grid_t),           intent(inout) :: gr
  type(geometry_t),       intent(in)    :: geo
  type(hamiltonian_t),    intent(in)    :: hm
  character(len=*),       intent(in)    :: dir
  type(output_t),         intent(in)    :: outp

  integer :: ik, ist, idim, idir, is, ierr, ip, nspin
  character(len=MAX_PATH_LEN) :: fname
  type(unit_t) :: fn_unit
  FLOAT, allocatable :: dtmp(:), elf(:,:), polarization(:, :)
  CMPLX, allocatable :: ztmp(:)

  PUSH_SUB(output_states)

  if(bitand(outp%what, OPTION__OUTPUT__DENSITY) /= 0) then
    fn_unit = units_out%length**(-gr%mesh%sb%dim)
    do is = 1, st%d%nspin
      if(st%d%nspin == 1) then
        write(fname, '(a)') 'density'
      else
        write(fname, '(a,i1)') 'density-sp', is
      end if
      call dio_function_output(outp%how, dir, fname, gr%fine%mesh, &
        st%rho(:, is), fn_unit, ierr, geo = geo, grp = st%dom_st_kpt_mpi_grp)
    end do
  end if

  if(bitand(outp%what, OPTION__OUTPUT__POL_DENSITY) /= 0) then
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

  if(bitand(outp%what, OPTION__OUTPUT__WFS) /= 0) then
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
            end if
          end do
        end do
      end if
    end do

    SAFE_DEALLOCATE_A(dtmp)
    SAFE_DEALLOCATE_A(ztmp)
  end if

  if(bitand(outp%what, OPTION__OUTPUT__WFS_SQMOD) /= 0) then
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

  POP_SUB(output_states)
end subroutine output_states

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
