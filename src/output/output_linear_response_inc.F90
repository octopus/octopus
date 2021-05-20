!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade.
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
subroutine X(output_lr) (outp, namespace, space, dir, st, mesh, lr, idir, isigma, ions, pert_unit)
  type(output_t),       intent(in) :: outp
  type(namespace_t),    intent(in) :: namespace
  type(space_t),        intent(in) :: space
  character(len=*),     intent(in) :: dir
  type(states_elec_t),  intent(in) :: st
  type(mesh_t),         intent(in) :: mesh
  type(lr_t),           intent(in) :: lr
  integer,              intent(in) :: idir      !< direction of perturbation
  integer,              intent(in) :: isigma
  type(ions_t),         intent(in) :: ions
  type(unit_t),         intent(in) :: pert_unit !< unit for perturbation

  integer :: ik, ist, idim, ierr, is, idir2
  character(len=80) :: fname
  type(unit_t) :: fn_unit
  FLOAT, allocatable :: dtmp(:)
  R_TYPE, allocatable :: tmp(:)
  character :: sigma

  PUSH_SUB(X(output_lr))

  if(isigma == 1) then
    sigma = '+'
  else
    sigma = '-'
  end if

  if(isigma == 1) then ! the density, current, etc. are only defined for the + frequency

    if (outp%what(OPTION__OUTPUT__DENSITY)) then
      fn_unit = units_out%length**(-space%dim)
      do is = 1, st%d%nspin
        if(st%d%nspin == 1) then
          write(fname, '(2a)') 'lr_density-', index2axis(idir)
        else
          write(fname, '(a,i1,2a)') 'lr_density-sp', is, '-', index2axis(idir)
        end if
        call X(io_function_output)(outp%how(OPTION__OUTPUT__DENSITY), dir, fname, namespace, space, &
          mesh, lr%X(dl_rho)(:, is), fn_unit / pert_unit, ierr, ions = ions)
      end do
    end if

    if (outp%what(OPTION__OUTPUT__POL_DENSITY)) then
      fn_unit = units_out%length**(1 - space%dim)
      SAFE_ALLOCATE(tmp(1:mesh%np))
      do is = 1, st%d%nspin
        do idir2 = 1, space%dim
          tmp(1:mesh%np) = -mesh%x(1:mesh%np, idir2) * lr%X(dl_rho)(:, is)
          if(st%d%nspin == 1) then
            write(fname, '(4a)') 'alpha_density-', index2axis(idir2), '-', index2axis(idir)
          else
            write(fname, '(a,i1,4a)') 'alpha_density-sp', is, '-', index2axis(idir2), '-', index2axis(idir)
          end if
          call X(io_function_output)(outp%how(OPTION__OUTPUT__POL_DENSITY), dir, fname, namespace, space, &
            mesh, tmp, fn_unit / pert_unit, ierr, ions = ions)
        end do
      end do
      SAFE_DEALLOCATE_A(tmp)
    end if

    if (outp%what(OPTION__OUTPUT__CURRENT)) then
      if (states_are_complex(st)) then
        fn_unit = units_out%time**(-1) * units_out%length**(-space%dim)
        do is = 1, st%d%nspin
          do idir2 = 1, space%dim
            if(st%d%nspin == 1) then
              write(fname, '(4a)') 'lr_current-', index2axis(idir2), '-',  index2axis(idir)
            else
              write(fname, '(a,i1,4a)') 'lr_current-sp', is, '-', index2axis(idir2), '-',  index2axis(idir)
            end if
            call zio_function_output(outp%how(OPTION__OUTPUT__CURRENT), dir, fname, namespace, space, &
              mesh, lr%dl_j(:, idir2, is), fn_unit / pert_unit, ierr, ions = ions)
          end do
        end do
      else
        message(1) = 'No current density output for real states since it is identically zero.'
        call messages_warning(1, namespace=namespace)
      end if
    end if

    if (space%dim==3) then
      if (outp%what(OPTION__OUTPUT__ELF)) call lr_elf('lr_elf_D','lr_elf')
    end if

  end if ! isigma == 1


  if (outp%what(OPTION__OUTPUT__WFS)) then
    fn_unit = sqrt(units_out%length**(-space%dim))
    do ist = st%st_start, st%st_end
      if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = st%d%kpt%start, st%d%kpt%end
          do idim = 1, st%d%dim
            if(st%d%nik > 1) then
              if(st%d%dim > 1) then
                write(fname, '(a,i3.3,a,i3.3,a,i1,3a)') &
                  'lr_wf-k', ik, '-st', ist, '-sp', idim, '-', index2axis(idir), sigma
              else
                write(fname, '(a,i3.3,a,i3.3,3a)') &
                  'lr_wf-k', ik, '-st', ist, '-', index2axis(idir), sigma
              end if
            else
              if(st%d%dim > 1) then
                write(fname, '(a,i3.3,a,i1,3a)') &
                  'lr_wf-st', ist, '-sp', idim, '-', index2axis(idir), sigma
              else
                write(fname, '(a,i3.3,3a)') &
                  'lr_wf-st', ist, '-', index2axis(idir), sigma
              end if
            end if
            call X(io_function_output) (outp%how(OPTION__OUTPUT__WFS), dir, fname, namespace, space, &
              mesh, lr%X(dl_psi) (1:, idim, ist, ik), fn_unit  / pert_unit, ierr, ions = ions)
          end do
        end do
      end if
    end do
  end if

  if (outp%what(OPTION__OUTPUT__WFS_SQMOD)) then
    fn_unit = units_out%length**(-space%dim)
    SAFE_ALLOCATE(dtmp(1:mesh%np_part))
    do ist = st%st_start, st%st_end
      if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = st%d%kpt%start, st%d%kpt%end
          do idim = 1, st%d%dim
            if(st%d%nik > 1) then
              if(st%d%dim > 1) then
                write(fname, '(a,i3.3,a,i3.3,a,i1,3a)') &
                  'sqm_lr_wf-k', ik, '-st', ist, '-sp', idim, '-', index2axis(idir), sigma
              else
                write(fname, '(a,i3.3,a,i3.3,3a)') &
                  'sqm_lr_wf-k', ik, '-st', ist, '-', index2axis(idir), sigma
              end if
            else
              if(st%d%dim > 1) then
                write(fname, '(a,i3.3,a,i1,3a)') &
                  'sqm_lr_wf-st', ist, '-sp', idim, '-', index2axis(idir), sigma
              else
                write(fname, '(a,i3.3,3a)') &
                  'sqm_lr_wf-st', ist, '-', index2axis(idir), sigma
              end if
            end if

            dtmp = abs(lr%X(dl_psi) (:, idim, ist, ik))**2
            call dio_function_output (outp%how(OPTION__OUTPUT__WFS_SQMOD), dir, fname, namespace, space, &
              mesh, dtmp, fn_unit / pert_unit, ierr, ions = ions)
          end do
        end do
      end if
    end do
    SAFE_DEALLOCATE_A(dtmp)
  end if

POP_SUB(X(output_lr))

contains

  ! ---------------------------------------------------------
  subroutine lr_elf(filename1, filename2)
    character(len=*), intent(in) :: filename1
    character(len=*), intent(in) :: filename2
    
    integer :: is, ierr

    PUSH_SUB(X(output_lr).lr_elf)

    ! these quantities are dimensionless, before the perturbation
    do is = 1, st%d%nspin
      if(st%d%nspin == 1) then
        write(fname, '(3a)') trim(filename1), '-', index2axis(idir)
      else
        write(fname, '(2a,i1,2a)') trim(filename1), '-sp', is, '-', index2axis(idir)
      end if
      call X(io_function_output)(outp%how(OPTION__OUTPUT__ELF), dir, trim(fname), namespace, space, &
        mesh, lr%X(dl_de)(1:mesh%np,is), unit_one / pert_unit, ierr, ions = ions)
    end do

    do is = 1, st%d%nspin
      if(st%d%nspin == 1) then
        write(fname, '(3a)') trim(filename2), '-', index2axis(idir)
      else
        write(fname, '(2a,i1,2a)') trim(filename2), '-sp', is, '-', index2axis(idir)
      end if
      call X(io_function_output)(outp%how(OPTION__OUTPUT__ELF), dir, trim(fname), namespace, space, &
        mesh, lr%X(dl_elf)(1:mesh%np,is), unit_one / pert_unit, ierr, ions = ions)
    end do

    POP_SUB(X(output_lr).lr_elf)

  end subroutine lr_elf

end subroutine X(output_lr)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
