!! Copyright (C) 2007-2012 Xavier Andrade, David Strubbe
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

subroutine X(phonons_lr_infrared)(mesh, ions, st, lr, kdotp_lr, imat, iatom, idir, infrared)
  type(mesh_t),         intent(in)    :: mesh
  type(ions_t),         intent(in)    :: ions
  type(states_elec_t),  intent(in)    :: st
  type(lr_t),           intent(in)    :: lr
  type(lr_t),           intent(in)    :: kdotp_lr(:) !< (ndim)
  integer,              intent(in)    :: imat
  integer,              intent(in)    :: iatom
  integer,              intent(in)    :: idir
  FLOAT,                intent(inout) :: infrared(:,:) !< (nmat, nmat)

  integer :: jdir, ik, ist
  FLOAT :: term

  PUSH_SUB(X(phonons_lr_infrared))

  if(smear_is_semiconducting(st%smear)) then
    do jdir = 1, ions%space%periodic_dim
      infrared(imat, jdir) = M_ZERO
      do ik = 1, st%d%nik
        term = M_ZERO
        do ist = 1, st%nst
          term = term + &
            TOFLOAT(X(mf_dotp)(mesh, st%d%dim, lr%X(dl_psi)(:, :, ist, ik), kdotp_lr(jdir)%X(dl_psi)(:, :, ist, ik)))
        end do
        infrared(imat, jdir) = infrared(imat, jdir) + M_TWO * term * st%smear%el_per_state * st%d%kweights(ik)
      end do
    end do
  end if
  
  do jdir = ions%space%periodic_dim + 1, ions%space%dim
    infrared(imat, jdir) = dmf_dotp(mesh, mesh%x(:, jdir), TOFLOAT(lr%X(dl_rho)(:, 1)))
  end do
  infrared(imat, idir) = infrared(imat, idir) - species_zval(ions%atom(iatom)%species)
  
  POP_SUB(X(phonons_lr_infrared))
end subroutine X(phonons_lr_infrared)

! ---------------------------------------------------------
!> calculate the wavefunction associated with each normal mode
subroutine X(phonons_lr_wavefunctions)(lr, namespace, space, st, mesh, kpoints, vib, restart_load, restart_dump)
  type(lr_t),         intent(inout) :: lr
  type(namespace_t),  intent(in)    :: namespace
  type(space_t),      intent(in)    :: space
  type(states_elec_t),intent(inout) :: st !< not changed, just because of restart_read intent
  type(mesh_t),       intent(in)    :: mesh
  type(kpoints_t),    intent(in)    :: kpoints
  type(vibrations_t), intent(in)    :: vib
  type(restart_t),    intent(inout) :: restart_load
  type(restart_t),    intent(inout) :: restart_dump

  type(lr_t) :: lrtmp
  integer :: ik, ist, idim, inm, iatom, imat, ierr, idir

  PUSH_SUB(X(phonons_lr_wavefunctions))

  call lr_init(lrtmp)
  call lr_allocate(lrtmp, st, mesh)

  lr%X(dl_psi) = M_ZERO

  do inm = 1, vib%num_modes

    do iatom = 1, vib%natoms
      do idir = 1, vib%ndim

        imat = vibrations_get_index(vib, iatom, idir)

        call restart_open_dir(restart_load, wfs_tag_sigma(phn_wfs_tag(iatom, idir), 1), ierr)
        if (ierr == 0) then
          call states_elec_load(restart_load, namespace, space, st, mesh, kpoints, ierr, lr = lrtmp)
        end if
        call restart_close_dir(restart_load)

        if(ierr /= 0) then
          message(1) = "Unable to read response wavefunctions from '"//trim(wfs_tag_sigma(phn_wfs_tag(iatom, idir), 1))//"'."
          call messages_fatal(1)
        end if
            
        do ik = 1, st%d%nik
          do ist = st%st_start, st%st_end
            do idim = 1, st%d%dim

              call lalg_axpy(mesh%np, vib%normal_mode(imat, inm), &
                lrtmp%X(dl_psi)(:, idim, ist, ik), lr%X(dl_psi)(:, idim, ist, ik))
                  
            end do
          end do
        end do

      end do
    end do

    call restart_open_dir(restart_dump, phn_nm_wfs_tag(inm), ierr)
    if (ierr == 0) then
      call states_elec_dump(restart_dump, space, st, mesh, kpoints, ierr, lr = lr)
    end if
    if (ierr /= 0) then
      message(1) = "Unable to write response wavefunctions to '"//trim(phn_nm_wfs_tag(inm))//"'."
      call messages_warning(1)
    end if
    call restart_close_dir(restart_dump)

  end do

  call lr_dealloc(lrtmp)
  POP_SUB(phonons_lr_wavefunctions)
end subroutine X(phonons_lr_wavefunctions)
