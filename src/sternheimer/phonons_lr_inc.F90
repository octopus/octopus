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
!! $Id$

subroutine X(phonons_lr_infrared)(gr, geo, st, lr, kdotp_lr, imat, iatom, idir, infrared)
  type(grid_t),         intent(in)    :: gr
  type(geometry_t),     intent(in)    :: geo
  type(states_t),       intent(in)    :: st
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
    do jdir = 1, gr%sb%periodic_dim
      infrared(imat, jdir) = M_ZERO
      do ik = 1, st%d%nik
        term = M_ZERO
        do ist = 1, st%nst
          term = term + &
            TOFLOAT(X(mf_dotp)(gr%mesh, st%d%dim, lr%X(dl_psi)(:, :, ist, ik), kdotp_lr(jdir)%X(dl_psi)(:, :, ist, ik)))
        enddo
        infrared(imat, jdir) = infrared(imat, jdir) + M_TWO * term * st%smear%el_per_state * st%d%kweights(ik)
      enddo
    enddo
  endif
  
  do jdir = gr%sb%periodic_dim + 1, gr%sb%dim
    infrared(imat, jdir) = dmf_dotp(gr%mesh, gr%mesh%x(:, jdir), TOFLOAT(lr%X(dl_rho)(:, 1)))
  end do
  infrared(imat, idir) = infrared(imat, idir) - species_zval(geo%atom(iatom)%spec)
  
  POP_SUB(X(phonons_lr_infrared))
end subroutine X(phonons_lr_infrared)

! ---------------------------------------------------------
!> calculate the wavefunction associated with each normal mode
subroutine X(phonons_lr_wavefunctions)(lr, st, gr, vib, restart_load, restart_dump)
  type(lr_t),         intent(inout) :: lr
  type(states_t),     intent(inout) :: st !< not changed, just because of restart_read intent
  type(grid_t),       intent(in)    :: gr
  type(vibrations_t), intent(in)    :: vib
  type(restart_t),    intent(inout) :: restart_load
  type(restart_t),    intent(inout) :: restart_dump

  type(lr_t) :: lrtmp
  integer :: ik, ist, idim, inm, iatom, imat, ierr, idir

  PUSH_SUB(X(phonons_lr_wavefunctions))

  call lr_init(lrtmp)
  call lr_allocate(lrtmp, st, gr%mesh)

  lr%X(dl_psi) = M_ZERO

  do inm = 1, vib%num_modes

    do iatom = 1, vib%natoms
      do idir = 1, vib%ndim

        imat = vibrations_get_index(vib, iatom, idir)

        call restart_cd(restart_load, dirname=wfs_tag_sigma(phn_wfs_tag(iatom, idir), 1))
        call states_load(restart_load, st, gr, ierr, lr = lrtmp)
        call restart_cd(restart_load)

        if(ierr /= 0) then
          message(1) = "Failed to load response wavefunctions from '"//trim(wfs_tag_sigma(phn_wfs_tag(iatom, idir), 1))//"'"
          call messages_fatal(1)
        end if
            
        do ik = 1, st%d%nik
          do ist = st%st_start, st%st_end
            do idim = 1, st%d%dim

              call lalg_axpy(gr%mesh%np, vib%normal_mode(imat, inm), &
                lrtmp%X(dl_psi)(:, idim, ist, ik), lr%X(dl_psi)(:, idim, ist, ik))
                  
            end do
          end do
        end do

      end do
    end do

    call restart_cd(restart_dump, dirname=phn_nm_wfs_tag(inm))
    call states_dump(restart_dump, st, gr, ierr, lr = lr)
    call restart_cd(restart_dump)

  end do

  call lr_dealloc(lrtmp)
  POP_SUB(phonons_lr_wavefunctions)
end subroutine X(phonons_lr_wavefunctions)
