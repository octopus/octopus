!! Copyright (C) 2011 D. Strubbe
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
!! $Id: output_etsf_inc.F90 5880 2009-09-03 23:44:44Z dstrubbe $

subroutine X(bgw_vxc_dat)(bgw, dir, st, gr, vxc)
  type(output_bgw_t), intent(in) :: bgw
  character(len=*),   intent(in) :: dir
  type(states_t),     intent(in) :: st
  type(grid_t),       intent(in) :: gr
  FLOAT,              intent(in) :: vxc(:,:)

  integer :: iunit, ispin, ik, ikk, ist, ist2, idiag, ioff, ndiag, noffdiag, spin_index(st%d%nspin)
  integer, allocatable :: diag(:), off1(:), off2(:)
  FLOAT :: kpoint(3)
  R_TYPE, allocatable :: psi(:), psi2(:)
  CMPLX, allocatable :: mtxel(:,:)

  PUSH_SUB(X(bgw_vxc_dat))

#ifdef HAVE_BERKELEYGW

  if(st%parallel_in_states) call messages_not_implemented("BerkeleyGW output parallel in states")
  if(st%d%kpt%parallel) call messages_not_implemented("BerkeleyGW output parallel in k-points")

  if(mpi_grp_is_root(mpi_world)) iunit = io_open(trim(dir) // 'vxc.dat', action='write')
  ndiag = bgw%vxc_diag_nmax - bgw%vxc_diag_nmin + 1
  SAFE_ALLOCATE(psi(gr%mesh%np))
  if(bgw%vxc_offdiag_nmin < 1 .or. bgw%vxc_offdiag_nmax < 1) then
    noffdiag = 0
  else
    noffdiag = (bgw%vxc_offdiag_nmax - bgw%vxc_offdiag_nmin + 1)**2
  endif
  if(noffdiag > 0) then
    SAFE_ALLOCATE(psi2(gr%mesh%np))
    SAFE_ALLOCATE(off1(noffdiag))
    SAFE_ALLOCATE(off2(noffdiag))
  endif
  SAFE_ALLOCATE(mtxel(ndiag + noffdiag, st%d%nspin))

  message(1) = "BerkeleyGW output: vxc.dat"
  call messages_info(1)

  ! BerkeleyGW allows using only spin down, but we will not give that option here
  do ispin = 1, st%d%nspin
    spin_index(ispin) = ispin
  enddo

  SAFE_ALLOCATE(diag(ndiag))
  do idiag = 1, ndiag
    diag(idiag) = bgw%vxc_diag_nmin + idiag - 1
  enddo

  if(noffdiag > 0) then
    ioff = 1
    do ist = bgw%vxc_offdiag_nmin, bgw%vxc_offdiag_nmax
      do ist2 = bgw%vxc_offdiag_nmin, bgw%vxc_offdiag_nmax
        off1(ioff) = ist
        off2(ioff) = ist2
        ioff = ioff + 1
      enddo
    enddo
  endif

  ! in case of hybrids, we should apply exchange operator too here
  ! in that case, we can write x.dat file as well
  
  do ik = st%d%kpt%start, st%d%kpt%end, st%d%nspin
    kpoint(1:gr%sb%dim) = gr%sb%kpoints%reduced%red_point(1:gr%sb%dim, ik) ! crystal coordinates

    do ispin = 1, st%d%nspin
      ikk = ik + ispin - 1
      do idiag = 1, ndiag
        call states_get_state(st, gr%mesh, 1, diag(idiag), ikk, psi)
        mtxel(idiag, ispin) = X(mf_dotp)(gr%mesh, psi(:), psi(:) * vxc(:, ispin))
      enddo

      ! could do only upper or lower triangle here
      do ioff = 1, noffdiag
        call states_get_state(st, gr%mesh, 1, off1(ioff), ikk, psi)
        call states_get_state(st, gr%mesh, 1, off2(ioff), ikk, psi2)
        mtxel(ndiag + ioff, ispin) = X(mf_dotp)(gr%mesh, psi(:), psi2(:) * vxc(:, ispin))
      enddo
    enddo

    ! convert to eV
    mtxel(:,:) = M_TWO * P_Ry * mtxel(:,:)

    if(mpi_grp_is_root(mpi_world)) &
      call write_matrix_elements(iunit, kpoint, st%d%nspin, ndiag, noffdiag, spin_index, diag, off1, off2, mtxel) 
  enddo

  if(mpi_grp_is_root(mpi_world)) call io_close(iunit)
  SAFE_DEALLOCATE_A(diag)
  SAFE_DEALLOCATE_A(psi)
  if(noffdiag > 0) then
    SAFE_DEALLOCATE_A(off1)
    SAFE_DEALLOCATE_A(off2)
    SAFE_DEALLOCATE_A(psi2)
  endif
  SAFE_DEALLOCATE_A(mtxel)

#else
    message(1) = "Cannot do BerkeleyGW output: the library was not linked."
    call messages_fatal(1)
#endif

  POP_SUB(X(bgw_vxc_dat))

end subroutine X(bgw_vxc_dat)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
