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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

subroutine X(bgw_vxc_dat)(bgw, dir, st, gr, hm, vxc)
  type(output_bgw_t),     intent(in)    :: bgw
  character(len=*),       intent(in)    :: dir
  type(states_t), target, intent(in)    :: st
  type(grid_t),           intent(in)    :: gr
  type(hamiltonian_t),    intent(inout) :: hm
  FLOAT,                  intent(in)    :: vxc(:,:)

  integer :: iunit, iunit_x, ispin, ik, ikk, ist, ist2, idiag, ioff, ndiag, noffdiag, spin_index(st%d%nspin)
  integer, allocatable :: diag(:), off1(:), off2(:)
  FLOAT :: kpoint(3)
  R_TYPE, allocatable :: psi(:,:), psi2(:), xpsi(:,:)
  CMPLX, allocatable :: mtxel(:,:), mtxel_x(:,:)

  PUSH_SUB(X(bgw_vxc_dat))

#ifdef HAVE_BERKELEYGW

  if(st%parallel_in_states) call messages_not_implemented("BerkeleyGW output parallel in states")
  if(st%d%kpt%parallel) call messages_not_implemented("BerkeleyGW output parallel in k-points")

  if(mpi_grp_is_root(mpi_world)) iunit = io_open(trim(dir) // 'vxc.dat', action='write')
  SAFE_ALLOCATE(psi(1:gr%mesh%np, 1))

  ndiag = bgw%vxc_diag_nmax - bgw%vxc_diag_nmin + 1
  if(bgw%vxc_diag_nmin < 1 .or. bgw%vxc_diag_nmax < 1) then
    ndiag = 0
  else
    ndiag = bgw%vxc_diag_nmax - bgw%vxc_diag_nmin + 1
  endif

  if(bgw%vxc_offdiag_nmin < 1 .or. bgw%vxc_offdiag_nmax < 1) then
    noffdiag = 0
  else
    noffdiag = (bgw%vxc_offdiag_nmax - bgw%vxc_offdiag_nmin + 1)**2
  endif
  if(noffdiag > 0) then
    SAFE_ALLOCATE(psi2(1:gr%mesh%np))
    SAFE_ALLOCATE(off1(1:noffdiag))
    SAFE_ALLOCATE(off2(1:noffdiag))
  endif
  SAFE_ALLOCATE(mtxel(1:ndiag + noffdiag, 1:st%d%nspin))

  if(bgw%calc_exchange) then
    if(mpi_grp_is_root(mpi_world)) iunit_x = io_open(trim(dir) // 'x.dat', action='write')
    SAFE_ALLOCATE(xpsi(1:gr%mesh%np, 1))
    if(.not. associated(hm%hf_st)) hm%hf_st => st
    SAFE_ALLOCATE(mtxel_x(1:ndiag + noffdiag, 1:st%d%nspin))
  endif

  ! BerkeleyGW allows using only spin down, but we will not give that option here
  do ispin = 1, st%d%nspin
    spin_index(ispin) = ispin
  enddo

  SAFE_ALLOCATE(diag(1:ndiag))
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
        call states_get_state(st, gr%mesh, 1, diag(idiag), ikk, psi(:, 1))
        ! multiplying psi*vxc first might be more efficient
        mtxel(idiag, ispin) = X(mf_dotp)(gr%mesh, psi(:, 1), psi(:, 1) * vxc(:, ispin))
        if(bgw%calc_exchange) then
          xpsi(:,:) = M_ZERO
          call X(exchange_operator)(hm, gr%der, psi, xpsi, ist, ikk, M_ONE)
          mtxel_x(idiag, ispin) = X(mf_dotp)(gr%mesh, psi(:, 1), xpsi(:, 1))
        endif
      enddo

      ! could do only upper or lower triangle here
      do ioff = 1, noffdiag
        call states_get_state(st, gr%mesh, 1, off1(ioff), ikk, psi(:, 1))
        call states_get_state(st, gr%mesh, 1, off2(ioff), ikk, psi2)
        mtxel(ndiag + ioff, ispin) = X(mf_dotp)(gr%mesh, psi(:, 1), psi2(:) * vxc(:, ispin))
        ! FIXME: we should calc xpsi only for each state, not for each offdiag
        if(bgw%calc_exchange) then
          xpsi(:,:) = M_ZERO
          call X(exchange_operator)(hm, gr%der, psi, xpsi, ist, ikk, M_ONE)
          mtxel_x(ndiag + ioff, ispin) = R_CONJ(X(mf_dotp)(gr%mesh, psi2, xpsi(:, 1)))
        endif
      enddo
    enddo

    ! convert to eV
    mtxel(:,:) = M_TWO * P_Ry * mtxel(:,:)
    if(mpi_grp_is_root(mpi_world)) &
      call write_matrix_elements(iunit, kpoint, st%d%nspin, ndiag, noffdiag, spin_index, diag, off1, off2, mtxel)

    if(bgw%calc_exchange) then
      mtxel_x(:,:) = M_TWO * P_Ry * mtxel_x(:,:)
      if(mpi_grp_is_root(mpi_world)) &
        call write_matrix_elements(iunit_x, kpoint, st%d%nspin, ndiag, noffdiag, spin_index, diag, off1, off2, mtxel_x)
    endif
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

  if(bgw%calc_exchange) then
    if(mpi_grp_is_root(mpi_world)) call io_close(iunit_x)
    SAFE_DEALLOCATE_A(xpsi)
    SAFE_DEALLOCATE_A(mtxel_x)
  endif

#else
    message(1) = "Cannot do BerkeleyGW output: the library was not linked."
    call messages_fatal(1)
#endif

  POP_SUB(X(bgw_vxc_dat))

end subroutine X(bgw_vxc_dat)

! --------------------------------------------------------- 
!> Calculate 'vmtxel' file of dipole matrix elements for BerkeleyGW BSE
subroutine X(bgw_vmtxel)(bgw, dir, st, gr, ifmax)
  type(output_bgw_t), intent(in) :: bgw
  character(len=*),   intent(in) :: dir
  type(states_t),     intent(in) :: st
  type(grid_t),       intent(in) :: gr
  integer,            intent(in) :: ifmax(:,:)

  integer :: iunit, nmat, ik, ikk, ikcvs, is, ic, iv, ip
  R_TYPE, allocatable :: psi(:), rpsi(:)
  CMPLX, allocatable :: vmtxel(:) ! could be real
  FLOAT, allocatable :: rvec(:)

  PUSH_SUB(X(bgw_vmtxel))

  nmat = st%d%nik * bgw%vmtxel_ncband * bgw%vmtxel_nvband * st%d%nspin
  SAFE_ALLOCATE(psi(1:gr%mesh%np))
  SAFE_ALLOCATE(rpsi(1:gr%mesh%np))
  SAFE_ALLOCATE(vmtxel(1:nmat))
  SAFE_ALLOCATE(rvec(1:gr%mesh%np))

  forall(ip = 1:gr%mesh%np) &
    rvec(ip) = dot_product(gr%mesh%x(ip, 1:3), bgw%vmtxel_polarization(1:3))

  do ik = st%d%kpt%start, st%d%kpt%end, st%d%nspin
    do is = 1, st%d%nspin
      ikk = ik + is - 1
      do iv = 1, bgw%vmtxel_nvband
        call states_get_state(st, gr%mesh, 1, ifmax(ik, is) - iv + 1, ikk, psi(:))
        rpsi(:) = psi(:) * rvec(:)
        do ic = 1, bgw%vmtxel_ncband
          call states_get_state(st, gr%mesh, 1, ifmax(ik, is) + ic, ikk, psi(:))
          ikcvs = is + (iv - 1 + (ic - 1 + (ik - 1)*bgw%vmtxel_ncband)*bgw%vmtxel_nvband)*st%d%nspin
          vmtxel(ikcvs) = X(mf_dotp)(gr%mesh, psi(:), rpsi(:))
!          write(6,*) ikcvs, ikk, is, ifmax(ik, is) - iv + 1, ifmax(ik, is) + ic, vmtxel(ikcvs)
          ! NB: Casida eps_diff file has these values times sqrt(nspin)
        enddo
      enddo
    enddo
  enddo

  if(mpi_grp_is_root(mpi_world)) then
    iunit = io_open(trim(dir) // 'vmtxel.dat', action='write', form='formatted')
    write(iunit,*) st%d%nik/st%d%nspin,bgw%vmtxel_ncband,bgw%vmtxel_nvband,st%d%nspin,1
    write(iunit,*) (vmtxel(ikcvs),ikcvs=1,nmat)
    call io_close(iunit)
  endif

  if(mpi_grp_is_root(mpi_world)) then
    iunit = io_open(trim(dir) // 'vmtxel', action='write', form='unformatted')
    write(iunit) st%d%nik/st%d%nspin,bgw%vmtxel_ncband,bgw%vmtxel_nvband,st%d%nspin,1
    write(iunit) (vmtxel(ikcvs),ikcvs=1,nmat)
    call io_close(iunit)
  endif

  SAFE_DEALLOCATE_A(vmtxel)
  SAFE_DEALLOCATE_A(rpsi)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(rvec)

  POP_SUB(X(bgw_vmtxel))
  
end subroutine X(bgw_vmtxel)

! --------------------------------------------------------- 
subroutine X(bgw_write_fs)(iunit, field_r, field_g, shell, nspin, gr, cube, cf, is_wfn)
  integer,               intent(in)    :: iunit
  R_TYPE, target,        intent(in)    :: field_r(:,:)
  CMPLX,                 intent(inout) :: field_g(:,:)
  type(fourier_shell_t), intent(in)    :: shell
  integer,               intent(in)    :: nspin
  type(grid_t),          intent(in)    :: gr
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf
  logical,               intent(in)    :: is_wfn !< make false for RHO, VXC

  integer :: ig, ix, iy, iz, is
  FLOAT :: norm
  CMPLX, pointer :: zfield_r(:)

  PUSH_SUB(X(bgw_write_fs))

  ! We always need to use FFT`s from a complex function, since BerkeleyGW does not use
  ! the half-sphere Hermitian representation for real functions.

#ifdef R_TREAL
  SAFE_ALLOCATE(zfield_r(1:gr%mesh%np))
#endif

  do is = 1, nspin
#ifdef R_TREAL
    zfield_r(1:gr%mesh%np) = TOCMPLX(field_r(1:gr%mesh%np, is), M_ZERO)
#else
    zfield_r => field_r(:, is)
#endif
    call zmesh_to_cube(gr%mesh, zfield_r(:), cube, cf, local = .true.)
    call zcube_function_rs2fs(cube, cf)

!    if(is_wfn) then
!      ! norm in real space
!      norm = M_ZERO
!      do iz = 1, cube%rs_n_global(3)
!        do iy = 1, cube%rs_n_global(2)
!          do ix = 1, cube%rs_n_global(1) 
!            norm = norm + abs(cf%zrs(ix, iy, iz))**2
!          enddo
!        enddo
!      enddo
!      norm = sqrt(norm * gr%mesh%volume_element)
!
!      norm = M_ZERO
!      do iz = 1, cube%fs_n_global(3)
!        do iy = 1, cube%fs_n_global(2)
!          do ix = 1, cube%fs_n_global(1) 
!            norm = norm + abs(cf%fs(ix, iy, iz))**2
!          enddo
!        enddo
!      enddo
!      norm = sqrt(norm * gr%mesh%volume_element / product(cube%rs_n_global(1:3)))
!    endif
!    
    field_g(:,:) = M_ZERO
    norm = M_ZERO
    do ig = 1, shell%ngvectors
      ix = shell%coords(1, ig)
      iy = shell%coords(2, ig)
      iz = shell%coords(3, ig)
      if(is_wfn) then
        field_g(ig, is) = cf%fs(ix, iy, iz) * &
          sqrt(gr%mesh%volume_element / product(cube%rs_n_global(1:3)))
        norm = norm + abs(field_g(ig,is))**2
      else
        field_g(ig, is) = cf%fs(ix, iy, iz) * gr%mesh%volume_element
      endif
    enddo

    ! renormalize
    if(is_wfn) then
      field_g(:,:) = field_g(:,:) / sqrt(norm)
      if(abs(norm - M_ONE) > 0.01) then
        write(message(1), '(a,f12.6)') 'Wavefunction norm within G-sphere (before renormalization) is only ', norm
        call messages_warning(1)
      endif
    endif

  enddo

#ifdef R_TREAL
  SAFE_DEALLOCATE_P(zfield_r)
#endif

  ! do Gram-Schmidt here if appropriate



  if(mpi_grp_is_root(mpi_world)) then
    call write_binary_complex_data(iunit, shell%ngvectors, ubound(field_g, 1), nspin, field_g)
  endif

  POP_SUB(X(bgw_write_fs))
end subroutine X(bgw_write_fs)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
