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
  ndiag = bgw%vxc_diag_nmax - bgw%vxc_diag_nmin + 1
  SAFE_ALLOCATE(psi(gr%mesh%np, 1))

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

  if(bgw%calc_exchange) then
    if(mpi_grp_is_root(mpi_world)) iunit_x = io_open(trim(dir) // 'x.dat', action='write')
    SAFE_ALLOCATE(xpsi(gr%mesh%np, 1))
    if(.not. associated(hm%hf_st)) hm%hf_st => st
    SAFE_ALLOCATE(mtxel_x(ndiag + noffdiag, st%d%nspin))
  endif

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
        call states_get_state(st, gr%mesh, 1, diag(idiag), ikk, psi(:, 1))
        mtxel(idiag, ispin) = X(mf_dotp)(gr%mesh, psi(:, 1), psi(:, 1) * vxc(:, ispin))
        if(bgw%calc_exchange) then
          !        call X(derivatives_set_bc)(gr%der, psi(:, 1))
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
  SAFE_ALLOCATE(zfield_r(gr%mesh%np))
#endif

  do is = 1, nspin
#ifdef R_TREAL
    zfield_r(1:gr%mesh%np) = cmplx(field_r(1:gr%mesh%np, is), M_ZERO, REAL_PRECISION)
#else
    zfield_r => field_r(:, is)
#endif
    call zmesh_to_cube(gr%mesh, zfield_r(:), cube, cf, local = .true.)
    call zcube_function_rs2fs(cube, cf)

    norm = M_ZERO
    do iz = 1, cube%rs_n_global(3)
      do iy = 1, cube%rs_n_global(2)
        do ix = 1, cube%rs_n_global(1) 
          if(is_wfn) then
            norm = norm + abs(cf%zrs(ix, iy, iz))**2
          else
            norm = norm + cf%zrs(ix, iy, iz)
          endif
        enddo
      enddo
    enddo
    norm = norm * gr%mesh%volume_element
    if(is_wfn) norm = sqrt(norm)
    if(mpi_grp_is_root(mpi_world)) then
      write(0,*) 'norm in real space = ', norm
    endif

    norm = M_ZERO
    do iz = 1, cube%fs_n_global(3)
      do iy = 1, cube%fs_n_global(2)
        do ix = 1, cube%fs_n_global(1) 
          if(is_wfn) then
            norm = norm + abs(cf%fs(ix, iy, iz))**2
          else
            norm = norm + cf%fs(ix, iy, iz)
          endif
        enddo
      enddo
    enddo
    if(is_wfn) then
      norm = sqrt(norm * gr%mesh%volume_element / product(cube%rs_n_global(1:3)))
    else
      norm = norm / product(cube%rs_n_global(1:3))
    endif
    if(mpi_grp_is_root(mpi_world)) then
      write(0,*) 'norm in Fourier space = ', norm
    endif
    
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
        field_g(ig, is) = cf%fs(ix, iy, iz) / product(cube%rs_n_global(1:3))
        norm = norm + field_g(ig,is)
      endif
    enddo

    ! renormalize

    write(0,*) 'shell norm = ', norm
    if(is_wfn) then
      field_g(:,:) = field_g(:,:) / sqrt(norm)
      if(abs(norm - M_ONE) > 0.01) then
        write(message(1), '(a,f12.6)') 'Wavefunction norm within G-sphere (before renormalization) is only ', norm
        call messages_warning(1)
      endif
    endif

    if(.not. is_wfn) then
      write(0,*) 'shell%red_gvec(1) = ', shell%red_gvec(1:3, 1)
      write(0,*) 'average = ', cf%fs(1,1,1) / product(cube%rs_n_global(1:3))
    endif
  enddo

#ifdef R_TREAL
  SAFE_DEALLOCATE_P(zfield_r)
#endif

  ! do Gram-Schmidt here if appropriate



  if(mpi_grp_is_root(mpi_world)) then
    call write_binary_complex_data(iunit, shell%ngvectors, shell%ngvectors, nspin, field_g)
  endif

  POP_SUB(X(bgw_write_fs))
end subroutine X(bgw_write_fs)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
