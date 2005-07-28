!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id$


subroutine X(project)(m, p, psi, ppsi, periodic, ik)
    type(mesh_type), intent(in) :: m
    type(projector), intent(in) :: p
    R_TYPE, intent(in)          :: psi(:)
    R_TYPE, intent(inout)       :: ppsi(:)
    logical, intent(in)         :: periodic
    integer, intent(in)         :: ik

    integer :: i, j, n, idim
    R_TYPE, allocatable :: lpsi(:), plpsi(:)
    R_TYPE :: uvpsi

    call push_sub('project')

    n = p%n
    allocate(lpsi(n), plpsi(n))

      lpsi(1:n)  = psi(p%jxyz(1:n))
      if(periodic) lpsi(1:n)  = lpsi(1:n) * p%phases(1:n, ik)
      plpsi(1:n) = R_TOTYPE(M_ZERO)

      do i = 1, p%c
       do j = 1, p%c
          uvpsi = sum(lpsi(1:n)*p%b(1:n, j)*m%vol_pp(p%jxyz(1:n))) 
          if(periodic) then
             plpsi(1:n) = plpsi(1:n) + p%uvu(i, j) * uvpsi * p%a(1:n, i) * p%phases(1:n, ik)
          else
             plpsi(1:n) = plpsi(1:n) + p%uvu(i, j) * uvpsi * p%a(1:n, i)
          endif
       enddo
      enddo

      ppsi(p%jxyz(1:n)) = ppsi(p%jxyz(1:n)) + plpsi(1:n)

    deallocate(lpsi)
    call pop_sub()
end subroutine X(project)

subroutine X(epot_forces) (gr, ep, st, t, reduce_)
  type(grid_type), target, intent(in)    :: gr
  type(epot_type),     intent(in)    :: ep
  type(states_type),   intent(in)    :: st
  FLOAT,     optional, intent(in)    :: t
  logical,   optional, intent(in)    :: reduce_

  type(geometry_type), pointer :: geo
  integer :: i, j, l, m, idim, ist, ik, ii, jj, ivnl
  FLOAT :: d, r, zi, zj, x(3)
  R_TYPE :: uVpsi, p
  type(atom_type), pointer :: atm
  R_TYPE, allocatable :: ppsi(:, :) 

#if defined(HAVE_MPI)
  FLOAT :: f(3)
  integer :: ierr
#endif

  call push_sub('epot_forces')
  
  geo => gr%geo

  ! init to 0
  do i = 1, geo%natoms
    geo%atom(i)%f = M_ZERO
  end do

  ! non-local component of the potential.
  allocate(ppsi(gr%m%np, st%d%dim))
  atm_loop: do i = 1, geo%natoms
     atm => geo%atom(i)
     if(atm%spec%local) cycle
     do ivnl = 1, ep%nvnl
        if(ep%p(ivnl)%index .ne. i) cycle
        ik_loop: do ik = 1, st%d%nik
          st_loop: do ist = st%st_start, st%st_end

             do j = 1, NDIM
                ppsi = R_TOTYPE(M_ZERO)
                do idim = 1, st%d%dim
                   call X(project)(gr%m, ep%dp(j, ivnl), st%X(psi)(:, idim, ist, ik), ppsi(:, idim), &
                                   periodic = .false., ik = ik)
                enddo
                atm%f(j) = atm%f(j) + M_TWO * st%occ(ist, ik) * &
                           X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, ist, ik), ppsi)
             enddo

          end do st_loop
        end do ik_loop
     enddo
  end do atm_loop
  deallocate(ppsi)

#if defined(HAVE_MPI)
  do i = 1, geo%natoms
    atm => geo%atom(i)
    if(atm%spec%local) cycle
    if(present(reduce_)) then
      if(reduce_) then
        call MPI_ALLREDUCE(atm%f(1), f(1), NDIM, &
             MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
        atm%f = f
      end if
    end if
  enddo
#endif
  
  ! Now the ion, ion force term
  do i = 1, geo%natoms
    zi = geo%atom(i)%spec%Z_val
    do j = 1, geo%natoms
      if(i .ne. j) then
        zj = geo%atom(j)%spec%Z_val
        r = sqrt(sum((geo%atom(i)%x(1:NDIM) - geo%atom(j)%x(1:NDIM))**2))
        d = zi * zj/r**3
        
        geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + &
             d*(geo%atom(i)%x(1:NDIM) - geo%atom(j)%x(1:NDIM))
      end if
    end do
  end do
  
  ! now comes the local part of the PP
  if(.not.simul_box_is_periodic(gr%sb).or.geo%only_user_def) then ! Real space
    call local_RS()
#ifdef HAVE_FFT
  else ! Fourier space
    call local_FS()
#endif
  end if

  if(present(t).and.ep%no_lasers>0) then
    call laser_field(gr%sb, ep%no_lasers, ep%lasers, t, x)
    do i = 1, geo%natoms
      geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + &
           geo%atom(i)%spec%Z_val * x(1:NDIM)
    end do
  end if

  if(associated(ep%e)) then
    do i = 1, geo%natoms
      geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + &
           geo%atom(i)%spec%Z_val * ep%e(1:NDIM)
    end do
  end if

  !TODO: forces due to the magnetic fields (static and time-dependent)

  call pop_sub()
  
contains
  subroutine local_RS()
    FLOAT :: r, x(3), d, gv(3)
    integer  :: ns
    
    ns = min(2, st%d%nspin)
    
    do i = 1, geo%natoms
      atm => geo%atom(i)
      do j = 1, NP
        call mesh_r(gr%m, j, r, x=x, a=atm%x)
        if(r < r_small) cycle
        
        call specie_get_glocal(atm%spec, x, gv)
        d = sum(st%rho(j, 1:ns))*gr%m%vol_pp(j)
        atm%f(1:NDIM) = atm%f(1:NDIM) - d*gv(1:NDIM)
      end do
    end do
  end subroutine local_RS
  
#ifdef HAVE_FFT
  subroutine local_FS()
    type(dcf) :: cf_for
    FLOAT, allocatable :: force(:)
    
    allocate(force(NP))
    call dcf_new_from(cf_for, ep%local_cf(1)) ! at least one specie must exist
    call dcf_alloc_FS(cf_for)
    call dcf_alloc_RS(cf_for)
    
    do i = 1, geo%natoms
      atm => geo%atom(i)
      do j = 1, NDIM
        cf_for%FS = M_z0
        call cf_phase_factor(gr%sb, gr%m, atm%x, ep%local_cf(atm%spec%index), cf_for)
        
        call dcf_FS_grad(gr%sb, gr%m, cf_for, j)
        call dcf_FS2RS(cf_for)
        call dcf2mf(gr%m, cf_for, force)
        do l = 1, st%d%nspin
          atm%f(j) = atm%f(j) + sum(force(1:NP)*st%rho(1:NP, l)*gr%m%vol_pp(1:NP))
        end do
      end do
    end do
    
    call dcf_free(cf_for)
    deallocate(force)
  end subroutine local_FS
#endif

end subroutine X(epot_forces)
