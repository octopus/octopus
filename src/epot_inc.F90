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

subroutine X(epot_forces) (ep, mesh, st, geo, t, reduce_)
  type(epot_type),     intent(IN)    :: ep
  type(mesh_type),     intent(IN)    :: mesh
  type(geometry_type), intent(inout) :: geo
  type(states_type),   intent(IN)    :: st
  FLOAT,     optional, intent(in)    :: t
  logical,   optional, intent(in)    :: reduce_

  integer :: i, j, l, m, idim, ist, ik, ii, jj, ivnl
  FLOAT :: d, r, zi, zj, x(3)
  R_TYPE :: uVpsi, p
  type(atom_type), pointer :: atm

#if defined(HAVE_MPI)
  FLOAT :: f(3)
  integer :: ierr
#endif

  call push_sub('epot_forces')
  
  ! init to 0
  do i = 1, geo%natoms
    geo%atom(i)%f = M_ZERO
  end do


  ivnl = 1
  atm_loop: do i = 1, geo%natoms
    atm => geo%atom(i)
    if(atm%spec%local) cycle

    do l = 0, atm%spec%ps%l_max
      if(l == atm%spec%ps%l_loc) cycle
      do m = -l, l
        
        ik_loop: do ik = 1, st%d%nik
          st_loop: do ist = st%st_start, st%st_end
            dim_loop: do idim = 1, st%d%dim

              do ii = 1, ep%vnl(ivnl)%c
                do jj = 1, ep%vnl(ivnl)%c
                  
                  p = sum(ep%vnl(ivnl)%uv(:, ii) *  &
                     st%X(psi)(ep%vnl(ivnl)%jxyz(:), idim, ist, ik) * mesh%vol_pp(ep%vnl(ivnl)%jxyz(:))**2)
                  uvpsi =  st%occ(ist, ik) * p * ep%vnl(ivnl)%uvu(ii, jj)
                  
                  do j = 1, conf%dim
                    p = sum( ep%vnl(ivnl)%duv(j, :, jj) * &
                       R_CONJ(st%X(psi)(ep%vnl(ivnl)%jxyz(:), idim, ist, ik)) )
                    atm%f(j) = atm%f(j) + M_TWO * R_REAL(uvpsi * p)
                  end do
                end do
              end do
              
            end do dim_loop
          end do st_loop
        end do ik_loop
        
        ivnl = ivnl + 1
      end do
    end do
    
  end do atm_loop

#if defined(HAVE_MPI)
  do i = 1, geo%natoms
    atm => geo%atom(i)
    if(atm%spec%local) cycle
    if(present(reduce_)) then
      if(reduce_) then
        call MPI_ALLREDUCE(atm%f(1), f(1), conf%dim, &
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
        r = sqrt(sum((geo%atom(i)%x(1:conf%dim) - geo%atom(j)%x(1:conf%dim))**2))
        d = zi * zj/r**3
        
        geo%atom(i)%f(1:conf%dim) = geo%atom(i)%f(1:conf%dim) + &
             d*(geo%atom(i)%x(1:conf%dim) - geo%atom(j)%x(1:conf%dim))
      end if
    end do
  end do
  
  ! now comes the local part of the PP
  if(conf%periodic_dim==0.or.conf%only_user_def) then ! Real space
    call local_RS()
#ifdef HAVE_FFT
  else ! Fourier space
    call local_FS()
#endif
  end if

  if(present(t).and.ep%no_lasers>0) then
    call laser_field(ep%no_lasers, ep%lasers, t, x)
    do i = 1, geo%natoms
      geo%atom(i)%f(1:conf%dim) = geo%atom(i)%f(1:conf%dim) + &
           geo%atom(i)%spec%Z_val * x(1:conf%dim)
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
      do j = 1, mesh%np
        call mesh_r(mesh, j, r, x=x, a=atm%x)
        if(r < r_small) cycle
        
        call specie_get_glocal(atm%spec, x, gv)
        d = sum(st%rho(j, 1:ns))*mesh%vol_pp(j)
        atm%f(1:conf%dim) = atm%f(1:conf%dim) - d*gv(1:conf%dim)
      end do
    end do
  end subroutine local_RS
  
#ifdef HAVE_FFT
  subroutine local_FS()
    type(dcf) :: cf_for
    FLOAT, allocatable :: force(:)
    
    allocate(force(mesh%np))
    call dcf_new_from(cf_for, ep%local_cf(1)) ! at least one specie must exist
    call dcf_alloc_FS(cf_for)
    call dcf_alloc_RS(cf_for)
    
    do i = 1, geo%natoms
      atm => geo%atom(i)
      do j = 1, conf%dim
        cf_for%FS = M_z0
        call cf_phase_factor(mesh, atm%x, ep%local_cf(atm%spec%index), cf_for)
        
        call dcf_FS_grad(mesh, cf_for, j)
        call dcf_FS2RS(cf_for)
        call dcf2mf(mesh, cf_for, force)
        do l = 1, st%d%nspin
          atm%f(j) = atm%f(j) + sum(force(1:mesh%np)*st%rho(1:mesh%np, l)*mesh%vol_pp(1:mesh%np))
        end do
      end do
    end do
    
    call dcf_free(cf_for)
    deallocate(force)
  end subroutine local_FS
#endif

end subroutine X(epot_forces)
