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

subroutine X(epot_forces) (ep, sys, t, reduce_)
  type(epot_type),   intent(in)    :: ep
  type(system_type), intent(inout) :: sys
  FLOAT,          intent(in), optional :: t
  logical,           intent(in), optional :: reduce_

  integer :: i, j, l, m, add_lm, idim, ist, ik, ii, jj
  FLOAT :: d, r, zi, zj, vl, dvl, x(3)
  R_TYPE :: uVpsi, p
  type(atom_type), pointer :: atm

#if defined(HAVE_MPI) && defined(MPI_TD)
  FLOAT :: f(3)
  integer :: ierr
#endif 
  
  ! init to 0
  do i = 1, sys%natoms
    sys%atom(i)%f = M_ZERO
  end do

  ! the non-local contribution
  ! this comes first to do the reduce...
  atm_loop: do i = 1, sys%natoms
    atm => sys%atom(i)
    if(atm%spec%local) cycle
    
    ik_loop: do ik = 1, sys%st%nik
      st_loop: do ist = sys%st%st_start, sys%st%st_end
        dim_loop: do idim = 1, sys%st%dim
          add_lm = 1
          l_loop: do l = 0, atm%spec%ps%L_max
            if(l == atm%spec%ps%L_loc) then
              add_lm = add_lm + (2*l + 1)
              cycle l_loop
            end if
            
            m_loop: do m = -l, l
              do ii = 1, atm%spec%ps%kbc
                do jj = 1, atm%spec%ps%kbc
                  uVpsi = sum(atm%X(uV)(:, add_lm, ii) * sys%st%occ(ist, ik)  * &
                       sys%st%X(psi)(atm%Jxyz(:), idim, ist, ik)) * &
                       sys%m%vol_pp**2 * atm%X(uVu)(add_lm, ii, jj)
                  
                  do j = 1, 3
                    p = sum(atm%X(duV)(j, :, add_lm, jj) * R_CONJ(sys%st%X(psi) (atm%Jxyz(:), idim, ist, ik)))
                    atm%f(j) = atm%f(j) + M_TWO * R_REAL(uVpsi * p)
                  end do
                end do
              end do
              
              add_lm = add_lm + 1
            end do m_loop
          end do l_loop
          
        end do dim_loop
      end do st_loop
    end do ik_loop
    
#if defined(HAVE_MPI) && defined(MPI_TD)
    if(present(reduce_)) then
      if(reduce_) then
        call MPI_ALLREDUCE(atm%f(1), f(1), conf%dim, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        atm%f = f
      end if
    end if
#endif
    
  end do atm_loop
  
  ! Now the ion, ion force term
  do i = 1, sys%natoms
    zi = sys%atom(i)%spec%Z_val
    do j = 1, sys%natoms
      if(i .ne. j) then
        zj = sys%atom(j)%spec%Z_val
        r = sqrt(sum((sys%atom(i)%x(1:conf%dim) - sys%atom(j)%x(1:conf%dim))**2))
        d = zi * zj/r**3
        
        sys%atom(i)%f(1:conf%dim) = sys%atom(i)%f(1:conf%dim) + &
             d*(sys%atom(i)%x(1:conf%dim) - sys%atom(j)%x(1:conf%dim))
      end if
    end do
  end do
  
  ! now comes the local part of the PP
  if(ep%vpsl_space == 0) then ! Real space
    call local_RS()
#ifdef HAVE_FFT
  else ! Fourier space
    call local_FS()
#endif
  end if

  if(present(t).and.ep%no_lasers>0) then
    call epot_laser_field(ep, t, x)
    do i = 1, sys%natoms
      sys%atom(i)%f(1:conf%dim) = sys%atom(i)%f(1:conf%dim) + &
           sys%atom(i)%spec%Z_val * x(1:conf%dim)
    end do
  end if
  
contains
  subroutine local_RS()
    FLOAT :: r, x(3), d, gv(3)
    integer  :: ns
    
    ns = min(2, sys%st%d%nspin)
    
    do i = 1, sys%natoms
      atm => sys%atom(i)
      do j = 1, sys%m%np
        call mesh_r(sys%m, j, r, x=x, a=sys%atom(i)%x)
        if(r < r_small) cycle
        
        call specie_get_glocal(atm%spec, x, gv)
        d = sum(sys%st%rho(j, 1:ns))*sys%m%vol_pp
        atm%f(:) = atm%f(:) - d*gv(:)
      end do
    end do
  end subroutine local_RS
  
#ifdef HAVE_FFT
  subroutine local_FS()
    type(dcf) :: cf_for
    FLOAT, allocatable :: force(:)
    
    allocate(force(sys%m%np))
    call dcf_new_from(cf_for, ep%local_cf(1)) ! at least one specie must exist
    call dcf_alloc_FS(cf_for)
    call dcf_alloc_RS(cf_for)
    
    do i = 1, sys%natoms
      atm => sys%atom(i)
      do j = 1, conf%dim
        cf_for%FS = M_z0
        call cf_phase_factor(sys%m, atm%x, ep%local_cf(atm%spec%index), cf_for)
        
        call dcf_FS_grad(sys%m, cf_for, j)
        call dcf_FS2RS(cf_for)
        call dcf2mf(sys%m, cf_for, force)
        do l = 1, sys%st%d%nspin
          atm%f(j) = atm%f(j) + sum(force(:)*sys%st%rho(:, l))*sys%m%vol_pp
        end do
      end do
    end do
    
    deallocate(force)
  end subroutine local_FS
#endif

end subroutine X(epot_forces)
