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

subroutine R_FUNC(forces) (h, sys, t, reduce)
  type(hamiltonian_type), intent(IN) :: h
  type(system_type), intent(inout) :: sys
  real(r8), intent(in), optional :: t
  logical, intent(in), optional :: reduce

  integer :: i, j, l, m, add_lm, idim, ist, ik, ii, jj
  real(r8) :: d, r, x(3), zi, zj, vl, dvl
  R_TYPE :: uVpsi, p
  type(atom_type), pointer :: atm

#if defined(HAVE_MPI) && defined(MPI_TD)
  real(r8) :: f(3)
  integer :: ierr
#endif 

  call push_sub('forces')

  ! init to 0
  do i = 1, sys%natoms
    sys%atom(i)%f = M_ZERO
  end do

  ! And now the non-local part...
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
                 uVpsi = sum(atm%R_FUNC(uV)(:, add_lm, ii) * sys%st%occ(ist, ik)  * &
                         sys%st%R_FUNC(psi)(atm%Jxyz(:), idim, ist, ik)) * &
                         sys%m%vol_pp**2 * atm%R_FUNC(uVu)(add_lm, ii, jj)

                 do j = 1, 3
                    p = sum(atm%R_FUNC(duV)(j, :, add_lm, jj) * R_CONJ(sys%st%R_FUNC(psi) (atm%Jxyz(:), idim, ist, ik)))
                   atm%f(j) = atm%f(j) + 2._r8 * R_REAL(uVpsi * p)
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
    if(present(reduce)) then
    if(reduce) then
      call MPI_ALLREDUCE(atm%f(1), f(1), conf%dim, &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      atm%f = f
    end if
    end if
#endif

  end do atm_loop
  
  if(present(t).and.h%no_lasers>0) then
    call laser_field(h%no_lasers, h%lasers, t, x)
    do i = 1, sys%natoms
      sys%atom(i)%f(1:conf%dim) = sys%atom(i)%f(1:conf%dim) + &
           sys%atom(i)%spec%Z_val * x(1:conf%dim)
    end do
  end if

  ! first the ion, ion force term
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
  if(h%vpsl_space == 0) then ! Real space
    call local_RS()
  else ! Fourier space
    call local_FS()
  end if

  call pop_sub()

  contains
    subroutine local_RS()
      real(r8) :: r, x(3), d, gv(3)
      integer  :: ns

      ns = min(2, sys%st%nspin)

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

    ! WARNING: this is still not working
    subroutine local_FS()
      complex(r8), allocatable :: fw(:,:,:)
      real(r8), allocatable :: fr(:,:,:), force(:)
      integer :: db(3), dbc(3)

      call fft_getdim_real   (h%fft, db)
      call fft_getdim_complex(h%fft, dbc)
      
      allocate(fw(dbc(1), dbc(2), dbc(3)), fr(db(1), db(2), db(3)), force(sys%m%np))

      do i = 1, sys%natoms
        atm => sys%atom(i)
        do j = 1, conf%dim
          fw = M_z0
          call phase_factor(sys%m, db, atm%x, atm%spec%local_fw, fw)
          call dmesh_gradq(sys%m, fw, j, db)
          call rfft_backward(h%fft, fw, fr)
          force = M_ZERO
          call dcube_to_mesh(sys%m, fr, force, db)
          do l = 1, sys%st%nspin
            atm%f(j) = atm%f(j) + sum(force(:)*sys%st%rho(:, l))*sys%m%vol_pp
          end do
        end do
      end do

      deallocate(fw, fr, force)
    end subroutine local_FS

end subroutine R_FUNC(forces)
