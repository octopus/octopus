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

subroutine R_FUNC(forces) (h, sys, t, no_lasers, lasers, reduce)
  type(hamiltonian_type), intent(IN) :: h
  type(system_type), intent(inout) :: sys
  integer, intent(in), optional :: no_lasers
  type(laser_type), intent(IN),optional :: lasers(:)
  real(r8), intent(in), optional :: t
  logical, intent(in), optional :: reduce

  integer :: i, j, l, m, add_lm, idim, ist, ik
  real(r8) :: d, r, x(3), zi, zj, vl, dvl
  R_TYPE :: uVpsi, p
  type(atom_type), pointer :: atm
#if defined(THREE_D)
  complex(r8), allocatable :: fw1(:,:,:), fw2(:,:,:)
  real(r8), allocatable :: fr(:,:,:), force(:)
#elif defined(ONE_D)
  complex(r8), allocatable :: fw1(:), fw2(:)
  real(r8), allocatable :: fr(:), force(:)
#endif

#if defined(HAVE_MPI) && defined(MPI_TD)
  real(r8) :: f(3)
  integer :: ierr
#endif 

  sub_name = 'forces'; call push_sub()

  ! And now the non-local part...
  ! this comes first to do the reduce...
  atm_loop: do i = 1, sys%natoms
    atm => sys%atom(i)
    atm%f(:) = 0._r8

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
              uVpsi = sum(atm%uV(:, add_lm) * sys%st%occ(ist, ik)  * &
                   sys%st%R_FUNC(psi)(atm%Jxyz(:), idim, ist, ik)) * &
                   sys%m%vol_pp**2 * atm%uVu(add_lm)
            
              do j = 1, 3
                p = sum(atm%duV(j, :, add_lm) * R_CONJ(sys%st%R_FUNC(psi) (atm%Jxyz(:), idim, ist, ik)))
                atm%f(j) = atm%f(j) + 2._r8 * R_REAL(uVpsi * p)
              end do

              add_lm = add_lm + 1
            end do m_loop
          end do l_loop

        end do dim_loop
      end do st_loop
    end do ik_loop

#if defined(HAVE_MPI) && defined(MPI_TD)
    if(present(reduce) .and. reduce) then
      call MPI_ALLREDUCE(atm%f(1), f(1), 3, &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      atm%f = f
    end if
#endif

  end do atm_loop
  
  if(present(no_lasers) .and. present(lasers)) then
    if(no_lasers>0) then
      call laser_field(no_lasers, lasers, t, x)
      do i = 1, sys%natoms
        sys%atom(i)%f(:) = sys%atom(i)%f(:) + sys%atom(i)%spec%Z_val * x(:)
      end do
    end if
  end if

  ! first the ion, ion force term
  do i = 1, sys%natoms
    zi = sys%atom(i)%spec%Z_val
    do j = 1, sys%natoms
      if(i .ne. j) then
        zj = sys%atom(j)%spec%Z_val
        r = sqrt(sum((sys%atom(i)%x(:) - sys%atom(j)%x(:))**2))
        d = zi * zj/r**3
 
        sys%atom(i)%f = sys%atom(i)%f + d*(sys%atom(i)%x(:) - sys%atom(j)%x(:))
      end if
    end do
  end do

  ! now comes the local part of the PP
#if defined(THREE_D)
  if(.not.atm%spec%local) then
#endif
    if(h%vpsl_space == 0) then ! Real space
      do i = 1, sys%natoms
        atm => sys%atom(i)
        
        do j = 1, sys%m%np
          call mesh_r(sys%m, j, r, x=x, a=sys%atom(i)%x)
          if(r < r_small) cycle
          
          vl  = splint(atm%spec%ps%vlocal, r)
          dvl = splint(atm%spec%ps%dvlocal, r)
          
          d = sum(sys%st%rho(j, :)) * sys%m%vol_pp* &
               (dvl - (vl - atm%spec%Z_val)/r)/r**2
          atm%f(:) = atm%f(:) + d * x(:)
        end do
      end do
    else ! Fourier space
#if defined(THREE_D)
      allocate( &
           fw1(sys%m%hfft_n2,   sys%m%fft_n2(2), sys%m%fft_n2(3)), &
           fw2(sys%m%hfft_n2,   sys%m%fft_n2(2), sys%m%fft_n2(3)), &
           fr (sys%m%fft_n2(1), sys%m%fft_n2(2), sys%m%fft_n2(3)), &
           force(sys%m%np))
      
      do i = 1, sys%natoms
        atm => sys%atom(i)
        do j = 1, 3
          fw1 = M_z0
          call phase_factor(sys%m, sys%m%fft_n2, atm%x, atm%spec%local_fw, fw1)
          call mesh_gradient_in_FS(sys%m, sys%m%hfft_n2, sys%m%fft_n2, fw1, fw2, j)
          
          call rfftwnd_f77_one_complex_to_real(sys%m%dplanb2, fw2, fr)
          force = 0._r8
          call dcube_to_mesh(sys%m, fr, force, t=2)
          do l = 1, sys%st%nspin
            atm%f(j) = atm%f(j) + sum(force(:)*sys%st%rho(:, l))*sys%m%vol_pp
          end do
        end do
      end do
      deallocate(fw1, fw2, fr, force)
#elif defined(ONE_D)
      allocate( &
           fw1(sys%m%hfft_n2), &
           fw2(sys%m%hfft_n2), &
           fr (sys%m%fft_n2(1)), &
           force(sys%m%np))
      
      do i = 1, sys%natoms
        atm => sys%atom(i)
        fw1 = M_z0
        call phase_factor(sys%m, sys%m%fft_n2, atm%x, atm%spec%local_fw, fw1)
        call mesh_gradient_in_FS(sys%m, sys%m%hfft_n2, sys%m%fft_n2, fw1, fw2)
        call rfftwnd_f77_one_complex_to_real(sys%m%dplanb2, fw2, fr)
        force = 0._r8
        call dcube_to_mesh(sys%m, fr, force, t=2)
        do l = 1, sys%st%nspin
          atm%f(1) = atm%f(1) + sum(force(:)*sys%st%rho(:, l))*sys%m%vol_pp
        end do
      end do
      deallocate(fw1, fw2, fr, force)
#endif
    end if
#if defined(THREE_D)
  end if
#endif

  call pop_sub()
end subroutine R_FUNC(forces)
