!! Suzuki-Trotter split-operator method

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

!WARNING - velocity gauge is not properly implemented!

subroutine td_rti3(sys, h, td, tt)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  type(td_type), intent(inout) :: td
  real(r8), intent(in) :: tt

  integer :: ik, ist, idim
  real(r8) :: t, temp(3), las(3)

  temp(:) = 2.0_r8*M_PI/(sys%m%fft_n(:)*sys%m%h(:))

  ! get reference time right
  t = tt - td%dt

  ! get lasers for velocity gauge
  if(td%no_lasers > 0 .and. td%gauge == 2) then 
    call laser_vector_field(td%no_lasers, td%lasers, t, las)
  end if

  ! propagate with T/2
  do ik = 1, sys%st%nik
    do ist = sys%st%st_start, sys%st%st_end
      do idim = 1, sys%st%dim
        ! kinetic contribution (in momentum space)
        call kinetic(sys%m, td%dt/2._r8)

        ! non-local part of the pseudopotential
        call non_local_pp(sys%m, td%dt/2._r8, .true.)
      end do
    end do
  end do

  ! setup hamiltonian
  call zcalcdens(sys%st, sys%m%np, sys%st%rho, .true.)
  call zhamiltonian_setup(h, sys)

  ! propagate local part of the potential
  call local_part(sys%m, td%dt)

  ! propagate with T/2
  do ik = 1, sys%st%nik
    do ist = sys%st%st_start, sys%st%st_end
      do idim = 1, sys%st%dim
        ! non-local part of the pseudopotential
        call non_local_pp(sys%m, td%dt/2._r8, .false.)

        ! kinetic contribution (in momentum space)
        call kinetic(sys%m, td%dt/2._r8)
      end do
    end do
  end do

contains
  subroutine kinetic(m, dt)
    type(mesh_type), intent(IN) :: m
    real(r8), intent(in) :: dt

    integer :: ix, iy, iz, ixx(3)
    real(r8) :: vec
    complex(r8), allocatable :: wf_r(:,:,:), wf_k(:,:,:)

    allocate(&
         wf_r(sys%m%fft_n(1), sys%m%fft_n(2), sys%m%fft_n(3)), &
         wf_k(sys%m%fft_n(1), sys%m%fft_n(2), sys%m%fft_n(3))  &
         )

    wf_r = M_z0
    call zmesh_to_cube(sys%m, sys%st%zpsi(1:,idim, ist, ik), wf_r)
    call fftwnd_f77_one(sys%m%zplanf, wf_r, wf_k)
    
    do ix = 1, m%fft_n(1)
      ixx(1) = pad_feq(ix, m%fft_n(1), .true.)
      do iy = 1, m%fft_n(2)
        ixx(2) = pad_feq(iy, m%fft_n(2), .true.)
        do iz = 1, m%fft_n(3)
          ixx(3) = pad_feq(iz, m%fft_n(3), .true.)
          vec = sum((temp(:)*ixx(:))**2)/2._r8
          
          ! TODO add vector potential
          wf_k(ix, iy, iz) = exp(- M_zI*dt*vec) * wf_k(ix, iy, iz)
        end do
      end do
    end do
    
    call fftwnd_f77_one(sys%m%zplanb, wf_k, wf_r)
    wf_r = wf_r/real(sys%m%fft_n(1)*sys%m%fft_n(2)*sys%m%fft_n(3), r8)

    sys%st%zpsi(1:,idim, ist, ik)= M_z0
    call zcube_to_mesh(sys%m, wf_r, sys%st%zpsi(1:,idim, ist, ik))
    
    deallocate(wf_r, wf_k)

!!$      if(td%no_lasers > 0 .and. td%gauge == 2) then
!!$        do ix = 1, td%fft_n
!!$          do iy = 1, td%fft_n
!!$            do iz = 1, td%fft_n
!!$              ixx = pad_feq(ix, td%fft_n, .true.)
!!$              iyy = pad_feq(iy, td%fft_n, .true.)
!!$              izz = pad_feq(iz, td%fft_n, .true.)
!!$              
!!$              wf(ix, iy, iz) = wf(ix, iy, iz) * exp(- M_zI*td%dt/2._r8*( &
!!$                   2._r8*temp*(las(1)*real(ixx,r8) + las(2)*real(iyy,r8) + las(3)*real(izz,r8)) &
!!$                   + sum(f**2)))
!!$            end do
!!$          end do
!!$        end do
!!$      end if
  end subroutine kinetic

  subroutine non_local_pp(m, dt, order)
    type(mesh_type), intent(IN) :: m
    real(r8), intent(in) :: dt
    logical, intent(in) :: order

    integer :: step, ia, ia_start, ia_end, l, l_start, l_end, lm, add_lm
    integer :: ikbc, jkbc, kbc_start, kbc_end
    complex(r8) :: uVpsi, ctemp
    type(atom_type), pointer :: atm
    type(specie_type), pointer :: spec

    if(order) then
      step = 1;  ia_start = 1; ia_end = sys%natoms
    else
      step = -1; ia_start = sys%natoms; ia_end = 1
    end if

    do_atm: do ia = ia_start, ia_end, step
      atm => sys%atom(ia)
      spec => atm%spec
      
      if(spec%local) cycle do_atm
      
      if(order) then
        l_start   = 0; l_end   = spec%ps%L_max
        kbc_start = 1; kbc_end = spec%ps%kbc
        add_lm = 1
      else
        l_start   = spec%ps%L_max; l_end = 0
        kbc_start = spec%ps%kbc; kbc_end = 1
        add_lm = (spec%ps%L_max + 1)**2
      end if

      do_l: do l = l_start, l_end, step
        if (l == spec%ps%L_loc) then
          add_lm = add_lm + step*(2*l + 1)
          cycle do_l
        end if

        do_lm: do lm = -l*step, l*step, step
          do ikbc = kbc_start, kbc_end, step
            do jkbc = kbc_start, kbc_end, step

              uVpsi = sum(atm%uV(:, add_lm, ikbc)*sys%st%zpsi(atm%Jxyz(:), idim, ist, ik))*sys%m%vol_pp
              ctemp = uVpsi * (exp(-M_zI*dt*atm%uVu(add_lm, ikbc, jkbc)) - 1.0_r8)
          
              sys%st%zpsi(atm%Jxyz(:), idim, ist, ik) = sys%st%zpsi(atm%Jxyz(:), idim, ist, ik) + &
                   ctemp * atm%uV(:, add_lm, jkbc)
            end do
          end do
          
          add_lm = add_lm + step
        end do do_lm
      end do do_l

    end do do_atm

  end subroutine non_local_pp

  subroutine local_part(m, dt)
    type(mesh_type), intent(IN) :: m
    real(r8), intent(in) :: dt

    integer :: j, ist, ik, ix, iy, iz
    real(r8) :: r, x(3)
    
    if(td%no_lasers > 0 .and. td%gauge == 1) then 
      call laser_field(td%no_lasers, td%lasers, t + dt/2._r8, las)
    end if

    ! Propagation of the local part and external field
    do ik = 1, sys%st%nik
      do ist = sys%st%st_start, sys%st%st_end

        do j = 1, m%np
          ix = m%Lx(j) + (m%fft_n(1) - 1)/2 + 1
          iy = m%Ly(j) + (m%fft_n(2) - 1)/2 + 1
          iz = m%Lz(j) + (m%fft_n(3) - 1)/2 + 1

          r = h%Vpsl(j) + h%Vhartree(j)
          select case(sys%st%ispin)
          case(1) ! dim = 1
            r = r + h%Vxc(j, 1)
          case(2) ! dim = 1
            if(modulo(ik, 2) == 0) then ! we have a spin down
              r = r + h%Vxc(j, 1)
            else ! spin down
              r = r + h%Vxc(j, 2)
            end if
          case(3) ! dim = 2
            message(1) = "Suzuki-Trotter propagation not implemented for ispin=3"
            call write_fatal(1)
          end select

          if(td%ab .eq. 1) then
            r = r + td%ab_pot(j)
          end if

          if(td%no_lasers > 0 .and. td%gauge == 1) then
            call mesh_xyz(m, j, x)
            r = r + sum(x(:)*las(:))
          end if

          ! WARNING: dimension is hard-wired to 1
          sys%st%zpsi(j, 1, ist, ik) = sys%st%zpsi(j, 1, ist, ik) * &
               exp(- M_zI*td%dt*r)

        end do
      end do
    end do
  end subroutine local_part

end subroutine td_rti3
