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

! calculates the eigenvalues of the real orbitals
subroutine R_FUNC(hamiltonian_eigenval)(h, st, sys)
  type(hamiltonian_type), intent(IN) :: h
  type(states_type), intent(inout) :: st
  type(system_type), intent(in) :: sys

  R_TYPE, allocatable :: Hpsi(:,:)
  R_TYPE :: e
  integer :: ik, ist

  call push_sub('hamiltonian_eigenval')
  allocate(Hpsi(sys%m%np, st%dim))

  do ik = 1, st%nik
    do ist = st%st_start, st%st_end
      call R_FUNC(Hpsi) (h, sys%m, st, sys, ik, st%R_FUNC(psi)(:, :, ist, ik), Hpsi)
      e = R_FUNC(states_dotp)(sys%m, st%dim, st%R_FUNC(psi)(1:, :, ist, ik), Hpsi)
      st%eigenval(ist, ik) = R_REAL(e)
    end do
  end do

  deallocate(Hpsi)
  call pop_sub()
end subroutine R_FUNC(hamiltonian_eigenval)

subroutine R_FUNC(Hpsi) (h, m, st, sys, ik, psi, Hpsi, t)
  type(hamiltonian_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  type(system_type), intent(in) :: sys ! this is necessary due to the nl part of the PP
  integer, intent(in) :: ik
  R_TYPE, intent(IN) :: psi(m%np, st%dim)
  R_TYPE, intent(out) :: Hpsi(m%np, st%dim)
  real(r8), intent(in), optional :: t

  call push_sub('Hpsi')

  call R_FUNC(kinetic) (h, ik, m, st, psi, Hpsi)
  call R_FUNC(vlpsi)   (h, m, st, ik, psi, Hpsi)
  if(sys%nlpp) call R_FUNC(vnlpsi)  (ik, m, st, sys, psi, Hpsi)

  ! Relativistic corrections...
  select case(h%reltype)
  case(NOREL)
#if defined(COMPLEX_WFNS) && defined(R_TCOMPLEX)
  case(SPIN_ORBIT)
    call zso (h, m, sys%natoms, sys%atom, st%dim, ik, psi, Hpsi)
#endif
  case default
    message(1) = 'Error: Internal.'
    call write_fatal(1)
  end select

  if(present(t)) then
    call R_FUNC(vlasers)  (h, m, st, psi, Hpsi, t)
    call R_FUNC(vborders) (h, m, st, psi, Hpsi)
  endif

  call pop_sub()
end subroutine R_FUNC(Hpsi)

subroutine R_FUNC(kinetic) (h, ik, m, st, psi, Hpsi)
  type(hamiltonian_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  R_TYPE, intent(IN) :: psi(m%np, st%dim)
  R_TYPE, intent(out) :: Hpsi(m%np, st%dim)
  integer :: ik

  R_TYPE, allocatable :: grad(:,:)
  real(r8) :: k2
  integer :: idim, i, np, dim

  call push_sub('kinetic')
  np = m%np
  dim = st%dim

  if(conf%periodic_dim>0) then
#if defined(COMPLEX_WFNS)
    allocate(grad(3, m%np))
    k2 = sum(st%kpoints(:, ik)**2)
    do idim = 1, dim
      call R_FUNC(mf_laplacian) (m, psi(:, idim), Hpsi(:, idim), cutoff_ = h%cutoff)
      call R_FUNC(mf_gradient)  (m, psi(:, idim), grad(:, :))
      do i = 1, m%np
        Hpsi(i, idim) = -M_HALF*(Hpsi(i, idim) &
             + M_TWO*M_zI*sum(st%kpoints(:, ik)*grad(:, i)) &
             - k2*psi(i, idim))
      end do
    end do
    deallocate(grad)
#else
    message(1) = "Real wavefunction for ground state not yet implemented for polymers:"
    message(2) = "Reconfigure with --enable-complex, and remake"
    call write_fatal(2)
#endif
  
  else
    do idim = 1, dim
      call R_FUNC(mf_laplacian) (m, psi(:, idim), Hpsi(:, idim), cutoff_ = h%cutoff)
    end do
    call R_FUNC(scal)(m%np*dim, R_TOTYPE(-M_HALF), Hpsi, 1)
  end if


  call pop_sub()
end subroutine R_FUNC(kinetic)

subroutine R_FUNC(exp_kinetic) (h, ik, m, st, psi, factor)
  type(hamiltonian_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  R_TYPE :: psi(:, :)
  integer :: ik
  R_TYPE, intent(in) :: factor

  complex(r8), allocatable :: fw(:,:,:)
  integer :: i, n(3), idim, np, dim
  real(r8) :: k2

  call push_sub('exp_kinetic')
  np = m%np
  dim = st%dim

  if(conf%periodic_dim>0) then
     message(1) = 'Internal bug in exp_kinetic'
     call write_fatal(1)
  endif

#ifdef T_REAL
  call fft_getdim_complex(m%dfft, n)
#else
  call fft_getdim_real   (m%dfft, n)
#endif
  allocate(fw(n(1), n(2), n(3)))
  do idim = 1, dim
     call R_FUNC(mesh_rs2fs)(m, psi(:, idim), fw)
     call R_FUNC(mesh_laplq)(m, fw, n, exponential = -factor/M_TWO, cutoff = h%cutoff)
     call R_FUNC(mesh_fs2rs)(m, fw, psi(:, idim))
  enddo
  deallocate(fw)

  call pop_sub()
end subroutine R_FUNC(exp_kinetic)

subroutine R_FUNC(vnlpsi) (ik, m, st, sys, psi, Hpsi)
  integer, intent(in) :: ik
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  type(system_type), intent(IN) :: sys
  R_TYPE, intent(IN) :: psi(m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)

  integer :: i, is, idim, ia, ikbc, jkbc, k, l, lm, add_lm
  real(r8) :: x(conf%dim)
  R_TYPE :: phase
  R_TYPE :: uVpsi
  R_TYPE, allocatable :: lpsi(:), lHpsi(:)
  type(atom_type), pointer :: atm
  type(specie_type), pointer :: spec
  R_TYPE, external :: R_DOT
  

  call push_sub('vnlpsi')

  ! Ionic pseudopotential
  do_atm: do ia = 1, sys%natoms
    atm => sys%atom(ia)
    spec => atm%spec
    ! do we have a pseudopotential, or a local pot?
    if(spec%local) cycle do_atm

    ! calculate the phase factors exp(ik*x) 
    ! they are necessary to get the correct periodicity 
    ! for the nonlocal part of the pseudopotential
    if (conf%periodic_dim/=0) then
      allocate(atm%phases(atm%mps,st%nik))
      do i=1,atm%mps
        call mesh_xyz(m, atm%jxyz(i), x)
        do k=1,st%nik
          atm%phases(i,k)=exp(M_zI*sum(st%kpoints(:,k)*x(:)))
        end do
      end do
    end if

    do_dim: do idim = 1, st%dim
      allocate(lpsi(atm%mps), lHpsi(atm%mps))
      if (conf%periodic_dim==0) then
        lpsi(:) = psi(atm%jxyz(:), idim)
      else
        lpsi(:) = atm%phases(:,ik)*psi(atm%jxyz(:), idim)
      end if
      lHpsi(:) = M_z0
      add_lm = 1
      do_l: do l = 0, spec%ps%l_max
        if (l == spec%ps%L_loc) then
          add_lm = add_lm + (2*l + 1)
          cycle do_l
        end if
        
        do_m: do lm = -l, l
          do ikbc = 1, spec%ps%kbc
            do jkbc = 1, spec%ps%kbc
              uvpsi = R_DOT(atm%mps, atm%R_FUNC(uv)(:, add_lm, ikbc), 1, lpsi(:), 1) * &
                   m%vol_pp*atm%R_FUNC(uvu)(add_lm, ikbc, jkbc) 
              if (conf%periodic_dim==0) then
                call R_FUNC(axpy) (atm%mps, uvpsi, atm%R_FUNC(uv)(:, add_lm, jkbc), 1, lHpsi(:), 1)
              else
                call R_FUNC(axpy) (atm%mps, uvpsi, R_CONJ(atm%phases(:,ik))*atm%R_FUNC(uv)(:, add_lm, jkbc), 1, lHpsi(:), 1)
              end if
            end do
          end do
          
          add_lm = add_lm + 1
        end do do_m
      end do do_l
      Hpsi(atm%jxyz(:), idim) = Hpsi(atm%jxyz(:), idim) + lHpsi(:)
      deallocate(lpsi, lHpsi)
    end do do_dim
  end do do_atm

  call pop_sub()
end subroutine R_FUNC(vnlpsi)

! For the moment, let us define this function only in the complex case...
#if defined(R_TCOMPLEX)
subroutine R_FUNC(exp_vnlpsi) (ik, m, st, sys, psi, factor, order)
  integer, intent(in) :: ik
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  type(system_type), intent(in) :: sys
  R_TYPE :: psi(m%np, st%dim)
  R_TYPE, intent(in) :: factor
  logical, intent(in) :: order

  integer :: is, idim, ia, ikbc, jkbc, l, lm, add_lm, &
             ia_start, ia_end, step, l_start, l_end, kbc_start, kbc_end
  complex(r8) :: uvpsi, p2, ctemp
  complex(r8), allocatable :: lpsi(:), lhpsi(:), initzpsi(:, :)
  type(atom_type), pointer :: atm
  type(specie_type), pointer :: spec
  complex(r8), external :: zdotc

  call push_sub('vnlpsi')

    if(order) then
      step = 1;  ia_start = 1; ia_end = sys%natoms
    else
      step = -1; ia_start = sys%natoms; ia_end = 1
    end if

    allocate(initzpsi(m%np, 1:sys%st%dim))
    initzpsi = psi
    ! Ionic pseudopotential
    do_atm: do ia = ia_start, ia_end, step
      atm => sys%atom(ia)
      spec => atm%spec
      ! do we have a pseudopotential, or a local pot?
      if(spec%local) cycle do_atm

      do_dim: do idim = 1, sys%st%dim
        allocate(lpsi(atm%mps), lHpsi(atm%mps))
        lpsi(:) = initzpsi(atm%jxyz(:), idim)
        lHpsi(:) = M_z0
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
            add_lm = add_lm + (2*l + 1)*step
            cycle do_l
          end if

          do_m: do lm = -l*step, l*step, step
            do ikbc = kbc_start, kbc_end, step
              do jkbc = kbc_start, kbc_end, step
                 p2 = zdotc(atm%mps, atm%zuv(:, add_lm, ikbc), 1, atm%zuv(:, add_lm, ikbc), 1)*m%vol_pp
                 ctemp = atm%zuvu(add_lm, ikbc, jkbc)*p2*factor
                 uvpsi = zdotc(atm%mps, atm%zuv(:, add_lm, ikbc), 1, lpsi(:), 1) * m%vol_pp* &
                       (exp(ctemp) - M_z1)/p2
                 call zaxpy (atm%mps, uvpsi, atm%zuv(:, add_lm, jkbc), 1, lHpsi(:), 1)
              end do
            end do
            add_lm = add_lm + step
          end do do_m
        end do do_l
        psi(atm%jxyz(:), idim) = psi(atm%jxyz(:), idim) + lhpsi(:)
        deallocate(lpsi, lHpsi)
      end do do_dim
    end do do_atm

    deallocate(initzpsi)
    call pop_sub()
end subroutine R_FUNC(exp_vnlpsi)
#endif

subroutine R_FUNC(vlpsi) (h, m, st, ik, psi, Hpsi)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  integer, intent(in) :: ik
  R_TYPE, intent(in) :: psi(m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)

  integer :: is, idim, np, dim

  call push_sub('vlpsi')

  np = m%np
  dim = st%dim

    do idim = 1, dim
      Hpsi(:, idim) = Hpsi(:, idim) + (h%Vpsl + h%Vhartree)*psi(1:,idim)
    end do

    select case(st%ispin)
    case(UNPOLARIZED)
      hpsi(:, 1) = hpsi(:, 1) + h%vxc(:, 1)*psi(1:, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        hpsi(:, 1) = hpsi(:, 1) + h%vxc(:, 1)*psi(1:, 1)
      else
        hpsi(:, 1) = hpsi(:, 1) + h%vxc(:, 2)*psi(1:, 1)
      end if
    case(SPINORS)
      hpsi(:, 1) = hpsi(:, 1) + h%vxc(:, 1)*psi(1:, 1) + &
                   (h%vxc(:, 3) - M_zI*h%vxc(:, 4))*psi(1:, 2)
      hpsi(:, 2) = hpsi(:, 2) + h%vxc(:, 2)*psi(1:, 2) + &
                   (h%vxc(:, 3) + M_zI*h%vxc(:, 4))*psi(1:, 1)
    end select

  call pop_sub()
end subroutine R_FUNC(vlpsi)

subroutine R_FUNC(exp_vlpsi) (h, m, st, ik, psi, t, factor)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  integer, intent(in) :: ik
  R_TYPE, intent(inout) :: psi(m%np, st%dim)
  !R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)
  real(r8), intent(in) :: t
  R_TYPE, intent(in) :: factor

  integer :: is, idim, np, dim, k
  real(r8) :: x(3), f(3)

  call push_sub('vlpsi')

  !!!!WARNING: spinors not yet supported.
  np = m%np
  dim = st%dim

  do idim = 1, dim
     psi(:, idim) = exp(factor*(h%vpsl(:)+h%vhartree(:)))*psi(:, idim)
  end do

  select case(st%ispin)
    case(UNPOLARIZED)
      psi(:, 1) = exp(factor*h%vxc(:, 1))*psi(:, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        psi(:, 1) = exp(factor*h%vxc(:, 1))*psi(:, 1)
      else
        psi(:, 1) = exp(factor*h%vxc(:, 2))*psi(:, 1)
      end if
    case(SPINORS)
      message(1) = 'Internal error in exp_vlpsi'
      call write_fatal(1)
  end select

  if(h%no_lasers > 0) then
      select case(h%gauge)
      case(1) ! length gauge
        call laser_field(h%no_lasers, h%lasers, t, f)
        do k = 1, m%np
          call mesh_xyz(m, k, x)
          psi(k,:) = exp(factor*sum(x*f)) * psi(k,:)
        end do
      case(2) ! velocity gauge
        message(1) = 'Internal error in exp_vlpsi'
        call write_fatal(1)
      end select
  end if

  call pop_sub()
end subroutine R_FUNC(exp_vlpsi)

subroutine R_FUNC(vlasers) (h, m, st, psi, Hpsi, t)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  R_TYPE, intent(in) :: psi(m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)
  real(r8), intent(in) :: t

  integer :: is, i, k, idim, np, dim
  real(r8) :: x(3), f(3)
  R_TYPE, allocatable :: grad(:,:)

  call push_sub('vlasers')

  if(h%no_lasers > 0) then
    select case(h%gauge)
    case(1) ! length gauge
      call laser_field(h%no_lasers, h%lasers, t, f)
      
      do k = 1, m%np
        call mesh_xyz(m, k, x)
        hpsi(k,:) = hpsi(k,:) + sum(x*f) * psi(k,:)
      end do
      
    case(2) ! velocity gauge
      call laser_vector_field(h%no_lasers, h%lasers, t, f)
      allocate(grad(3, m%np))
      do idim = 1, st%dim
        call R_FUNC(mf_gradient)(m, psi(:, idim), grad)
        do k = 1, m%np
          hpsi(k, idim) = hpsi(k, idim) - M_zI * sum(f(:)*grad(:, k)) + &
               sum(f**2)/2._r8 * psi(k, idim)
        end do
      end do
      deallocate(grad)
    end select
  end if
  
  call pop_sub()
end subroutine R_FUNC(vlasers)

subroutine R_FUNC(vborders) (h, m, st, psi, Hpsi)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  R_TYPE, intent(in) :: psi(m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)

  integer :: idim

  call push_sub('vborders')

  if(h%ab .eq. IMAGINARY_ABSORBING) then
    do idim = 1, st%dim
       hpsi(:, idim) = hpsi(:, idim) + M_zI*h%ab_pot(:)*psi(1:, idim)
    end do
  end if

  call pop_sub()
end subroutine R_FUNC(vborders)

subroutine R_FUNC(hamiltonian_setup)(h, m, st, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(inout) :: st
  type(system_type), intent(in), optional :: sys

  integer :: is

  call push_sub('hamiltonian_setup')

  if(.not. h%ip_app) then ! No Hartree or xc if independent electrons.
    call hartree_solve(h%hart, m, h%vhartree, st%rho(:, 1:st%spin_channels))
    h%epot = M_ZERO
    do is = 1, st%spin_channels
       h%epot = h%epot - M_HALF*dmf_dotp(m, st%rho(:, is), h%vhartree)
    enddo

    call R_FUNC(xc_pot)(h%xc, m, st, h%hart, h%vxc, h%ex, h%ec, &
                        -minval(st%eigenval(st%nst, :)), st%qtot)

    select case(h%ispin)
      case(UNPOLARIZED)
        h%epot = h%epot - dmf_dotp(m, st%rho(:, 1), h%vxc(:, 1))
      case(SPIN_POLARIZED)
        h%epot = h%epot - dmf_dotp(m, st%rho(:, 1), h%vxc(:, 1)) &
                        - dmf_dotp(m, st%rho(:, 2), h%vxc(:, 2))
      case(SPINORS)
        h%epot = h%epot - dmf_dotp(m, st%rho(:, 1), h%vxc(:, 1)) &
                        - dmf_dotp(m, st%rho(:, 2), h%vxc(:, 2))
        h%epot = h%epot - M_TWO*zmf_dotp(m, st%rho(:, 3) + M_zI*st%rho(:, 4), &
                                             h%vxc(:, 3) - M_zI* h%vxc(:, 4))
    end select
  endif

  ! this, I think, belongs here
  if(present(sys)) call R_FUNC(hamiltonian_eigenval) (h, st, sys)

  call pop_sub()
end subroutine R_FUNC(hamiltonian_setup)
