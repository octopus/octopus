#include "config_F90.h"

module lcao
  use global
  use liboct
  use spline
  use mesh
  use system
  use hamiltonian
  use states
  use mix

  implicit none

contains

!builds a density which is the sum of the atomic densities
subroutine lcao_dens(sys, nspin, rho)
  type(system_type), intent(IN) :: sys
  integer, intent(in) :: nspin
  real(r8), intent(out) :: rho(sys%m%np, nspin)
  
  integer :: ia
  real(r8) :: r
  type(specie_type), pointer :: s
  type(atom_type),   pointer :: a
  
  rho = 0._r8
  do ia = 1, sys%natoms
    a => sys%atom(ia) ! shortcuts
    s => a%spec
    
    select case(s%label(1:5))
    case('jelli', 'point')
      call from_jellium(sys%m, rho(:, 1))
    case default
      call from_pseudopotential(sys%m, rho(:, 1))
    end select
  end do

  ! we now renormalize the density (necessary if we have a charged system)
  ! if spin polarized, we start with paramagnetic density
  r = dmesh_integrate(sys%m, rho(:, 1))
  write(message(1),'(a,f13.6)')'Info: Unnormalized total charge = ', r
  call write_info(1)

  r = sys%st%qtot/(r*real(nspin, r8))
  do ia = nspin, 1, -1
    rho(:, ia) = r*rho(:, 1)
  end do

contains
  subroutine from_jellium(m, rho)
    type(mesh_type), intent(in) :: m
    real(r8), intent(inout) :: rho(m%np)

    integer :: i, in_points
    real(r8) :: r

    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
      if(r <= s%jradius) then
        in_points = in_points + 1
      end if
    end do
    
    if(in_points > 0) then
      do i = 1, m%np
        call mesh_r(m, i, r, a=a%x)
        if(r <= s%jradius) then
          rho(i) = rho(i) + real(s%Z_val, r8)/(in_points*m%vol_pp)
        end if
      end do
    end if
  end subroutine from_jellium
  
  subroutine from_pseudopotential(m, rho)
    type(mesh_type), intent(in) :: m
    real(r8), intent(inout) :: rho(m%np)

    integer :: i, l
    real(r8) :: r, zel, zval
    R_TYPE :: psi

    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
      zval = s%Z_val
#if defined(THREE_D)
      do l = 0 , s%ps%L_max
        if(r >= r_small) then
          psi = splint(s%ps%Ur(l), r)
          zel = min(zval, 2.0_r8*(2*l+1))
          zval = max(0._r8, zval - zel)
          rho(i) = rho(i) + zel*psi*psi*(r**(2*l))/(4*M_PI)
        end if
      end do
#elif defined(ONE_D)
      !rho(i) = rho(i) + s%z_val*exp(-r**2)/sqrt(m_pi)
      call R_FUNC(calcdens)(sys%st, m%np, rho)
#endif
    end do
    
  end subroutine from_pseudopotential
end subroutine lcao_dens

subroutine lcao_wf(sys, h)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(in) :: h
  
  integer, parameter :: orbs_local = 2

  integer :: i, ik, n1, n2, i1, i2, l1, l2, lm1, lm2, d1, d2
  integer :: norbs, mode
  real(r8), allocatable :: hamilt(:,:), s(:,:), hpsi(:,:)
  real(r8), allocatable, target :: psis(:,:,:)
  real(r8), pointer :: psi1(:,:), psi2(:,:)

  ! variables for dsyev (LAPACK)
  character(len=1) :: jobz, uplo
  integer :: lwork, info
  real(r8), allocatable :: work(:), w(:)

  ! Counting...
  norbs = 0
  do i = 1, sys%natoms
    if(sys%atom(i)%spec%local) then
      norbs = norbs + orbs_local
    else
      norbs  = norbs + (sys%atom(i)%spec%ps%L_max + 1)**2
    end if
  end do
  if(sys%st%dim == 2) norbs = norbs*2 ! need twice the number of functions
  write(message(1),'(6x,i5,a)') norbs, ' functions to diagonalize.'
  call write_info(1)

  call oct_parse_int(C_string("LCAOMode"), 0, mode)
  if(mode < 0 .or. mode > 1) then
    message(1) = "LCAOMode not valid"
    message(2) = "LCAOMode = 0 (memory intensive) | 1 (cpu intensive)"
  end if

  ! Allocation of variables
  if(mode == 0) then
    allocate(psis(0:sys%m%np, sys%st%dim, norbs))
    psis = 0._r8
    n1 = 1
    do i1 = 1, sys%natoms
      do l1 = 0, sys%atom(i1)%spec%ps%L_max
        do lm1 = -l1, l1
          do d1 = 1, sys%st%dim
            call get_wf(i1, l1, lm1, d1, psis(:,:, n1))
            n1 = n1 + 1
          end do
        end do
      end do
    end do
  else
    allocate(psi1(0:sys%m%np, sys%st%dim), psi2(0:sys%m%np, sys%st%dim))
  end if

  ! Hamiltonian and overlap matrices, etc...
  allocate(hpsi(sys%m%np, sys%st%dim))
  allocate(hamilt(norbs, norbs), s(norbs, norbs))
  hamilt = 0.0_r8; s = 0.0_r8

  ik_loop : do ik = 1, sys%st%nik
    n1 = 1
    atoms1_loop: do i1 = 1, sys%natoms
      l1_loop: do l1 = 0, sys%atom(i1)%spec%ps%L_max
        lm1_loop: do lm1 = -l1, l1
          d1_loop: do d1 = 1, sys%st%dim
            if(mode == 0) then
              psi1 => psis(:,:, n1)
            else
              call get_wf(i1, l1, lm1, d1, psi1)
            end if

            call dhpsi(h, sys, ik, psi1, hpsi)
            
            n2 = 1
            atoms2_loop: do i2 = 1, sys%natoms
              l2_loop: do l2 = 0, sys%atom(i2)%spec%ps%L_max
                lm2_loop: do lm2 = -l2, l2
                  d2_loop: do d2 = 1, sys%st%dim
                    if(mode == 0) then
                      psi2 => psis(:,:, n2)
                    else
                      call get_wf(i2, l2, lm2, d2, psi2)
                    end if
                    
                    hamilt(n1, n2) = dstates_dotp(sys%m, sys%st%dim, hpsi, psi2(1:,:))
                    s(n1, n2) = dstates_dotp(sys%m, sys%st%dim, psi1(1:,:), psi2(1:,:))
                    hamilt(n2, n1) = hamilt(n1, n2)
                    s(n2, n1) = s(n1, n2)

                    if(n1 == n2 .and. conf%verbose >= 999) then
                      if(abs(s(n1, n2) - 1.0_r8) > 0.25_r8) then
                        write(message(1),'(a,i4,a,i2,a,i2,a,a,f12.6)') &
                             'Pseudo-wave function ', n1, ',', l1, ',', lm1,' out of box:', &
                             '|Phi|^2 = ', s(n1, n2)
                        call write_warning(1)
                      end if
                    end if
                    n2 = n2 + 1
                    if(n2 > n1) exit atoms2_loop
                  end do d2_loop
                end do lm2_loop
              end do l2_loop
            end do atoms2_loop

            n1 = n1 + 1
          end do d1_loop
        end do lm1_loop
      end do l1_loop
    end do atoms1_loop
    
    ! Setting variables for dsyev
    jobz = 'v'
    uplo = 'u'
    lwork = 3*norbs - 1
    allocate(work(lwork), w(norbs))

    call dsygv(1, jobz, uplo, norbs, hamilt, norbs, s, norbs, w, work, lwork, info)
    if(info.ne.0) then
      write(message(1),'(a,i5)') 'LAPACK "dsygv" returned error code ', info
      call write_fatal(1)
    endif
    deallocate(work, w)

    sys%st%R_FUNC(psi)(:,:,:, ik) = 0.0_r8

    n1 = 1
    do i1 = 1, sys%natoms
      do l1 = 0, sys%atom(i1)%spec%ps%L_max
        do lm1 = -l1, l1
          do d1 = 1, sys%st%dim
            if(mode == 0) then
              psi1 => psis(:,:, n1)
            else
              call get_wf(i1, l1, lm1, d1, psi1)
            end if

            do n2 = 1, sys%st%nst
              sys%st%R_FUNC(psi) (:,:, n2, ik) = sys%st%R_FUNC(psi) (:,:, n2, ik) + &
                   hamilt(n1, n2)*psi1(:,:)
            enddo
            n1 = n1 + 1
          end do
        end do
      end do
    end do

  enddo ik_loop

  deallocate(hamilt, s, hpsi)
  if(mode == 0) then
    deallocate(psis)
  else
    deallocate(psi1, psi2)
  end if

contains

  subroutine get_wf(i, l, lm, d, psi)
    integer, intent(in)   :: i, l, lm, d
    real(r8), intent(out) :: psi(0:sys%m%np, sys%st%dim)
    
    integer :: j
    real(r8) :: x(3), a(3), r, p, ylm, g(3)
    type(spline_type), pointer :: s
    
    a = sys%atom(i)%x

    psi(0,:) = 0.0_r8
    if(sys%atom(i)%spec%local) then
      ! add a couple of harmonic oscilator functions
    else
      s => sys%atom(i)%spec%ps%Ur(l)
      do j = 1, sys%m%np
        call mesh_r(sys%m, j, r, x=x, a=a)
        r = r/2._r8 ! enlarge wfs ??????
        p = splint(s, r)
        ylm = oct_ylm(x(1), x(2), x(3), l, lm)
        if(r > 0._r8) then
          psi(j, d) = p * ylm * r**l
        end if
      end do
    end if

    ! if spin channels are mixed we have to be careful
    if(sys%st%dim == 2) then
      if(d2 == 1) then
        psi(:, 2) = 0._r8
      else
        psi(:, 1) = 0._r8
      end if
    end if
    
  end subroutine get_wf
end subroutine lcao_wf

end module lcao
