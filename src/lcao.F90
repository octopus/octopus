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

  private
  public :: lcao_dens, lcao_init, lcao_wf, lcao_end

type lcao_type
  integer           :: mode
  integer           :: dim
  real(r8), pointer :: psis(:,:,:)
  real(r8), pointer :: hamilt(:,:), k_plus_psv(:,:,:), s(:,:,:)
  logical, pointer  :: atoml(:,:)
end type

type(lcao_type) :: lcao_data

integer, parameter :: MEM_INTENSIVE = 0, &
                      CPU_INTENSIVE = 1

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

  sub_name = 'lcao_dens'; call push_sub()
  
  rho = 0._r8
  do ia = 1, sys%natoms
    a => sys%atom(ia) ! shortcuts
    s => a%spec
    
    select case(s%label(1:5))
    case('jelli', 'point')
      call from_jellium(sys%m, rho(:, 1))
    case('usdef')
      call from_userdef(sys%m, rho(:, 1))
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

  call pop_sub()
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
  
  subroutine from_userdef(m, rho)
    type(mesh_type), intent(in) :: m
    real(r8), intent(inout) :: rho(m%np)

    integer :: i

    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
      rho(i) = rho(i) + real(s%Z_val, r8)/(m%np*m%vol_pp)
    end do
  end subroutine from_userdef

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
      do l = 0 , s%ps%L_max_occ
        if(r >= r_small) then
          psi = splint(s%ps%Ur(l), r)
          zel = min(zval, 2.0_r8*(2*l+1))
          zval = max(0._r8, zval - zel)
          rho(i) = rho(i) + zel*psi*psi*(r**(2*l))/(4*M_PI)
        end if
      end do
#elif defined(ONE_D)
      rho(i) = rho(i) + s%z_val*exp(-r**2)/sqrt(m_pi)
      !call R_FUNC(calcdens)(sys%st, m%np, rho)
#endif
    end do
    
  end subroutine from_pseudopotential
end subroutine lcao_dens

subroutine lcao_init(sys, h)
  type(system_type), intent(IN)      :: sys
  type(hamiltonian_type), intent(IN) :: h

  integer :: norbs, i, a, ik, n1, i1, l1, lm, lm1, d1, n2, i2, l2, lm2, d2
  integer, parameter :: orbs_local = 2

  real(r8), allocatable :: psi1(:,:), psi2(:,:), hpsi(:,:)
  real(r8) :: s

  real(r8) :: uVpsi

  sub_name = 'lcao_init'; call push_sub

  ! Counting -- remove the l-components which are essentially out of bounds.
  allocate(lcao_data%atoml(sys%natoms, 0:3))
  lcao_data%atoml = .true.
  allocate(psi1(0:sys%m%np, sys%st%dim))
  norbs = 0
  atoms_loop: do i1 = 1, sys%natoms
      l_loop: do l1 = 0, sys%atom(i1)%spec%ps%L_max_occ
        lm_loop: do lm1 = -l1, l1
          d_loop: do d1 = 1, sys%st%dim
            call get_wf(sys, i1, l1, lm1, d1, psi1)
            s = dstates_dotp(sys%m, sys%st%dim, psi1(1:,:), psi1(1:,:))
            if(s < 0.1_r8) then
              lcao_data%atoml(i1, l1) = .false.
              norbs = norbs - (lm1 + l1)*sys%st%dim 
              cycle l_loop
            endif
            norbs = norbs + 1
          end do d_loop
        end do lm_loop
      end do l_loop
  enddo atoms_loop 
  deallocate(psi1)
  lcao_data%dim = norbs
  write(message(1), '(a,i6)') 'Info: LCAO basis dimension: ', lcao_data%dim
  call write_info(1)

  ! Gets the mode
  call oct_parse_int(C_string("LCAOMode"), 0, lcao_data%mode)
  if(lcao_data%mode < 0 .or. lcao_data%mode > 1) then
    message(1) = "LCAOMode not valid"
    message(2) = "LCAOMode = 0 (memory intensive) | 1 (cpu intensive)"
  end if

  ! Gets the wave-functions
  if(lcao_data%mode == 0) then
    allocate(lcao_data%psis(0:sys%m%np, sys%st%dim, norbs))
    lcao_data%psis = 0._r8
    n1 = 1
    do i1 = 1, sys%natoms
      do l1 = 0, sys%atom(i1)%spec%ps%L_max_occ
        if(.not. lcao_data%atoml(i1, l1)) cycle
        do lm1 = -l1, l1
          do d1 = 1, sys%st%dim
            call get_wf(sys, i1, l1, lm1, d1, lcao_data%psis(:,:, n1))
            n1 = n1 + 1
          end do
        end do
      end do
    end do
  end if

  ! Allocation of variables
  allocate(lcao_data%hamilt(norbs, norbs), lcao_data%s(sys%st%nik, norbs, norbs), &
           lcao_data%k_plus_psv(sys%st%nik, norbs, norbs))

  !Fixes, once and for all, the overlap matrix and kinetic + pseudopotential matrixes..
  if(lcao_data%mode == 1) allocate(psi1(0:sys%m%np, sys%st%dim), psi2(0:sys%m%np, sys%st%dim))
  allocate(hpsi(sys%m%np, sys%st%dim))
  hpsi = 0.0_r8
  ik_loop : do ik = 1, sys%st%nik
    n1 = 1
    atoms1_loop: do i1 = 1, sys%natoms
      l1_loop: do l1 = 0, sys%atom(i1)%spec%ps%L_max_occ
        if(.not. lcao_data%atoml(i1, l1))  cycle
        lm1_loop: do lm1 = -l1, l1
          d1_loop: do d1 = 1, sys%st%dim
            select case(lcao_data%mode)
            case(CPU_INTENSIVE)
              call get_wf(sys, i1, l1, lm1, d1, psi1)
              call R_FUNC(mesh_derivatives) (sys%m, psi1(:, d1), lapl=hpsi(:, d1))
              hpsi = -hpsi/2.0_r8
              hpsi(:, d1) = hpsi(:, d1) + h%vpsl*psi1(1:,d1)
              call dvnlpsi(sys, psi1(:,:), hpsi)
            case(MEM_INTENSIVE)
              call R_FUNC(mesh_derivatives) (sys%m, lcao_data%psis(:,d1,n1), lapl=hpsi(:, d1))
              hpsi = -hpsi/2.0_r8
              hpsi(:, d1) = hpsi(:, d1) + h%vpsl*lcao_data%psis(1:,d1,n1)
              call dvnlpsi(sys, lcao_data%psis(:,:,n1), hpsi)
            end select

            n2 = 1
            atoms2_loop: do i2 = 1, sys%natoms
              l2_loop: do l2 = 0, sys%atom(i2)%spec%ps%L_max_occ
                if(.not. lcao_data%atoml(i1, l1)) cycle
                lm2_loop: do lm2 = -l2, l2
                  d2_loop: do d2 = 1, sys%st%dim
                    select case(lcao_data%mode)
                    case(CPU_INTENSIVE)
                      call get_wf(sys, i2, l2, lm2, d2, psi2)
                      lcao_data%s(ik, n1, n2) = dstates_dotp(sys%m, sys%st%dim, psi1(1:,:), psi2(1:,:))
                      lcao_data%k_plus_psv(ik, n1, n2) = dstates_dotp(sys%m, sys%st%dim, hpsi, psi2(1:, :))
                    case(MEM_INTENSIVE)
                      lcao_data%s(ik, n1, n2) = dstates_dotp(sys%m, sys%st%dim, lcao_data%psis(1:,:,n1), lcao_data%psis(1:,:,n2))
                      lcao_data%k_plus_psv(ik, n1, n2) = dstates_dotp(sys%m, sys%st%dim, hpsi, lcao_data%psis(1:, :, n2))
                    end select
                    lcao_data%s(ik, n2, n1) = lcao_data%s(ik, n1, n2)
                    lcao_data%k_plus_psv(ik, n2, n1) = lcao_data%k_plus_psv(ik, n1, n2)

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
  enddo ik_loop
  if(lcao_data%mode==1) deallocate(psi1, psi2)
  deallocate(hpsi)

  call pop_sub()
end subroutine lcao_init

subroutine lcao_end
  sub_name = 'lcao_end'; call push_sub()

  if(lcao_data%mode == MEM_INTENSIVE) deallocate(lcao_data%psis)
  deallocate(lcao_data%hamilt, lcao_data%s, lcao_data%atoml)

  call pop_sub()
end subroutine lcao_end

subroutine lcao_wf(sys, h)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(in) :: h
  
  integer, parameter :: orbs_local = 2

  integer :: a, idim, i, lm, ik, n1, n2, i1, i2, l1, l2, lm1, lm2, d1, d2
  integer :: norbs, mode
  real(r8), allocatable :: hpsi(:,:), psi1(:,:), psi2(:,:)
  real(r8), allocatable :: s(:,:)
  real(r8) :: uvpsi

  ! variables for dsyev (LAPACK)
  character(len=1) :: jobz, uplo
  integer :: lwork, info
  real(r8), allocatable :: work(:), w(:)

  sub_name = 'lcao_wf'; call push_sub()

  norbs = lcao_data%dim
  mode  = lcao_data%mode

  ! Allocation of variables
  if(mode == 1) allocate(psi1(0:sys%m%np, sys%st%dim), psi2(0:sys%m%np, sys%st%dim))

  ! Hamiltonian and overlap matrices, etc...
  allocate(hpsi(sys%m%np, sys%st%dim))

  ik_loop : do ik = 1, sys%st%nik
    n1 = 1
    atoms1_loop: do i1 = 1, sys%natoms
      l1_loop: do l1 = 0, sys%atom(i1)%spec%ps%L_max_occ
        if(.not.lcao_data%atoml(i1, l1)) cycle l1_loop
        lm1_loop: do lm1 = -l1, l1
          d1_loop: do d1 = 1, sys%st%dim
            select case(lcao_data%mode)
            case(MEM_INTENSIVE)
              hpsi(:, d1) = h%vhartree*lcao_data%psis(1:, d1, n1)
              hpsi(:, d1) = hpsi(:, d1) + h%Vxc(:, d1)*lcao_data%psis(1:, d1, n1)
            case(CPU_INTENSIVE)
              call get_wf(sys, i1, l1, lm1, d1, psi1)
              hpsi(:, d1) = h%vhartree*psi1(1:, d1)
              hpsi(:, d1) = hpsi(:, d1) + h%Vxc(:, d1)*psi1(1:, d1)
            end select
            
            n2 = 1
            atoms2_loop: do i2 = 1, sys%natoms
              l2_loop: do l2 = 0, sys%atom(i2)%spec%ps%L_max_occ
                if(.not.lcao_data%atoml(i2, l2)) cycle l2_loop
                lm2_loop: do lm2 = -l2, l2
                  d2_loop: do d2 = 1, sys%st%dim
                    select case(lcao_data%mode)
                    case(MEM_INTENSIVE)
                      lcao_data%hamilt(n1, n2) = lcao_data%k_plus_psv(ik, n1, n2) + & 
                                                 dstates_dotp(sys%m, sys%st%dim, hpsi, lcao_data%psis(1:, : ,n2))
                    case(CPU_INTENSIVE)
                      call get_wf(sys, i2, l2, lm2, d2, psi2)
                      lcao_data%hamilt(n1, n2) = lcao_data%k_plus_psv(ik, n1, n2) + &
                                                 dstates_dotp(sys%m, sys%st%dim, hpsi, psi2(1:,:))
                    end select

                    lcao_data%hamilt(n2, n1) = lcao_data%hamilt(n1, n2)
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

!!$    ! For debugging, print out the hamiltonian and overlap matrixes...
!!$    write(*,'(a)') '  LCAO debugging:'
!!$    write(*,'(a)') '    Hamiltonian matrix:'
!!$    do n1=1, norbs
!!$       do n2=1, norbs
!!$          write(*,'(4x,f12.6)', advance='no') lcao_data%hamilt(n1,n2)/units_out%energy%factor
!!$       enddo
!!$       write(*,*)
!!$    enddo
!!$    write(*,'(a)') '    Overlap matrix:'
!!$    do n1=1, norbs
!!$       do n2=1, norbs
!!$          write(*,'(4x,f12.6)', advance='no') lcao_data%s(ik, n1,n2)
!!$       enddo
!!$       write(*,*)
!!$    enddo
!!$    ! End of debugging

    ! Setting variables for dsyev
    jobz = 'v'
    uplo = 'u'
    lwork = 3*norbs - 1
    allocate(work(lwork), w(norbs))
    allocate(s(norbs, norbs))
    s(1:norbs, 1:norbs) = lcao_data%s(ik, 1:norbs, 1:norbs)

    call dsygv(1, jobz, uplo, norbs, lcao_data%hamilt, norbs, s, norbs, w, work, lwork, info)
    if(info.ne.0) then
      write(message(1),'(a,i5)') 'LAPACK "dsygv" returned error code ', info
      call write_fatal(1)
    endif
    deallocate(work, w)

    sys%st%R_FUNC(psi)(:,:,:, ik) = 0.0_r8

    n1 = 1
    do i1 = 1, sys%natoms
      do l1 = 0, sys%atom(i1)%spec%ps%L_max_occ
        if(.not.lcao_data%atoml(i1, l1)) cycle
        do lm1 = -l1, l1
          do d1 = 1, sys%st%dim
            if(mode == 0) then
              do n2 = 1, sys%st%nst
                 sys%st%R_FUNC(psi) (:,:, n2, ik) = sys%st%R_FUNC(psi) (:,:, n2, ik) + &
                   lcao_data%hamilt(n1, n2)*lcao_data%psis(:,:,n1)
              enddo
            else
              call get_wf(sys, i1, l1, lm1, d1, psi1)
              do n2 = 1, sys%st%nst
                 sys%st%R_FUNC(psi) (:,:, n2, ik) = sys%st%R_FUNC(psi) (:,:, n2, ik) + &
                    lcao_data%hamilt(n1, n2)*psi1(:,:)
              enddo
            end if

            n1 = n1 + 1
          end do
        end do
      end do
    end do

  enddo ik_loop

  deallocate(hpsi); if(mode == 1) deallocate(psi1, psi2)

  call pop_sub()
end subroutine lcao_wf

subroutine get_wf(sys, i, l, lm, d, psi)
  type(system_type), intent(IN) :: sys
  integer, intent(in)   :: i, l, lm, d
  real(r8), intent(out) :: psi(0:sys%m%np, sys%st%dim)
    
  integer :: j, d2
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
      p = splint(s, r)
      ylm = oct_ylm(x(1), x(2), x(3), l, lm)
      if(r > 0._r8) then
        psi(j, d) = p * ylm * r**l
      end if
    end do
  end if

  ! if spin channels are mixed we have to be careful
  if(sys%st%dim == 2) then
    if(d == 1) then
      psi(:, 2) = 0._r8
    else
      psi(:, 1) = 0._r8
    end if
  end if
    
end subroutine get_wf

end module lcao
