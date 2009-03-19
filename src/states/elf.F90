!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: states.F90 2515 2006-10-24 17:13:30Z acastro $

#include "global.h"

module elf_m
  use datasets_m
  use derivatives_m
  use fourier_space_m
  use io_m
  use loct_parser_m
  use cube_function_m
  use global_m
  use grid_m
  use mesh_m
  use messages_m
  use mpi_m
  use profiling_m
  use states_m
  use states_dim_m

  implicit none

  private
  public :: elf_init,               &
            elf_calc,               &
            elf_calc_fs

  logical :: with_current_term = .true.

contains

  subroutine elf_init
    !%Variable ElfWithCurrentTerm
    !%Type logical
    !%Default true
    !%Section Output
    !%Description
    !% The ELF, when calculated for complex wave functions, should contain
    !% a term dependent on the current. This term is properly calculated by
    !% default; however, for research purposes it may be useful not to add it.
    !% If this feature proves to be useless, this option should go away.
    !%End
    call loct_parse_logical(datasets_check('ElfWithCurrentTerm'), .true., with_current_term)
  end subroutine elf_init

  ! ---------------------------------------------------------
  ! (time-dependent) electron localization function, (TD)ELF.
  ! ---------------------------------------------------------
  subroutine elf_calc(st, gr, elf, de)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    FLOAT,            intent(inout) :: elf(:,:) ! elf(NP, 1) if st%d%ispin = 1, elf(NP, 3) otherwise.
                                                ! On output, it should contain the global ELF if st%d%ispin = 1,
                                                ! otherwise elf(:, 3) contains the global elf, and 
                                                ! elf(:, 1) and elf(:, 2) the spin resolved ELF.
    FLOAT,   optional, intent(inout):: de(:,:)

    FLOAT :: f, D0, dens
    integer :: i, is, nelfs
    FLOAT, allocatable :: rho(:,:), grho(:,:,:), jj(:,:,:)
    FLOAT, allocatable :: kappa(:,:)

    FLOAT, parameter :: dmin = CNST(1e-10)

    call push_sub('states.states_calc_elf')

    ! We may or may not want the total elf. 
    ! If we want it, the argument elf should have three components.c
    nelfs = size(elf, 2)

    ALLOCATE(rho(NP, st%d%nspin), NP)
    ALLOCATE(kappa(NP, st%d%nspin), NP)
    rho = M_ZERO
    kappa = M_ZERO
    call states_calc_dens(st, NP, rho)

    ALLOCATE(grho(NP, gr%mesh%sb%dim, st%d%nspin), NP*gr%mesh%sb%dim*st%d%nspin)
    ALLOCATE(  jj(NP, gr%mesh%sb%dim, st%d%nspin), NP*gr%mesh%sb%dim*st%d%nspin)

    call states_calc_tau_jp_gn(gr, st, kappa, jj, grho)

    ! spin-dependent quantities
    if(st%d%ispin == UNPOLARIZED) then
      rho = rho/M_TWO
      kappa = kappa/M_TWO
      jj = jj/M_TWO
      grho = grho/M_TWO
    end if
    
    if(.not.with_current_term) jj = M_ZERO

    ! kapp will contain rho * D
    do_is: do is = 1, st%d%nspin
      do i = 1, NP
        kappa(i, is) = kappa(i, is)*rho(i, is)        &    ! + tau * rho
          - M_FOURTH*sum(grho(i, 1:gr%mesh%sb%dim, is)**2)      &    ! - | nabla rho |^2 / 4
          - sum(jj(i, 1:gr%mesh%sb%dim, is)**2)                      ! - j^2
      end do

      deallocate(grho, jj)   ! these are no longer needed
    
      ! pass this information to the caller if requested
      if(present(de)) de(1:NP,is) = kappa(1:NP,is)

    end do do_is

    select case(gr%sb%dim)
      case(3); f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
      case(2); f = M_TWO * M_Pi
    end select

    select case(st%d%ispin)
    case(UNPOLARIZED)
      do i = 1, NP
        if(rho(i, 1) >= dmin) then
          select case(gr%sb%dim)
            case(3); D0 = f * rho(i, 1)**(M_EIGHT/M_THREE)
            case(2); D0 = f * rho(i, 1)**3
          end select
          elf(i, 1) = D0*D0/(D0*D0 + kappa(i, 1)**2)
        else
          elf(i, 1) = M_ZERO
        endif
      end do

    case(SPIN_POLARIZED, SPINORS)
      if(nelfs .eq. 3) then
        do i = 1, NP
          dens = rho(i, 1) + rho(i, 2)
          if( dens >= dmin ) then
            select case(gr%sb%dim)
              case(3); D0 = f * dens ** (M_FIVE/M_THREE) * rho(i, 1) * rho(i, 2)
              case(2); D0 = f * dens ** 2 * rho(i, 1) * rho(i, 2)
            end select
            elf(i, 3) = D0*D0/(D0*D0 + (kappa(i, 1)*rho(i, 2) + kappa(i,2)*rho(i, 1))**2)
          else
            elf(i, 3) = M_ZERO
          endif
        end do
      end if
      do i = 1, NP
        do is = 1, st%d%spin_channels
          if(rho(i, is) >= dmin) then
            select case(gr%sb%dim)
              case(3); D0 = f * rho(i, is)**(M_EIGHT/M_THREE)
              case(2); D0 = f * rho(i, is)**3
            end select
            elf(i, is) = D0*D0/(D0*D0 + kappa(i,is)**2)
          else
            elf(i, is) = M_ZERO
          endif
        end do
      end do
    end select

    deallocate(rho, kappa)

    call pop_sub()
  end subroutine elf_calc


  ! ---------------------------------------------------------
  ! ELF function in Fourier space. Not tested.
  ! ---------------------------------------------------------
  subroutine elf_calc_fs(st, gr, elf, de)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    FLOAT,             intent(inout):: elf(:,:)
    FLOAT,   optional, intent(inout):: de(:,:)

    FLOAT :: f, d, s
    integer :: i, is, ik, ist, idim
    CMPLX, allocatable :: psi_fs(:), gpsi(:,:)
    FLOAT, allocatable :: r(:), gradr(:,:), j(:,:)
    type(dcf_t) :: dcf_tmp
    type(zcf_t) :: zcf_tmp

    FLOAT, parameter :: dmin = CNST(1e-10)
    
    ! single or double occupancy
    if(st%d%nspin == 1) then
      s = M_TWO
    else
      s = M_ONE
    end if
 
    if (st%wfs_type == M_REAL) then
      call dcf_new(gr%mesh%idx%ll, dcf_tmp)
      call dcf_fft_init(dcf_tmp, gr%sb)
    else
      call zcf_new(gr%mesh%idx%ll, zcf_tmp)
      call zcf_fft_init(zcf_tmp, gr%sb)
    end if

    do_is: do is = 1, st%d%nspin
      ALLOCATE(    r(NP),       NP)
      ALLOCATE(gradr(NP, gr%mesh%sb%dim), NP*gr%mesh%sb%dim)
      ALLOCATE(    j(NP, gr%mesh%sb%dim), NP*gr%mesh%sb%dim)
      r = M_ZERO; gradr = M_ZERO; j  = M_ZERO

      elf(1:NP,is) = M_ZERO

      ALLOCATE(psi_fs(gr%mesh%np_part),  gr%mesh%np_part)
      ALLOCATE(gpsi  (NP, gr%mesh%sb%dim), NP*gr%mesh%sb%dim)
      do ik = is, st%d%nik, st%d%nspin
        do ist = 1, st%nst
          do idim = 1, st%d%dim
            
            if (st%wfs_type == M_REAL) then
              call dmf2mf_RS2FS(gr%mesh, st%dpsi(:, idim, ist, ik), psi_fs(:), dcf_tmp)
              call zderivatives_grad(gr%der, psi_fs(:), gpsi)
              do i = 1, gr%mesh%sb%dim
                gpsi(:,i) = gpsi(:,i) * gr%mesh%h(i)**2 * real(dcf_tmp%n(i), REAL_PRECISION) / (M_TWO*M_PI)
              end do
            else
              call zmf2mf_RS2FS(gr%mesh, st%zpsi(:, idim, ist, ik), psi_fs(:), zcf_tmp)
              call zderivatives_grad(gr%der, psi_fs(:), gpsi)
              do i = 1, gr%mesh%sb%dim
                gpsi(:,i) = gpsi(:,i) * gr%mesh%h(i)**2 * real(zcf_tmp%n(i), REAL_PRECISION) / (M_TWO*M_PI)
              end do
            end if

            r(:) = r(:) + st%d%kweights(ik)*st%occ(ist, ik) * abs(psi_fs(:))**2
            do i = 1, gr%mesh%sb%dim
              gradr(:,i) = gradr(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   M_TWO * real(conjg(psi_fs(:))*gpsi(:,i))
              j (:,i) =  j(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   aimag(conjg(psi_fs(:))*gpsi(:,i))
            end do

            do i = 1, NP
              if(r(i) >= dmin) then
                elf(i,is) = elf(i,is) + st%d%kweights(ik)*st%occ(ist, ik)/s * &
                     sum(abs(gpsi(i, 1:gr%mesh%sb%dim))**2)
              end if
            end do
          end do
        end do
      end do
      deallocate(psi_fs, gpsi)

      do i = 1, NP
        if(r(i) >= dmin) then
          elf(i,is) = elf(i,is) - (M_FOURTH*sum(gradr(i, 1:gr%mesh%sb%dim)**2) + sum(j(i, 1:gr%mesh%sb%dim)**2))/(s*r(i))
        end if
      end do
    
      if(present(de)) de(1:NP,is)=elf(1:NP,is)

      ! normalization
      f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
      do i = 1, NP
        if(abs(r(i)) >= dmin) then
          d    = f*(r(i)/s)**(M_FIVE/M_THREE)
          elf(i,is) = M_ONE/(M_ONE + (elf(i,is)/d)**2)
        else
          elf(i,is) = M_ZERO
        end if
      end do

      deallocate(r, gradr, j)

    end do do_is

    if (st%wfs_type == M_REAL) then
      call dcf_free(dcf_tmp)
    else
      call zcf_free(zcf_tmp)
    end if

  contains

    subroutine dmf2mf_RS2FS(m, fin, fout, c)
      type(mesh_t),  intent(in)    :: m
      FLOAT,         intent(in)    :: fin(:)
      CMPLX,         intent(out)   :: fout(:)
      type(dcf_t),   intent(inout) :: c
    
      call dcf_alloc_RS(c)
      call dcf_alloc_FS(c)
      call dmesh_to_cube(m, fin, c)
      call dcf_RS2FS(c)
      call dfourier_to_mesh(m, c, fout)
      call dcf_free_RS(c)
      call dcf_free_FS(c)
    end subroutine dmf2mf_RS2FS

    subroutine zmf2mf_RS2FS(m, fin, fout, c)
      type(mesh_t),  intent(in)    :: m
      CMPLX,         intent(in)    :: fin(:)
      CMPLX,         intent(out)   :: fout(:)
      type(zcf_t),   intent(inout) :: c
    
      call zcf_alloc_RS(c)
      call zcf_alloc_FS(c)
      call zmesh_to_cube(m, fin, c)
      call zcf_RS2FS(c)
      call zfourier_to_mesh(m, c, fout)
      call zcf_free_RS(c)
      call zcf_free_FS(c)
    end subroutine zmf2mf_RS2FS

  end subroutine elf_calc_fs

end module elf_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
