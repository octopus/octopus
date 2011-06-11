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
  use cube_function_m
  use density_m
  use datasets_m
  use derivatives_m
  use fourier_space_m
  use global_m
  use grid_m
  use io_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
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
    PUSH_SUB(elf_init)

    !%Variable ELFWithCurrentTerm
    !%Type logical
    !%Default true
    !%Section Output
    !%Description
    !% The ELF, when calculated for complex wavefunctions, should contain
    !% a term dependent on the current. This term is properly calculated by
    !% default; however, for research purposes it may be useful not to add it.
    !% If this feature proves to be useless, this option should go away.
    !%End
    call parse_logical(datasets_check('ELFWithCurrentTerm'), .true., with_current_term)

    POP_SUB(elf_init)
  end subroutine elf_init

  ! ---------------------------------------------------------
  ! (time-dependent) electron localization function, (TD)ELF.
  ! ---------------------------------------------------------
  subroutine elf_calc(st, gr, elf, de)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    FLOAT,            intent(inout) :: elf(:,:) ! elf(gr%mesh%np, 1) if st%d%ispin = 1, elf(gr%mesh%np, 3) otherwise.
                                                ! On output, it should contain the global ELF if st%d%ispin = 1,
                                                ! otherwise elf(:, 3) contains the global ELF, and 
                                                ! elf(:, 1) and elf(:, 2) the spin-resolved ELF.
    FLOAT,  optional, intent(inout):: de(:,:)

    FLOAT :: factor, D0, dens
    integer :: ip, is, nelfs
    FLOAT, allocatable :: rho(:,:), grho(:,:,:), jj(:,:,:)
    FLOAT, allocatable :: kappa(:,:)
    FLOAT, parameter :: dmin = CNST(1e-10)

    PUSH_SUB(elf_calc)

    ASSERT(gr%sb%dim == 2 .or. gr%sb%dim == 3)

    ! We may or may not want the total ELF. 
    ! If we want it, the argument elf should have three components.
    nelfs = size(elf, 2)

    SAFE_ALLOCATE(  rho(1:gr%mesh%np, 1:st%d%nspin))
    SAFE_ALLOCATE(kappa(1:gr%mesh%np, 1:st%d%nspin))
    rho = M_ZERO
    kappa = M_ZERO
    call density_calc(st, gr, rho)

    SAFE_ALLOCATE(grho(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%nspin))
    SAFE_ALLOCATE(  jj(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%nspin))

    call states_calc_quantities(gr%der, st, kinetic_energy_density = kappa, paramagnetic_current = jj, density_gradient = grho)

    ! spin-dependent quantities
    if(st%d%ispin == UNPOLARIZED) then
      rho = rho/M_TWO
      kappa = kappa/M_TWO
      jj = jj/M_TWO
      grho = grho/M_TWO
    end if
    
    if(.not.with_current_term) jj = M_ZERO

    ! kappa will contain rho * D
    do_is: do is = 1, st%d%nspin
      do ip = 1, gr%mesh%np
        kappa(ip, is) = kappa(ip, is)*rho(ip, is)        &    ! + tau * rho
          - M_FOURTH*sum(grho(ip, 1:gr%mesh%sb%dim, is)**2)      &    ! - | nabla rho |^2 / 4
          - sum(jj(ip, 1:gr%mesh%sb%dim, is)**2)                      ! - j^2
      end do

      ! pass this information to the caller if requested
      if(present(de)) de(1:gr%mesh%np,is) = kappa(1:gr%mesh%np,is)

    end do do_is

    SAFE_DEALLOCATE_A(grho)
    SAFE_DEALLOCATE_A(jj)   ! these are no longer needed

    select case(gr%sb%dim)
      case(3); factor = M_THREE / M_FIVE * (M_SIX * M_PI**2)**M_TWOTHIRD
      case(2); factor = M_TWO * M_Pi
    end select

    select case(st%d%ispin)
    case(UNPOLARIZED)
      do ip = 1, gr%mesh%np
        if(rho(ip, 1) >= dmin) then
          select case(gr%sb%dim)
            case(3); D0 = factor * rho(ip, 1)**(M_EIGHT / M_THREE)
            case(2); D0 = factor * rho(ip, 1)**3
          end select
          elf(ip, 1) = D0 * D0 / (D0 * D0 + kappa(ip, 1)**2)
        else
          elf(ip, 1) = M_ZERO
        endif
      end do

    case(SPIN_POLARIZED, SPINORS)
      if(nelfs .eq. 3) then
        do ip = 1, gr%mesh%np
          dens = rho(ip, 1) + rho(ip, 2)
          if( dens >= dmin ) then
            select case(gr%sb%dim)
              case(3); D0 = factor * dens ** (M_FIVE/M_THREE) * rho(ip, 1) * rho(ip, 2)
              case(2); D0 = factor * dens ** 2 * rho(ip, 1) * rho(ip, 2)
            end select
            elf(ip, 3) = D0 * D0 / (D0 * D0 + (kappa(ip, 1) * rho(ip, 2) + kappa(ip,2) * rho(ip, 1))**2)
          else
            elf(ip, 3) = M_ZERO
          endif
        end do
      end if
      do ip = 1, gr%mesh%np
        do is = 1, st%d%spin_channels
          if(rho(ip, is) >= dmin) then
            select case(gr%sb%dim)
              case(3); D0 = factor * rho(ip, is)**(M_EIGHT/M_THREE)
              case(2); D0 = factor * rho(ip, is)**3
            end select
            elf(ip, is) = D0 * D0 / (D0 * D0 + kappa(ip,is)**2)
          else
            elf(ip, is) = M_ZERO
          endif
        end do
      end do
    end select

    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(kappa)

    POP_SUB(elf_calc)
  end subroutine elf_calc


  ! ---------------------------------------------------------
  ! ELF function in Fourier space. Not tested.
  ! ---------------------------------------------------------
  subroutine elf_calc_fs(st, gr, elf, de)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    FLOAT,             intent(inout):: elf(:,:)
    FLOAT,   optional, intent(inout):: de(:,:)

    FLOAT :: factor, dd, sp
    integer :: ip, is, ik, ist, idim, idir
    CMPLX, allocatable :: psi_fs(:), gpsi(:,:), zpsi(:)
    FLOAT, allocatable :: rr(:), gradr(:,:), jj(:,:), dpsi(:)
    type(cube_function_t) :: cube_function_tmp
    FLOAT, parameter :: dmin = CNST(1e-10)

    PUSH_SUB(elf_calc_fs)
    
    ! single or double occupancy
    if(st%d%nspin == 1) then
      sp = M_TWO
    else
      sp = M_ONE
    end if
 

    call cube_function_init(cube_function_tmp, gr%mesh%idx%ll)

    if (states_are_real(st)) then
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np))
      call dcube_function_fft_init(cube_function_tmp, gr%sb)
    else
      SAFE_ALLOCATE(zpsi(1:gr%mesh%np))
      call zcube_function_fft_init(cube_function_tmp, gr%sb)
    end if

    do_is: do is = 1, st%d%nspin
      SAFE_ALLOCATE(   rr(1:gr%mesh%np))
      SAFE_ALLOCATE(gradr(1:gr%mesh%np, 1:gr%mesh%sb%dim))
      SAFE_ALLOCATE(   jj(1:gr%mesh%np, 1:gr%mesh%sb%dim))
      rr = M_ZERO; gradr = M_ZERO; jj = M_ZERO

      elf(1:gr%mesh%np,is) = M_ZERO

      SAFE_ALLOCATE(psi_fs(1:gr%mesh%np_part))
      SAFE_ALLOCATE(gpsi  (1:gr%mesh%np, gr%mesh%sb%dim))
      do ik = is, st%d%nik, st%d%nspin
        do ist = 1, st%nst
          do idim = 1, st%d%dim
            
            if (states_are_real(st)) then
              call states_get_state(st, gr%mesh, idim, ist, ik, dpsi)
              call dmf2mf_RS2FS(gr%mesh, dpsi, psi_fs(:), cube_function_tmp)
              call zderivatives_grad(gr%der, psi_fs(:), gpsi)
              do idir = 1, gr%mesh%sb%dim
                gpsi(:,idir) = gpsi(:, idir)*gr%mesh%spacing(idir)**2 * &
                     real(cube_function_tmp%n(idir), REAL_PRECISION)/(M_TWO*M_PI)
              end do
            else
              call states_get_state(st, gr%mesh, idim, ist, ik, zpsi)
              call zmf2mf_RS2FS(gr%mesh, zpsi, psi_fs(:), cube_function_tmp)
              call zderivatives_grad(gr%der, psi_fs(:), gpsi)
              do idir = 1, gr%mesh%sb%dim
                gpsi(:, idir) = gpsi(:, idir)*gr%mesh%spacing(idir)**2 * &
                     real(cube_function_tmp%n(idir), REAL_PRECISION)/(M_TWO*M_PI)
              end do
            end if

            rr(:) = rr(:) + st%d%kweights(ik)*st%occ(ist, ik) * abs(psi_fs(:))**2
            do idir = 1, gr%mesh%sb%dim
              gradr(:, idir) = gradr(:, idir) + st%d%kweights(ik) * st%occ(ist, ik) *  &
                   M_TWO * real(conjg(psi_fs(:))*gpsi(:, idir))
              jj(:, idir) = jj(:, idir) + st%d%kweights(ik) * st%occ(ist, ik) *  &
                   aimag(conjg(psi_fs(:)) * gpsi(:, idir))
            end do

            do ip = 1, gr%mesh%np
              if(rr(ip) >= dmin) then
                elf(ip, is) = elf(ip, is) + st%d%kweights(ik) * st%occ(ist, ik)/sp * &
                     sum(abs(gpsi(ip, 1:gr%mesh%sb%dim))**2)
              end if
            end do
          end do
        end do
      end do
      SAFE_DEALLOCATE_A(psi_fs)
      SAFE_DEALLOCATE_A(gpsi)

      do ip = 1, gr%mesh%np
        if(rr(ip) >= dmin) then
          elf(ip, is) = elf(ip, is) - &
            (M_FOURTH * sum(gradr(ip, 1:gr%mesh%sb%dim)**2) + sum(jj(ip, 1:gr%mesh%sb%dim)**2)) / (sp * rr(ip))
        end if
      end do
    
      if(present(de)) de(1:gr%mesh%np, is) = elf(1:gr%mesh%np, is)

      ! normalization
      factor = M_THREE / M_FIVE * (M_SIX * M_PI**2)**M_TWOTHIRD
      do ip = 1, gr%mesh%np
        if(abs(rr(ip)) >= dmin) then
          dd = factor * (rr(ip) / sp)**(M_FIVE / M_THREE)
          elf(ip, is) = M_ONE / (M_ONE + (elf(ip, is) / dd)**2)
        else
          elf(ip, is) = M_ZERO
        end if
      end do

      SAFE_DEALLOCATE_A(rr)
      SAFE_DEALLOCATE_A(gradr)
      SAFE_DEALLOCATE_A(jj)

    end do do_is

    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(zpsi)

    call cube_function_fft_end(cube_function_tmp)
    call cube_function_end(cube_function_tmp)

    POP_SUB(elf_calc_fs)

  contains

    ! ---------------------------------------------------------
    subroutine dmf2mf_RS2FS(mesh, fin, fout, cc)
      type(mesh_t),  intent(in)    :: mesh
      FLOAT,         intent(in)    :: fin(:)
      CMPLX,         intent(out)   :: fout(:)
      type(cube_function_t),   intent(inout) :: cc
    
      PUSH_SUB(elf_calc_fs.dmf2mf_RS2FS)

      call dcube_function_alloc_RS(cc)
      call dcube_function_alloc_FS(cc)
      call dmesh_to_cube(mesh, fin, cc)
      call dcube_function_RS2FS(cc)
      call dfourier_to_mesh(mesh, cc, fout)
      call dcube_function_free_RS(cc)
      call dcube_function_free_FS(cc)

      POP_SUB(elf_calc_fs.dmf2mf_RS2FS)
    end subroutine dmf2mf_RS2FS

    ! ---------------------------------------------------------
    subroutine zmf2mf_RS2FS(mesh, fin, fout, cc)
      type(mesh_t),  intent(in)    :: mesh
      CMPLX,         intent(in)    :: fin(:)
      CMPLX,         intent(out)   :: fout(:)
      type(cube_function_t),   intent(inout) :: cc
    
      PUSH_SUB(elf_calc_fs.zmf2mf_RS2FS)

      call zcube_function_alloc_RS(cc)
      call zcube_function_alloc_FS(cc)
      call zmesh_to_cube(mesh, fin, cc)
      call zcube_function_RS2FS(cc)
      call zfourier_to_mesh(mesh, cc, fout)
      call zcube_function_free_RS(cc)
      call zcube_function_free_FS(cc)

      POP_SUB(elf_calc_fs.zmf2mf_RS2FS)
    end subroutine zmf2mf_RS2FS

  end subroutine elf_calc_fs

end module elf_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
