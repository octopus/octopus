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
  use io_m
  use lib_oct_parser_m
  use cube_function_m
  use functions_m
  use global_m
  use grid_m
  use mesh_m
  use messages_m
  use mpi_m
  use states_m
  use basins_m
  use output_m

  implicit none

  private
  public :: elf_init,               &
            elf_calc,               &
            kinetic_energy_density
#if defined(HAVE_FFT)
  public :: elf_calc_fs
#endif

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
    call loct_parse_logical(check_inp('ElfWithCurrentTerm'), .true., with_current_term)
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

    FLOAT :: s, f, D0, dens
    integer :: i, is, ik, ist, idim, nelfs
    CMPLX, allocatable :: wf_psi(:), gwf_psi(:,:)
    FLOAT, allocatable :: rho(:, :), grho(:,:), jj(:,:)
    FLOAT, allocatable :: kappa(:, :)

    FLOAT, parameter :: dmin = CNST(1e-10)
#if defined(HAVE_MPI)
    FLOAT, allocatable :: reduce_elf(:)
#endif

    call push_sub('states.states_calc_elf')

    ! We may or may not want the total elf. 
    ! If we want it, the argument elf should have three components.c
    nelfs = size(elf, 2)

    ! single or double occupancy
    if(st%d%nspin == 1) then
      s = M_TWO
    else
      s = M_ONE
    end if

    ALLOCATE(rho(NP, st%d%nspin), NP)
    ALLOCATE(kappa(NP, st%d%nspin), NP)
    rho = M_ZERO
    kappa = M_ZERO
    call states_calc_dens(st, NP, rho)
    rho = rho/s

    do_is: do is = 1, st%d%nspin
      ! gradient of the spin density
      ALLOCATE(grho(NP, NDIM), NP*NDIM)      
      ! spin current
      ALLOCATE(  jj(NP, NDIM), NP*NDIM)      
      ! to store the wave-functions
      ALLOCATE( wf_psi(NP_PART),  NP_PART)   
      ! the gradients of the wave-functions
      ALLOCATE(gwf_psi(NP, NDIM), NP*NDIM)   
      grho = M_ZERO; jj  = M_ZERO; wf_psi = M_ZERO; gwf_psi = M_ZERO

      do ik = is, st%d%nik, st%d%nspin
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim

            ! all calculations will be done with complex wave-functions
            if (st%d%wfs_type == M_REAL) then
              wf_psi(:) = cmplx(st%dpsi(:, idim, ist, ik), KIND=REAL_PRECISION)
            else
              wf_psi(:) = st%zpsi(:, idim, ist, ik)
            end if

            ! calculate gradient of the wave-function
            call zf_gradient(gr%sb, gr%f_der, wf_psi(:), gwf_psi)

            do i = 1, NDIM
              grho(:,i) = grho(:,i) + st%d%kweights(ik)*st%occ(ist, ik)/s *  &
                 M_TWO * real(conjg(wf_psi(:))*gwf_psi(:,i))
              jj(:,i)   =  jj(:,i)  + st%d%kweights(ik)*st%occ(ist, ik)/s *  &
                 aimag(conjg(wf_psi(:))*gwf_psi(:,i))
            end do

            ! elf will now contain the spin-kinetic energy density, tau
            do i = 1, NP
              kappa(i, is) = kappa(i, is) + st%d%kweights(ik)*st%occ(ist, ik)/s * &
                 sum(abs(gwf_psi(i, 1:NDIM))**2)
            end do

          end do
        end do
      end do
      deallocate(wf_psi, gwf_psi)

#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        ALLOCATE(reduce_elf(1:NP), NP)
        call MPI_Allreduce(kappa(1, is), reduce_elf(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
        kappa(1:NP, is) = reduce_elf(1:NP)

        do i = 1, NDIM
          call MPI_Allreduce(grho(1, i), reduce_elf(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
          grho(1:NP, i) = reduce_elf(1:NP)
          
          call MPI_Allreduce(jj(1, i), reduce_elf(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
          jj(1:NP, i) = reduce_elf(1:NP)
        end do
        deallocate(reduce_elf)
      end if
#endif

      if(.not.with_current_term) jj = M_ZERO

      ! kapp will contain rho * D
      do i = 1, NP
        kappa(i, is) = kappa(i, is)*rho(i, is)        &    ! + tau * rho
           - M_FOURTH*sum(grho(i, 1:NDIM)**2) &            ! - | nabla rho |^2 / 4
           - sum(jj(i, 1:NDIM)**2)                         ! - j^2
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
              case(2); D0 = f * rho(i, 1)**3
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


#if defined(HAVE_FFT)
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
 
    if (st%d%wfs_type == M_REAL) then
      call dcf_new(gr%m%l, dcf_tmp)
      call dcf_fft_init(dcf_tmp, gr%sb)
    else
      call zcf_new(gr%m%l, zcf_tmp)
      call zcf_fft_init(zcf_tmp, gr%sb)
    end if

    do_is: do is = 1, st%d%nspin
      ALLOCATE(    r(NP),       NP)
      ALLOCATE(gradr(NP, NDIM), NP*NDIM)
      ALLOCATE(    j(NP, NDIM), NP*NDIM)
      r = M_ZERO; gradr = M_ZERO; j  = M_ZERO

      elf(1:NP,is) = M_ZERO

      ALLOCATE(psi_fs(NP_PART),  NP_PART)
      ALLOCATE(gpsi  (NP, NDIM), NP*NDIM)
      do ik = is, st%d%nik, st%d%nspin
        do ist = 1, st%nst
          do idim = 1, st%d%dim
            
            if (st%d%wfs_type == M_REAL) then
              call dmf2mf_RS2FS(gr%m, st%dpsi(:, idim, ist, ik), psi_fs(:), dcf_tmp)
              call zf_gradient(gr%sb, gr%f_der, psi_fs(:), gpsi)
              do i = 1, NDIM
                gpsi(:,i) = gpsi(:,i) * gr%m%h(i)**2 * real(dcf_tmp%n(i), REAL_PRECISION) / (M_TWO*M_PI)
              end do
            else
              call zmf2mf_RS2FS(gr%m, st%zpsi(:, idim, ist, ik), psi_fs(:), zcf_tmp)
              call zf_gradient(gr%sb, gr%f_der, psi_fs(:), gpsi)
              do i = 1, NDIM
                gpsi(:,i) = gpsi(:,i) * gr%m%h(i)**2 * real(zcf_tmp%n(i), REAL_PRECISION) / (M_TWO*M_PI)
              end do
            end if

            r(:) = r(:) + st%d%kweights(ik)*st%occ(ist, ik) * abs(psi_fs(:))**2
            do i = 1, NDIM
              gradr(:,i) = gradr(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   M_TWO * real(conjg(psi_fs(:))*gpsi(:,i))
              j (:,i) =  j(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   aimag(conjg(psi_fs(:))*gpsi(:,i))
            end do

            do i = 1, NP
              if(r(i) >= dmin) then
                elf(i,is) = elf(i,is) + st%d%kweights(ik)*st%occ(ist, ik)/s * &
                     sum(abs(gpsi(i, 1:NDIM))**2)
              end if
            end do
          end do
        end do
      end do
      deallocate(psi_fs, gpsi)

      do i = 1, NP
        if(r(i) >= dmin) then
          elf(i,is) = elf(i,is) - (M_FOURTH*sum(gradr(i, 1:NDIM)**2) + sum(j(i, 1:NDIM)**2))/(s*r(i))
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

    if (st%d%wfs_type == M_REAL) then
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
      call dmf2cf(m, fin, c)
      call dcf_RS2FS(c)
      call dcf_FS2mf(m, c, fout)
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
      call zmf2cf(m, fin, c)
      call zcf_RS2FS(c)
      call zcf_FS2mf(m, c, fout)
      call zcf_free_RS(c)
      call zcf_free_FS(c)
    end subroutine zmf2mf_RS2FS

  end subroutine elf_calc_fs
#endif

  ! ---------------------------------------------------------
  ! Calculates the kinetic energy density, tau, defined as:
  !
  !  tau(r, is) = sum_{i=1}^{Nst} | \nabla \phi_{i,is}(r) |^2
  !
  ! If the calculation is unpolarized, the subroutine gives
  ! back the sum of the two spin-resolved densities, and the
  ! argument is ignored:
  !
  !  tau(r, 1) = sum_{i=1}^{Nst} 2 | \nabla \phi_{i,is}(r) |^2
  !
  ! The previous expression assumes full or null occupations.
  ! If fractional occupation numbers, each term in the sum
  ! is weigthed by the occupation. Also, if we are working
  ! with an infinite system, all k-points are summed up, with
  ! its corresponding weight.
  !
  ! WARNING: It is not properly programmed for spinors.
  ! ---------------------------------------------------------
  subroutine kinetic_energy_density(st, gr, tau)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    FLOAT,            intent(out)   :: tau(: ,:)

    integer :: i, is, ik, ist, idim
    CMPLX, allocatable :: wf_psi(:), gwf_psi(:,:)
#if defined(HAVE_MPI)
    FLOAT, allocatable :: reduce_elf(:)
#endif

    call push_sub('elf.kinetic_energy_density')

    do is = 1, st%d%nspin

      ALLOCATE( wf_psi(NP_PART),  NP_PART)   
      ALLOCATE(gwf_psi(NP, NDIM), NP*NDIM)   
      wf_psi = M_z0; gwf_psi = M_z0

      tau(:, is) = M_ZERO

      do ik = is, st%d%nik, st%d%nspin
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim

            ! all calculations will be done with complex wave-functions
            if (st%d%wfs_type == M_REAL) then
              wf_psi(:) = cmplx(st%dpsi(:, idim, ist, ik), KIND=REAL_PRECISION)
            else
              wf_psi(:) = st%zpsi(:, idim, ist, ik)
            end if

            ! calculate gradient of the wave-function
            call zf_gradient(gr%sb, gr%f_der, wf_psi(:), gwf_psi)

            ! tau will now contain the spin-kinetic energy density
            do i = 1, NP
              tau(i, is) = tau(i, is) + st%d%kweights(ik)*st%occ(ist, ik) * &
                 sum(abs(gwf_psi(i, 1:NDIM))**2)
            end do

          end do
        end do
      end do
      deallocate(wf_psi, gwf_psi)

#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        ALLOCATE(reduce_elf(1:NP), NP)
        call MPI_Allreduce(tau(1, is), reduce_elf(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
        tau(1:NP, is) = reduce_elf(1:NP)
        deallocate(reduce_elf)
      end if
#endif

    end do


    call pop_sub()
  end subroutine kinetic_energy_density

end module elf_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
