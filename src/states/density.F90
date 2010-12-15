!! Copyright (C) 2002-2010 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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
!! $Id$

#include "global.h"

module density_m
  use batch_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use io_m
  use kpoints_m
  use loct_m
  use mpi_m ! if not before parser_m, ifort 11.072 can`t compile with MPI2
  use math_m
  use mesh_m
  use messages_m
  use multigrid_m
  use multicomm_m
  use mpi_m
  use mpi_lib_m
  use profiling_m
  use simul_box_m
  use smear_m
  use states_m
  use states_dim_m
  use symmetrizer_m
  use types_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m

  implicit none

  private

  public ::                           &
    states_dens_accumulate_batch,     &
    states_dens_reduce,               &
    states_calc_dens,                 &
    states_freeze_orbitals,           &
    states_total_density

  contains

  ! ---------------------------------------------------

  subroutine states_dens_accumulate_batch(st, gr, ik, psib, rho)
    type(states_t), intent(in)    :: st
    type(grid_t),   intent(in)    :: gr
    integer,        intent(in)    :: ik
    type(batch_t),  intent(inout) :: psib
    FLOAT, target,  intent(inout) :: rho(:,:)

    integer :: ist, ist2, ip, ispin
    CMPLX   :: term, psi1, psi2
    FLOAT, pointer :: dpsi(:, :)
    CMPLX, pointer :: zpsi(:, :)
    FLOAT, pointer :: crho(:)
    FLOAT, allocatable :: frho(:)
    type(profile_t), save :: prof

    PUSH_SUB(states_dens_accumulate_batch)
    call profiling_in(prof, "CALC_DENSITY")

    ASSERT(ubound(rho, dim = 1) == gr%fine%mesh%np .or. ubound(rho, dim = 1) == gr%fine%mesh%np_part)

    ispin = states_dim_get_spin_index(st%d, ik)

    if(st%d%ispin /= SPINORS) then 

      if(gr%have_fine_mesh) then
        SAFE_ALLOCATE(crho(1:gr%mesh%np_part))
        crho = M_ZERO
      else
        crho => rho(:, ispin)
      end if

      select case(batch_status(psib))
      case(BATCH_NOT_PACKED, BATCH_CL_PACKED)
        call batch_sync(psib)
        if(states_are_real(st)) then
          do ist = 1, psib%nst
            ist2 = psib%states(ist)%ist
            dpsi => psib%states(ist)%dpsi

            forall(ip = 1:gr%mesh%np)
              crho(ip) = crho(ip) + st%d%kweights(ik) * st%occ(ist2, ik) * dpsi(ip, 1)**2
            end forall
          end do
        else
          do ist = 1, psib%nst
            ist2 = psib%states(ist)%ist
            zpsi => psib%states(ist)%zpsi

            forall(ip = 1:gr%mesh%np)
              crho(ip) = crho(ip) + st%d%kweights(ik) * st%occ(ist2, ik) * &
                (real(zpsi(ip, 1), REAL_PRECISION)**2 + aimag(zpsi(ip, 1))**2)
            end forall
          end do
        end if
      case(BATCH_PACKED)
        if(states_are_real(st)) then
          do ip = 1, gr%mesh%np
            do ist = 1, psib%nst
              ist2 = psib%states(ist)%ist
              crho(ip) = crho(ip) + st%d%kweights(ik)*st%occ(ist2, ik)*psib%pack%dpsi(ist, ip)**2
            end do
          end do
        else
          do ip = 1, gr%mesh%np
            do ist = 1, psib%nst
              ist2 = psib%states(ist)%ist
              crho(ip) = crho(ip) + st%d%kweights(ik)*st%occ(ist2, ik)* &
                (real(psib%pack%zpsi(ist, ip), REAL_PRECISION)**2 + aimag(psib%pack%zpsi(ist, ip))**2)
            end do
          end do
        end if
      end select

      if(gr%have_fine_mesh) then
        SAFE_ALLOCATE(frho(1:gr%fine%mesh%np))
        call dmultigrid_coarse2fine(gr%fine%tt, gr%der, gr%fine%mesh, crho, frho, order = 2)
        ! some debugging output that I will keep here for the moment, XA
        !      call doutput_function(1, "./", "n_fine", gr%fine%mesh, frho, unit_one, ierr)
        !      call doutput_function(1, "./", "n_coarse", gr%mesh, crho, unit_one, ierr)
        forall(ip = 1:gr%fine%mesh%np) rho(ip, ispin) = rho(ip, ispin) + frho(ip)
        SAFE_DEALLOCATE_P(crho)
        SAFE_DEALLOCATE_A(frho)
      end if

    else !SPINORS

      ! in this case wavefunctions are always complex
      ASSERT(.not. gr%have_fine_mesh)
      call batch_sync(psib)

      do ist = 1, psib%nst
        ist2 = psib%states(ist)%ist
        zpsi => psib%states(ist)%zpsi

        do ip = 1, gr%fine%mesh%np

          psi1 = zpsi(ip, 1)
          psi2 = zpsi(ip, 2)

          rho(ip, 1) = rho(ip, 1) + &
            st%d%kweights(ik) * st%occ(ist2, ik) * (real(psi1, REAL_PRECISION)**2 + aimag(psi1)**2)
          rho(ip, 2) = rho(ip, 2) + &
            st%d%kweights(ik) * st%occ(ist2, ik) * (real(psi2, REAL_PRECISION)**2 + aimag(psi2)**2)

          term = st%d%kweights(ik) * st%occ(ist2, ik) * psi1 * conjg(psi2)
          rho(ip, 3) = rho(ip, 3) + real(term, REAL_PRECISION)
          rho(ip, 4) = rho(ip, 4) + aimag(term)

        end do
      end do
      
    end if

    call profiling_out(prof)

    POP_SUB(states_dens_accumulate_batch)
  end subroutine states_dens_accumulate_batch

  ! ---------------------------------------------------

  subroutine states_dens_reduce(st, gr, rho)
    type(states_t), intent(in)    :: st
    type(grid_t),   intent(in)    :: gr
    FLOAT,          intent(inout) :: rho(:,:)

    type(symmetrizer_t) :: symmetrizer
    FLOAT,  allocatable :: symmrho(:)
    integer :: ispin, np
#ifdef HAVE_MPI
    FLOAT,  allocatable :: reduce_rho(:)
    type(profile_t), save :: reduce_prof
#endif

    PUSH_SUB(states_dens_reduce)

    np = gr%fine%mesh%np

#ifdef HAVE_MPI
    ! reduce over states
    if(st%parallel_in_states) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
#ifndef HAVE_MPI2
      SAFE_ALLOCATE(reduce_rho(1:np))
#endif
      do ispin = 1, st%d%nspin
#ifndef HAVE_MPI2
        call blas_copy(np, rho(1, ispin), 1, reduce_rho(1), 1)
#endif
        call MPI_Allreduce(MPI_IN_PLACE_OR(reduce_rho(1)), rho(1, ispin), np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      end do
      SAFE_DEALLOCATE_A(reduce_rho)
      call profiling_out(reduce_prof)
    end if

    ! reduce over k-points
    if(st%d%kpt%parallel) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
      SAFE_ALLOCATE(reduce_rho(1:np))
      do ispin = 1, st%d%nspin
        call MPI_Allreduce(rho(1, ispin), reduce_rho(1), np, MPI_FLOAT, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
        call blas_copy(np, reduce_rho(1), 1, rho(1, ispin), 1)
      end do
      SAFE_DEALLOCATE_A(reduce_rho)
      call profiling_out(reduce_prof)
    end if
#endif

    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symmrho(1:np))
      call symmetrizer_init(symmetrizer, gr%fine%mesh)

      do ispin = 1, st%d%nspin
        call dsymmetrizer_apply(symmetrizer, rho(:, ispin), symmrho)
        rho(1:np, ispin) = symmrho(1:np)
      end do

      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmrho)
    end if

    POP_SUB(states_dens_reduce)
  end subroutine states_dens_reduce


  ! ---------------------------------------------------------
  !> Computes the density from the orbitals in st. If rho is
  !! present, the density is placed there; if it is not present,
  !! the density is placed in st%rho.
  ! ---------------------------------------------------------
  subroutine states_calc_dens(st, gr, rho)
    type(states_t),          intent(inout)  :: st
    type(grid_t),            intent(in)     :: gr
    FLOAT, optional, target, intent(out)    :: rho(:,:)

    integer :: ik
    FLOAT, pointer :: dens(:, :)
    type(batch_t)  :: psib

    PUSH_SUB(states_calc_dens)

    if(present(rho)) then
      dens => rho
    else
      dens => st%rho
    end if

    ASSERT(ubound(dens, dim = 1) == gr%fine%mesh%np .or. ubound(dens, dim = 1) == gr%fine%mesh%np_part)

    dens(1:gr%fine%mesh%np, 1:st%d%nspin) = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end
      if(states_are_real(st)) then
        call batch_init(psib, st%d%dim, st%st_start, st%st_end, st%dpsi(:, :, st%st_start:, ik))
      else
        call batch_init(psib, st%d%dim, st%st_start, st%st_end, st%zpsi(:, :, st%st_start:, ik))
      end if

      call states_dens_accumulate_batch(st, gr, ik, psib, dens)
      
      call batch_end(psib)
    end do

    call states_dens_reduce(st, gr, dens)

    nullify(dens)
    POP_SUB(states_calc_dens)
  end subroutine states_calc_dens

  ! ---------------------------------------------------------

  subroutine states_freeze_orbitals(st, gr, mc, n)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(multicomm_t), intent(in)    :: mc
    integer,           intent(in)    :: n

    integer :: ist, ik
    type(states_t) :: staux
    type(batch_t)  :: psib

    PUSH_SUB(states_freeze_orbitals)

    if(n >= st%nst) then
      write(message(1),'(a)') 'Attempting to freeze a number of orbitals which is larger or equal to'
      write(message(2),'(a)') 'the total number. The program has to stop.'
      call write_fatal(2)
    end if

    ASSERT(.not. st%parallel_in_states)

    if(.not.associated(st%frozen_rho)) then
      SAFE_ALLOCATE(st%frozen_rho(1:gr%mesh%np, 1:st%d%dim))
    end if

    st%frozen_rho = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end
      if(n < st%st_start .or. n > st%st_end) cycle

      if(states_are_real(st)) then
        call batch_init(psib, st%d%dim, st%st_start, n, st%dpsi(:, :, ist:ist, ik))
      else
        call batch_init(psib, st%d%dim, st%st_start, n, st%zpsi(:, :, ist:ist, ik))
      end if
      
      call states_dens_accumulate_batch(st, gr, ik, psib, st%frozen_rho)
      
      call batch_end(psib)
    end do

    call states_dens_reduce(st, gr, st%frozen_rho)

    call states_copy(staux, st)

    st%nst = st%nst - n

    call states_deallocate_wfns(st)
    call states_distribute_nodes(st, mc)
    call states_allocate_wfns(st, gr%mesh, TYPE_CMPLX)

#if defined(HAVE_MPI) 

    if(staux%parallel_in_states) then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = staux%st_start, staux%st_end
          if(ist <= n) cycle
          if(.not.state_is_local(st, ist-n)) then
            call mpi_send(staux%zpsi(1, 1, ist, ik), gr%mesh%np_part*st%d%dim, MPI_CMPLX, staux%node(ist), &
              ist, st%mpi_grp%comm, mpi_err)

            call mpi_recv(st%zpsi(1, 1, ist-n, ik), gr%mesh%np_part*st%d%dim, MPI_CMPLX, st%node(ist-n), &
              ist, st%mpi_grp%comm, mpi_err)
          else
            st%zpsi(:, :, ist-n, ik) = staux%zpsi(:, :, ist, ik)
          end if
   
        end do
      end do
   else
     do ik = st%d%kpt%start, st%d%kpt%end
       do ist = staux%st_start, staux%st_end
         if(ist <= n) cycle
         st%zpsi(:, :, ist-n, ik) = staux%zpsi(:, :, ist, ik)
       end do
     end do
   end if

#else

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%zpsi(:, :, ist, ik) = staux%zpsi(:, :, n + ist, ik)
      end do
    end do

#endif

    SAFE_DEALLOCATE_P(st%occ)
    SAFE_DEALLOCATE_P(st%eigenval)
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)
    st%occ      = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%occ(ist, ik) = staux%occ(n+ist, ik)
        st%eigenval(ist, ik) = staux%eigenval(n+ist, ik)
      end do
    end do

    call states_end(staux)
    POP_SUB(states_freeze_orbitals)
  end subroutine states_freeze_orbitals


  ! ---------------------------------------------------------
  !> this routine calculates the total electronic density,
  !! which is the sum of the part coming from the orbitals, the
  !! non-linear core corrections and the frozen orbitals
  subroutine states_total_density(st, mesh, rho)
    type(states_t), intent(in)  :: st
    type(mesh_t),   intent(in)  :: mesh
    FLOAT,          intent(out) :: rho(:,:)

    integer :: is, ip

    PUSH_SUB(states_total_density)

    forall(is = 1:st%d%nspin, ip = 1:mesh%np)
      rho(ip, is) = st%rho(ip, is)
    end forall

    if(associated(st%rho_core)) then
      forall(is = 1:st%d%spin_channels, ip = 1:mesh%np)
        rho(ip, is) = rho(ip, is) + st%rho_core(ip)/st%d%nspin
      end forall
    end if

    ! Add, if it exists, the frozen density from the inner orbitals.
    if(associated(st%frozen_rho)) then
      forall(is = 1:st%d%spin_channels, ip = 1:mesh%np)
        rho(ip, is) = rho(ip, is) + st%frozen_rho(ip, is)
      end forall
    end if

    POP_SUB(states_total_density)
  end subroutine states_total_density

end module density_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
