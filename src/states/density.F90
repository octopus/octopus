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
  use blas_m
  use batch_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use io_m
  use kpoints_m
  use loct_m
  use math_m
  use mesh_m
  use messages_m
  use multigrid_m
  use multicomm_m
  use mpi_m ! if not before parser_m, ifort 11.072 can`t compile with MPI2
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
    density_calc_t,                   &
    density_calc_init,               &
    density_calc_accumulate,          &
    density_calc_end,                 &
    density_calc,                     &
    states_freeze_orbitals,           &
    states_total_density

  type density_calc_t
    FLOAT,          pointer :: density(:, :)
    type(states_t), pointer :: st
    type(grid_t),   pointer :: gr
  end type density_calc_t

contains
  
  subroutine density_calc_init(this, st, gr, density)
    type(density_calc_t),           intent(out)   :: this
    type(states_t),       target,   intent(in)    :: st
    type(grid_t),         target,   intent(in)    :: gr
    FLOAT,                target,   intent(out)   :: density(:, :)

    PUSH_SUB(density_calc_init)

    this%density => density
    this%st => st
    this%gr => gr

    this%density = M_ZERO

    POP_SUB(density_calc_init)
  end subroutine density_calc_init

  ! ---------------------------------------------------

  subroutine density_calc_accumulate(this, ik, psib)
    type(density_calc_t), intent(inout) :: this
    integer,              intent(in)    :: ik
    type(batch_t),        intent(inout) :: psib

    integer :: ist, ist2, ip, ispin
    CMPLX   :: term, psi1, psi2
    FLOAT, pointer :: dpsi(:, :)
    CMPLX, pointer :: zpsi(:, :)
    FLOAT, pointer :: crho(:)
    FLOAT, allocatable :: frho(:)
    type(profile_t), save :: prof
    logical :: correct_size

    PUSH_SUB(density_calc_accumulate)
    call profiling_in(prof, "CALC_DENSITY")

    correct_size = ubound(this%density, dim = 1) == this%gr%fine%mesh%np &
      .or. ubound(this%density, dim = 1) == this%gr%fine%mesh%np_part

    ispin = states_dim_get_spin_index(this%st%d, ik)

    if(this%st%d%ispin /= SPINORS) then 

      if(this%gr%have_fine_mesh) then
        SAFE_ALLOCATE(crho(1:this%gr%mesh%np_part))
        crho = M_ZERO
      else
        crho => this%density(:, ispin)
      end if

      select case(batch_status(psib))
      case(BATCH_NOT_PACKED, BATCH_CL_PACKED)
        call batch_sync(psib)
        if(states_are_real(this%st)) then
          do ist = 1, psib%nst
            ist2 = psib%states(ist)%ist
            dpsi => psib%states(ist)%dpsi

            forall(ip = 1:this%gr%mesh%np)
              crho(ip) = crho(ip) + this%st%d%kweights(ik) * this%st%occ(ist2, ik) * dpsi(ip, 1)**2
            end forall
          end do
        else
          do ist = 1, psib%nst
            ist2 = psib%states(ist)%ist
            zpsi => psib%states(ist)%zpsi

            forall(ip = 1:this%gr%mesh%np)
              crho(ip) = crho(ip) + this%st%d%kweights(ik) * this%st%occ(ist2, ik) * &
                (real(zpsi(ip, 1), REAL_PRECISION)**2 + aimag(zpsi(ip, 1))**2)
            end forall
          end do
        end if
      case(BATCH_PACKED)
        if(states_are_real(this%st)) then
          do ip = 1, this%gr%mesh%np
            do ist = 1, psib%nst
              ist2 = psib%states(ist)%ist
              crho(ip) = crho(ip) + this%st%d%kweights(ik)*this%st%occ(ist2, ik)*psib%pack%dpsi(ist, ip)**2
            end do
          end do
        else
          do ip = 1, this%gr%mesh%np
            do ist = 1, psib%nst
              ist2 = psib%states(ist)%ist
              crho(ip) = crho(ip) + this%st%d%kweights(ik)*this%st%occ(ist2, ik)* &
                (real(psib%pack%zpsi(ist, ip), REAL_PRECISION)**2 + aimag(psib%pack%zpsi(ist, ip))**2)
            end do
          end do
        end if
      end select

      if(this%gr%have_fine_mesh) then
        SAFE_ALLOCATE(frho(1:this%gr%fine%mesh%np))
        call dmultigrid_coarse2fine(this%gr%fine%tt, this%gr%der, this%gr%fine%mesh, crho, frho, order = 2)
        ! some debugging output that I will keep here for the moment, XA
        !      call doutput_function(1, "./", "n_fine", this%gr%fine%mesh, frho, unit_one, ierr)
        !      call doutput_function(1, "./", "n_coarse", this%gr%mesh, crho, unit_one, ierr)
        forall(ip = 1:this%gr%fine%mesh%np) this%density(ip, ispin) = this%density(ip, ispin) + frho(ip)
        SAFE_DEALLOCATE_P(crho)
        SAFE_DEALLOCATE_A(frho)
      end if

    else !SPINORS

      ! in this case wavefunctions are always complex
      ASSERT(.not. this%gr%have_fine_mesh)
      call batch_sync(psib)

      do ist = 1, psib%nst
        ist2 = psib%states(ist)%ist
        zpsi => psib%states(ist)%zpsi

        do ip = 1, this%gr%fine%mesh%np

          psi1 = zpsi(ip, 1)
          psi2 = zpsi(ip, 2)

          this%density(ip, 1) = this%density(ip, 1) + &
            this%st%d%kweights(ik) * this%st%occ(ist2, ik) * (real(psi1, REAL_PRECISION)**2 + aimag(psi1)**2)
          this%density(ip, 2) = this%density(ip, 2) + &
            this%st%d%kweights(ik) * this%st%occ(ist2, ik) * (real(psi2, REAL_PRECISION)**2 + aimag(psi2)**2)

          term = this%st%d%kweights(ik) * this%st%occ(ist2, ik) * psi1 * conjg(psi2)
          this%density(ip, 3) = this%density(ip, 3) + real(term, REAL_PRECISION)
          this%density(ip, 4) = this%density(ip, 4) + aimag(term)

        end do
      end do
      
    end if

    call profiling_out(prof)

    POP_SUB(density_calc_accumulate)
  end subroutine density_calc_accumulate

  ! ---------------------------------------------------

  subroutine density_calc_end(this)
    type(density_calc_t), intent(inout) :: this

    type(symmetrizer_t) :: symmetrizer
    FLOAT,  allocatable :: symmrho(:)
    integer :: ispin, np
#ifdef HAVE_MPI
    FLOAT,  allocatable :: reduce_rho(:)
    type(profile_t), save :: reduce_prof
#endif

    PUSH_SUB(density_calc_end)

    np = this%gr%fine%mesh%np

#ifdef HAVE_MPI
    ! reduce over states
    if(this%st%parallel_in_states) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
#ifndef HAVE_MPI2
      SAFE_ALLOCATE(reduce_rho(1:np))
#endif
      do ispin = 1, this%st%d%nspin
#ifndef HAVE_MPI2
        call blas_copy(np, this%density(1, ispin), 1, reduce_rho(1), 1)
#endif
        call MPI_Allreduce(MPI_IN_PLACE_OR(reduce_rho(1)), &
          this%density(1, ispin), np, MPI_FLOAT, MPI_SUM, this%st%mpi_grp%comm, mpi_err)
      end do
      SAFE_DEALLOCATE_A(reduce_rho)
      call profiling_out(reduce_prof)
    end if

    ! reduce over k-points
    if(this%st%d%kpt%parallel) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
      SAFE_ALLOCATE(reduce_rho(1:np))
      do ispin = 1, this%st%d%nspin
        call MPI_Allreduce(this%density(1, ispin), reduce_rho(1), np, MPI_FLOAT, MPI_SUM, this%st%d%kpt%mpi_grp%comm, mpi_err)
        call blas_copy(np, reduce_rho(1), 1, this%density(1, ispin), 1)
      end do
      SAFE_DEALLOCATE_A(reduce_rho)
      call profiling_out(reduce_prof)
    end if
#endif

    if(this%st%symmetrize_density) then
      SAFE_ALLOCATE(symmrho(1:np))
      call symmetrizer_init(symmetrizer, this%gr%fine%mesh)

      do ispin = 1, this%st%d%nspin
        call dsymmetrizer_apply(symmetrizer, this%density(:, ispin), symmrho)
        this%density(1:np, ispin) = symmrho(1:np)
      end do

      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmrho)
    end if

    POP_SUB(density_calc_end)
  end subroutine density_calc_end


  ! ---------------------------------------------------------
  !> Computes the density from the orbitals in st. 
  ! ---------------------------------------------------------
  subroutine density_calc(st, gr, density)
    type(states_t),          intent(inout)  :: st
    type(grid_t),            intent(in)     :: gr
    FLOAT,                   intent(out)    :: density(:, :)

    integer :: ik
    type(batch_t)        :: psib
    type(density_calc_t) :: dens_calc

    PUSH_SUB(density_calc)

    ASSERT(ubound(density, dim = 1) == gr%fine%mesh%np .or. ubound(density, dim = 1) == gr%fine%mesh%np_part)

    call density_calc_init(dens_calc, st, gr, density)

    do ik = st%d%kpt%start, st%d%kpt%end
      if(states_are_real(st)) then
        call batch_init(psib, st%d%dim, st%st_start, st%st_end, st%dpsi(:, :, st%st_start:, ik))
      else
        call batch_init(psib, st%d%dim, st%st_start, st%st_end, st%zpsi(:, :, st%st_start:, ik))
      end if

      call density_calc_accumulate(dens_calc, ik, psib)
      
      call batch_end(psib)
    end do

    call density_calc_end(dens_calc)

    POP_SUB(density_calc)
  end subroutine density_calc

  ! ---------------------------------------------------------

  subroutine states_freeze_orbitals(st, gr, mc, n)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(multicomm_t), intent(in)    :: mc
    integer,           intent(in)    :: n

    integer :: ist, ik
    type(states_t) :: staux
    type(batch_t)  :: psib
    type(density_calc_t) :: dens_calc

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

    call density_calc_init(dens_calc, st, gr, st%frozen_rho)

    do ik = st%d%kpt%start, st%d%kpt%end
      if(n < st%st_start .or. n > st%st_end) cycle

      if(states_are_real(st)) then
        call batch_init(psib, st%d%dim, st%st_start, n, st%dpsi(:, :, ist:ist, ik))
      else
        call batch_init(psib, st%d%dim, st%st_start, n, st%zpsi(:, :, ist:ist, ik))
      end if
      
      call density_calc_accumulate(dens_calc, ik, psib)
      
      call batch_end(psib)
    end do

    call density_calc_end(dens_calc)

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
