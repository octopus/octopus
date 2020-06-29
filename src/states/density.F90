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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module density_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use iso_c_binding
  use comm_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use kpoints_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multigrid_oct_m
  use multicomm_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use symmetrizer_oct_m
  use types_oct_m
  use wfs_elec_oct_m

  implicit none

  private

  public ::                           &
    density_calc_t,                   &
    density_calc_init,                &
    density_calc_accumulate,          &
    density_calc_end,                 &
    density_calc,                     &
    states_elec_freeze_orbitals,           &
    states_elec_total_density,             &
    ddensity_accumulate_grad,         &
    zdensity_accumulate_grad,         &
    states_elec_freeze_redistribute_states,&
    states_elec_freeze_adjust_qtot

  type density_calc_t
    private
    FLOAT,                pointer :: density(:, :)
    type(states_elec_t),  pointer :: st
    type(grid_t),         pointer :: gr
    type(accel_mem_t)            :: buff_density
    integer                       :: pnp
    logical                       :: packed
  end type density_calc_t

contains
  
  subroutine density_calc_init(this, st, gr, density)
    type(density_calc_t),           intent(out)   :: this
    type(states_elec_t),  target,   intent(in)    :: st
    type(grid_t),         target,   intent(in)    :: gr
    FLOAT,                target,   intent(out)   :: density(:, :)

    logical :: correct_size

    PUSH_SUB(density_calc_init)

    this%st => st
    this%gr => gr

    this%density => density
    this%density = M_ZERO

    this%packed = .false.

    correct_size = ubound(this%density, dim = 1) == this%gr%fine%mesh%np .or. &
         ubound(this%density, dim = 1) == this%gr%fine%mesh%np_part
    ASSERT(correct_size)

    POP_SUB(density_calc_init)
  end subroutine density_calc_init

  ! ---------------------------------------------------

  subroutine density_calc_pack(this)
    type(density_calc_t),           intent(inout)   :: this

    PUSH_SUB(density_calc_pack)
    
    this%packed = .true.
    this%pnp = accel_padded_size(this%gr%mesh%np)
    call accel_create_buffer(this%buff_density, ACCEL_MEM_READ_WRITE, TYPE_FLOAT, this%pnp*this%st%d%nspin)
    
    ! set to zero
    call accel_set_buffer_to_zero(this%buff_density, TYPE_FLOAT, this%pnp*this%st%d%nspin)
    
    POP_SUB(density_calc_pack)
  end subroutine density_calc_pack

  ! ---------------------------------------------------

  subroutine density_calc_accumulate(this, psib, istin)
    type(density_calc_t),         intent(inout) :: this
    type(wfs_elec_t),             intent(in)    :: psib
    integer, optional,            intent(in)    :: istin

    integer :: ist, ip, ispin, istin_
    FLOAT   :: nrm
    CMPLX   :: term, psi1, psi2
    CMPLX, allocatable :: psi(:), fpsi(:)
    FLOAT, allocatable :: weight(:), sqpsi(:)
    type(profile_t), save :: prof
    integer            :: wgsize
    type(accel_mem_t) :: buff_weight
    type(accel_kernel_t), pointer :: kernel
    logical  :: select_state

    PUSH_SUB(density_calc_accumulate)
    call profiling_in(prof, "CALC_DENSITY")
    
    istin_ = optional_default(istin, -1)
    select_state = (istin_ > 0)

    ispin = states_elec_dim_get_spin_index(this%st%d, psib%ik)

    SAFE_ALLOCATE(weight(1:psib%nst))
    do ist = 1, psib%nst
      weight(ist) = this%st%d%kweights(psib%ik)*this%st%occ(psib%ist(ist), psib%ik)
    end do

    if (.not. this%gr%have_fine_mesh) then 

      select case(psib%status())
      case(BATCH_NOT_PACKED)
        select case (this%st%d%ispin)
        case (UNPOLARIZED, SPIN_POLARIZED)
          if(states_are_real(this%st)) then
            do ist = 1, psib%nst
              if (select_state .and. psib%ist(ist) /= istin_) cycle
              if(abs(weight(ist)) <= M_EPSILON) cycle
              !$omp parallel do simd schedule(static)
              do ip = 1, this%gr%mesh%np
                this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*psib%dff(ip, 1, ist)**2
              end do
            end do
          else
            do ist = 1, psib%nst
              if (select_state .and. psib%ist(ist) /= istin_) cycle
              if(abs(weight(ist)) <= M_EPSILON) cycle
              !$omp parallel do schedule(static)
              do ip = 1, this%gr%mesh%np
                this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)* &
                  TOFLOAT(conjg(psib%zff(ip, 1, ist))*psib%zff(ip, 1, ist))
              end do
            end do
          end if
        case (SPINORS)
          do ist = 1, psib%nst
            if (select_state .and. psib%ist(ist) /= istin_) cycle
            if(abs(weight(ist)) <= M_EPSILON) cycle
            !$omp parallel do schedule(static) private(psi1, psi2, term)
            do ip = 1, this%gr%mesh%np          
              psi1 = psib%zff(ip, 1, ist)
              psi2 = psib%zff(ip, 2, ist)
              this%density(ip, 1) = this%density(ip, 1) + weight(ist)*TOFLOAT(conjg(psi1)*psi1)
              this%density(ip, 2) = this%density(ip, 2) + weight(ist)*TOFLOAT(conjg(psi2)*psi2)
              term = weight(ist)*psi1*conjg(psi2)
              this%density(ip, 3) = this%density(ip, 3) + TOFLOAT(term)
              this%density(ip, 4) = this%density(ip, 4) + aimag(term)
            end do
          end do
        end select

      case(BATCH_PACKED)

        select case (this%st%d%ispin)
        case (UNPOLARIZED, SPIN_POLARIZED)
          if(states_are_real(this%st)) then
            !$omp parallel do schedule(static)
            do ip = 1, this%gr%mesh%np
              do ist = 1, psib%nst
                if (select_state .and. psib%ist(ist) /= istin_) cycle
                this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*psib%dff_pack(ist, ip)**2
              end do
            end do
          else
            !$omp parallel do schedule(static)
            do ip = 1, this%gr%mesh%np
              do ist = 1, psib%nst
                if (select_state .and. psib%ist(ist) /= istin_) cycle
                this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)* &
                  TOFLOAT(conjg(psib%zff_pack(ist, ip))*psib%zff_pack(ist, ip))
              end do
            end do
          end if
        case (SPINORS)
          ASSERT(mod(psib%nst_linear, 2) == 0)
          !$omp parallel do schedule(static) private(ist, psi1, psi2, term)
          do ip = 1, this%gr%mesh%np
            do ist = 1, psib%nst
              if (select_state .and. psib%ist(ist) /= istin_) cycle
              psi1 = psib%zff_pack(2*ist - 1, ip)
              psi2 = psib%zff_pack(2*ist,     ip)
              term = weight(ist)*psi1*conjg(psi2)

              this%density(ip, 1) = this%density(ip, 1) + weight(ist)*TOFLOAT(conjg(psi1)*psi1)
              this%density(ip, 2) = this%density(ip, 2) + weight(ist)*TOFLOAT(conjg(psi2)*psi2)
              this%density(ip, 3) = this%density(ip, 3) + TOFLOAT(term)
              this%density(ip, 4) = this%density(ip, 4) + aimag(term)
            end do
          end do
        end select

      case(BATCH_DEVICE_PACKED)
        if(.not. this%packed) call density_calc_pack(this)

        call accel_create_buffer(buff_weight, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, psib%nst)
        call accel_write_buffer(buff_weight, psib%nst, weight)

        select case (this%st%d%ispin)
        case (UNPOLARIZED, SPIN_POLARIZED)        
          if(states_are_real(this%st)) then
            kernel => kernel_density_real
          else
            kernel => kernel_density_complex
          end if

          call accel_set_kernel_arg(kernel, 0, psib%nst)
          call accel_set_kernel_arg(kernel, 1, this%gr%mesh%np)
          call accel_set_kernel_arg(kernel, 2, this%pnp*(ispin - 1))
          call accel_set_kernel_arg(kernel, 3, buff_weight)
          call accel_set_kernel_arg(kernel, 4, psib%ff_device)
          call accel_set_kernel_arg(kernel, 5, log2(psib%pack_size(1)))
          call accel_set_kernel_arg(kernel, 6, this%buff_density)
        
        case (SPINORS)
          kernel => kernel_density_spinors

          call accel_set_kernel_arg(kernel, 0, psib%nst)
          call accel_set_kernel_arg(kernel, 1, this%gr%mesh%np)
          call accel_set_kernel_arg(kernel, 2, this%pnp)
          call accel_set_kernel_arg(kernel, 3, buff_weight)
          call accel_set_kernel_arg(kernel, 4, psib%ff_device)
          call accel_set_kernel_arg(kernel, 5, log2(psib%pack_size(1)))
          call accel_set_kernel_arg(kernel, 6, this%buff_density)
        end select

        wgsize = accel_kernel_workgroup_size(kernel)

        call accel_kernel_run(kernel, (/pad(this%gr%mesh%np, wgsize)/), (/wgsize/))

        call accel_finish()
        
        call accel_release_buffer(buff_weight)
        
      end select

    else if(this%gr%have_fine_mesh) then
      ! Non-collinear density not implemented with fine grid
      ASSERT(this%st%d%ispin /= SPINORS)

      SAFE_ALLOCATE(psi(1:this%gr%mesh%np_part))
      SAFE_ALLOCATE(fpsi(1:this%gr%fine%mesh%np))
      SAFE_ALLOCATE(sqpsi(1:this%gr%fine%mesh%np))

      do ist = 1, psib%nst

        if (select_state .and. psib%ist(ist) /= istin_) cycle
        if(abs(weight(ist)) <= M_EPSILON) cycle

        call batch_get_state(psib, ist, this%gr%mesh%np, psi)

        call zmultigrid_coarse2fine(this%gr%fine%tt, this%gr%der, this%gr%fine%mesh, psi, fpsi, order = 2)

        !$omp parallel do schedule(static)
        do ip = 1, this%gr%fine%mesh%np
          sqpsi(ip) = TOFLOAT(conjg(fpsi(ip))*fpsi(ip))
        end do

        nrm = dmf_integrate(this%gr%fine%mesh, sqpsi)
        
        !$omp parallel do schedule(static)
        do ip = 1, this%gr%fine%mesh%np
          this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*sqpsi(ip)/nrm
        end do
        
      end do


      ! For this to work again, a namespace has to be available. I don`t make
      ! the change now because this debugging output has been commented out for
      ! 4 years. - SO
      ! some debugging output that I will keep here for the moment, XA
      !      call dio_function_output(1, "./", "n_fine", this%gr%fine%mesh, frho, unit_one, ierr)
      !      call dio_function_output(1, "./", "n_coarse", this%gr%mesh, crho, unit_one, ierr)
      !        do ip = 1, this%gr%fine%mesh%np
      !            this%density(ip, ispin) = this%density(ip, ispin) + frho(ip)
      !        end do

      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(fpsi)
      SAFE_DEALLOCATE_A(sqpsi)

    end if

    SAFE_DEALLOCATE_A(weight)

    call profiling_out(prof)

    POP_SUB(density_calc_accumulate)
  end subroutine density_calc_accumulate

  ! ---------------------------------------------------

  subroutine density_calc_end(this)
    type(density_calc_t), intent(inout) :: this

    type(symmetrizer_t) :: symmetrizer
    FLOAT, allocatable :: tmpdensity(:)
    integer :: ispin, ip
    type(profile_t), save :: reduce_prof
    FLOAT, allocatable :: fdensity(:)

    PUSH_SUB(density_calc_end)

    if(this%packed) then
      SAFE_ALLOCATE(tmpdensity(1:this%gr%mesh%np_part))

      ! the density is in device memory
      do ispin = 1, this%st%d%nspin
        call accel_read_buffer(this%buff_density, this%gr%mesh%np, tmpdensity, offset = (ispin - 1)*this%pnp)

        if(this%gr%have_fine_mesh) then
           SAFE_ALLOCATE(fdensity(1:this%gr%fine%mesh%np))
           call dmultigrid_coarse2fine(this%gr%fine%tt, this%gr%der, this%gr%fine%mesh, tmpdensity, fdensity, order = 2)

           do ip = 1, this%gr%fine%mesh%np
             this%density(ip, ispin) = this%density(ip, ispin) + fdensity(ip)
           end do

           SAFE_DEALLOCATE_A(fdensity)
        else
          do ip = 1, this%gr%mesh%np
            this%density(ip, ispin) = this%density(ip, ispin) + tmpdensity(ip)
          end do
        end if

      end do

      this%packed = .false.
      call accel_release_buffer(this%buff_density)
      SAFE_DEALLOCATE_A(tmpdensity)
    end if

    ! reduce over states and k-points
    if(this%st%parallel_in_states .or. this%st%d%kpt%parallel) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
      call comm_allreduce(this%st%st_kpt_mpi_grp%comm, this%density, dim = (/this%gr%fine%mesh%np, this%st%d%nspin/))
      call profiling_out(reduce_prof)
    end if

    if(this%st%symmetrize_density) then
      SAFE_ALLOCATE(tmpdensity(1:this%gr%fine%mesh%np))
      call symmetrizer_init(symmetrizer, this%gr%fine%mesh)

      do ispin = 1, this%st%d%nspin
        call dsymmetrizer_apply(symmetrizer, this%gr%fine%mesh%np, field = this%density(:, ispin), &
                                 symmfield = tmpdensity)
        this%density(1:this%gr%fine%mesh%np, ispin) = tmpdensity(1:this%gr%fine%mesh%np)
      end do

      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(tmpdensity)
    end if

    POP_SUB(density_calc_end)
  end subroutine density_calc_end


  ! ---------------------------------------------------------
  !> Computes the density from the orbitals in st. 
  subroutine density_calc(st, gr, density, istin)
    type(states_elec_t),     intent(in)  :: st
    type(grid_t),            intent(in)  :: gr
    FLOAT,                   intent(out) :: density(:, :)
    integer, optional,       intent(in)  :: istin 

    integer :: ik, ib
    type(density_calc_t) :: dens_calc

    PUSH_SUB(density_calc)

    ASSERT(ubound(density, dim = 1) == gr%fine%mesh%np .or. ubound(density, dim = 1) == gr%fine%mesh%np_part)

    call density_calc_init(dens_calc, st, gr, density)
    
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        if (present(istin)) then
          call density_calc_accumulate(dens_calc, st%group%psib(ib, ik), istin = istin)
        else
          call density_calc_accumulate(dens_calc, st%group%psib(ib, ik))
        end if
      end do
    end do

    call density_calc_end(dens_calc)

    POP_SUB(density_calc)
  end subroutine density_calc

  ! ---------------------------------------------------------

  subroutine states_elec_freeze_orbitals(st, namespace, gr, mc, n, family_is_mgga)
    type(states_elec_t), intent(inout) :: st
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(in)    :: gr
    type(multicomm_t),   intent(in)    :: mc
    integer,             intent(in)    :: n
    logical,             intent(in)    :: family_is_mgga

    integer :: ist, ik, ib, nblock, st_min
    integer :: nodeto, nodefr, nsend, nreceiv
    type(states_elec_t) :: staux
    CMPLX, allocatable :: psi(:, :, :), rec_buffer(:,:)
    type(wfs_elec_t)  :: psib
    type(density_calc_t) :: dens_calc
#ifdef HAVE_MPI
    integer :: status(MPI_STATUS_SIZE)
#endif

    PUSH_SUB(states_elec_freeze_orbitals)

    if(n >= st%nst) then
      write(message(1),'(a)') 'Attempting to freeze a number of orbitals which is larger or equal to'
      write(message(2),'(a)') 'the total number. The program has to stop.'
      call messages_fatal(2, namespace=namespace)
    end if

    ASSERT(states_are_complex(st))

    if(.not.associated(st%frozen_rho)) then
      SAFE_ALLOCATE(st%frozen_rho(1:gr%mesh%np, 1:st%d%nspin))
    end if

    call density_calc_init(dens_calc, st, gr, st%frozen_rho)

    do ik = st%d%kpt%start, st%d%kpt%end
      if(n < st%st_start) cycle

      do ib =  st%group%block_start, st%group%block_end
        !We can use the full batch 
        if(states_elec_block_max(st, ib) <= n) then

          call density_calc_accumulate(dens_calc, st%group%psib(ib, ik))
          if(states_elec_block_max(st, ib) == n) exit

        else !Here we only use a part of this batch 

          nblock = n - states_elec_block_min(st, ib) + 1

          SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim, 1:nblock))

          do ist = 1, nblock
            call states_elec_get_state(st, gr%mesh, states_elec_block_min(st, ib) + ist - 1, ik, psi(:, :, ist)) 
          end do

          call wfs_elec_init(psib, st%d%dim, states_elec_block_min(st, ib), n, psi, ik)
          call density_calc_accumulate(dens_calc, psib)
          call psib%end()
          SAFE_DEALLOCATE_A(psi)
          
          exit

        end if

      end do
    end do

    call density_calc_end(dens_calc)

    if(family_is_mgga) then
      if(.not.associated(st%frozen_tau)) then
        SAFE_ALLOCATE(st%frozen_tau(1:gr%mesh%np, 1:st%d%nspin))
      end if    
      if(.not.associated(st%frozen_gdens)) then
        SAFE_ALLOCATE(st%frozen_gdens(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%nspin))
      end if
      if(.not.associated(st%frozen_ldens)) then
        SAFE_ALLOCATE(st%frozen_ldens(1:gr%mesh%np, 1:st%d%nspin))
      end if

      call states_elec_calc_quantities(gr%der, st, .true., kinetic_energy_density = st%frozen_tau, &
           density_gradient = st%frozen_gdens, density_laplacian = st%frozen_ldens, st_end = n) 
    end if 


    call states_elec_copy(staux, st)

    call states_elec_freeze_redistribute_states(st, namespace, gr, mc, n)

    SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim, 1))
    SAFE_ALLOCATE(rec_buffer(1:gr%mesh%np, 1:st%d%dim))

    if(staux%parallel_in_states) then

      do ik = st%d%kpt%start, st%d%kpt%end
        !We cound how many states we have to send, and how many we  will receive
        nsend = 0
        do ist = staux%st_start, staux%st_end
          if(ist > n) nsend = nsend + 1
        end do
        nreceiv = st%st_end-st%st_start+1

        st_min = min(max(staux%st_start-n, 1),st%st_start)
 
        do ist = st_min, st%st_end
          nodeto = -1
          nodefr = -1
          if(nsend > 0 .and. ist+n <= staux%st_end) nodeto = st%node(ist)
          if(nreceiv > 0 .and. ist >= st%st_start) nodefr = staux%node(ist+n)

          !Local copy
          if(nsend >0 .and. nreceiv>0 .and. nodeto == nodefr .and. nodefr == st%mpi_grp%rank) then
            call states_elec_get_state(staux, gr%mesh, ist+n, ik, psi(:, :, 1))
            call states_elec_set_state(st, gr%mesh, ist, ik, psi(:, :, 1))            
            nsend = nsend -1
            nreceiv= nreceiv-1
          else
            if(nsend > 0 .and. nodeto > -1 .and. nodeto /= st%mpi_grp%rank) then
              call states_elec_get_state(staux, gr%mesh, ist+n, ik, psi(:, :, 1))
#if defined(HAVE_MPI)
              call MPI_Send(psi(1, 1, 1), gr%mesh%np*st%d%dim, MPI_CMPLX, nodeto, ist, &
                    st%mpi_grp%comm, mpi_err)
#endif
              nsend = nsend -1
            end if

            if(nreceiv > 0 .and. nodefr > -1 .and. nodefr /= st%mpi_grp%rank) then
#if defined(HAVE_MPI)
              call MPI_Recv(rec_buffer(1, 1), gr%mesh%np*st%d%dim, MPI_CMPLX, nodefr, &
                 ist, st%mpi_grp%comm, status, mpi_err)
#endif
              call states_elec_set_state(st, gr%mesh, ist, ik, rec_buffer(:, :))
              nreceiv= nreceiv-1
            end if
          end if
        end do
      end do

      ! Add a barrier to ensure that the process are synchronized
#if defined(HAVE_MPI)
      call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
   
    else
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call states_elec_get_state(staux, gr%mesh, ist + n, ik, psi(:, :, 1))
          call states_elec_set_state(st, gr%mesh, ist, ik, psi(:, :, 1))
        end do
      end do
    end if

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(rec_buffer)
    
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        st%occ(ist, ik) = staux%occ(n+ist, ik)
        st%eigenval(ist, ik) = staux%eigenval(n+ist, ik)
      end do
    end do

    call states_elec_end(staux)
    POP_SUB(states_elec_freeze_orbitals)
  end subroutine states_elec_freeze_orbitals

  ! ---------------------------------------------------------
  subroutine states_elec_freeze_redistribute_states(st, namespace, gr, mc, nn)
    type(states_elec_t), intent(inout) :: st
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(in)    :: gr
    type(multicomm_t),   intent(in)    :: mc
    integer,             intent(in)    :: nn

    PUSH_SUB(states_elec_freeze_redistribute_states)

    st%nst = st%nst - nn

    call states_elec_deallocate_wfns(st)
    call states_elec_distribute_nodes(st, namespace, mc)
    call states_elec_allocate_wfns(st, gr%mesh, TYPE_CMPLX)

    SAFE_DEALLOCATE_P(st%eigenval)
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)

    SAFE_DEALLOCATE_P(st%occ)
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
    st%occ      = M_ZERO


    POP_SUB(states_elec_freeze_redistribute_states)
  end subroutine states_elec_freeze_redistribute_states

  ! ---------------------------------------------------------
  subroutine states_elec_freeze_adjust_qtot(st)
    type(states_elec_t), intent(inout) :: st

    integer :: ik, ist

    PUSH_SUB(states_elec_freeze_adjust_occs)

    ! Change the smearing method by fixing the occupations to  
    ! that of the ground-state such that the unfrozen states inherit 
    ! those values.
    st%smear%method = SMEAR_FIXED_OCC 

    ! Set total charge
    st%qtot = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%qtot = st%qtot + st%occ(ist, ik) * st%d%kweights(ik)
      end do
    end do

#if defined(HAVE_MPI)        
    if(st%parallel_in_states .or. st%d%kpt%parallel) then
      call comm_allreduce(st%st_kpt_mpi_grp%comm, st%qtot)
    end if
#endif  


    POP_SUB(states_elec_freeze_adjust_qtot)
  end subroutine states_elec_freeze_adjust_qtot


  ! ---------------------------------------------------------
  !> this routine calculates the total electronic density,
  !! which is the sum of the part coming from the orbitals, the
  !! non-linear core corrections and the frozen orbitals
  subroutine states_elec_total_density(st, mesh, total_rho)
    type(states_elec_t),  intent(in)  :: st
    type(mesh_t),         intent(in)  :: mesh
    FLOAT,                intent(out) :: total_rho(:,:)

    integer :: is, ip

    PUSH_SUB(states_elec_total_density)

    do is = 1, st%d%nspin
      do ip = 1, mesh%np
        total_rho(ip, is) = st%rho(ip, is)
      end do
    end do

    if(associated(st%rho_core)) then
      do is = 1, st%d%spin_channels
        do ip = 1, mesh%np
          total_rho(ip, is) = total_rho(ip, is) + st%rho_core(ip)/st%d%spin_channels
        end do
      end do
    end if

    ! Add, if it exists, the frozen density from the inner orbitals.
    if(associated(st%frozen_rho)) then
      do is = 1, st%d%nspin
        do ip = 1, mesh%np
          total_rho(ip, is) = total_rho(ip, is) + st%frozen_rho(ip, is)
        end do
      end do
    end if

    POP_SUB(states_elec_total_density)
  end subroutine states_elec_total_density

#include "undef.F90"
#include "real.F90"
#include "density_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "density_inc.F90"

end module density_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
