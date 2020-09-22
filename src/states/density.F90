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
  use profiling_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use states_oct_m
  use states_dim_oct_m
  use symmetrizer_oct_m
  use types_oct_m

  implicit none

  private

  public ::                           &
    density_calc_t,                   &
    density_calc_init,                &
    density_calc_accumulate,          &
    density_calc_end,                 &
    density_calc,                     &
    states_freeze_orbitals,           &
    states_total_density,             &
    ddensity_accumulate_grad,         &
    zdensity_accumulate_grad

  type density_calc_t
    FLOAT,                pointer :: density(:, :)
    type(states_t),       pointer :: st
    type(grid_t),         pointer :: gr
    type(accel_mem_t)            :: buff_density
    integer                       :: pnp
    logical                       :: packed
  end type density_calc_t

contains
  
  subroutine density_calc_init(this, st, gr, density)
    type(density_calc_t),           intent(out)   :: this
    type(states_t),       target,   intent(in)    :: st
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

  subroutine density_calc_accumulate(this, ik, psib, dim)
    type(density_calc_t),         intent(inout) :: this
    integer,                      intent(in)    :: ik
    type(batch_t),                intent(inout) :: psib
    integer, optional,            intent(in)    :: dim

    integer :: ist, ip, ispin, dim_
    FLOAT   :: nrm
    CMPLX   :: term, psi1, psi2
    CMPLX, allocatable :: psi(:), fpsi(:), zpsi(:, :)
    FLOAT, allocatable :: weight(:), sqpsi(:)
    type(profile_t), save :: prof
    integer            :: wgsize
    type(accel_mem_t) :: buff_weight
    type(accel_kernel_t), pointer :: kernel

    PUSH_SUB(density_calc_accumulate)
    call profiling_in(prof, "CALC_DENSITY")

    dim_ = optional_default(dim, 1)

    if(dim_ /= 1) then
       ASSERT(.not. this%gr%have_fine_mesh)
       ASSERT(.not.batch_status(psib)== BATCH_PACKED)
    end if
    
    ispin = states_dim_get_spin_index(this%st%d, ik)

    SAFE_ALLOCATE(weight(1:psib%nst))
    forall(ist = 1:psib%nst) weight(ist) = this%st%d%kweights(ik)*this%st%occ(psib%states(ist)%ist, ik)

    if(this%st%d%ispin /= SPINORS .and. .not. this%gr%have_fine_mesh) then 

      select case(batch_status(psib))
      case(BATCH_NOT_PACKED)
        if(states_are_real(this%st)) then
          do ist = 1, psib%nst
            forall(ip = 1:this%gr%mesh%np)
              this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*psib%states(ist)%dpsi(ip,dim_)**2
            end forall
          end do
        else
          do ist = 1, psib%nst
            forall(ip = 1:this%gr%mesh%np)
              this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)* &
                (real(psib%states(ist)%zpsi(ip, dim_), REAL_PRECISION)**2 + aimag(psib%states(ist)%zpsi(ip, dim_))**2)
            end forall
          end do
        end if
      case(BATCH_PACKED)
        if(states_are_real(this%st)) then
          do ip = 1, this%gr%mesh%np
            do ist = 1, psib%nst
              this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*psib%pack%dpsi(ist, ip)**2
            end do
          end do
        else
          do ip = 1, this%gr%mesh%np
            do ist = 1, psib%nst
              this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)* &
                (real(psib%pack%zpsi(ist, ip), REAL_PRECISION)**2 + aimag(psib%pack%zpsi(ist, ip))**2)
            end do
          end do
        end if
      case(BATCH_DEVICE_PACKED)
        if(.not. this%packed) call density_calc_pack(this)

        if(states_are_real(this%st)) then
          kernel => kernel_density_real
        else
          kernel => kernel_density_complex
        end if

        call accel_create_buffer(buff_weight, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, psib%nst)
        call accel_write_buffer(buff_weight, psib%nst, weight)

        call accel_set_kernel_arg(kernel, 0, psib%nst)
        call accel_set_kernel_arg(kernel, 1, this%gr%mesh%np)
        call accel_set_kernel_arg(kernel, 2, this%pnp*(ispin - 1))
        call accel_set_kernel_arg(kernel, 3, buff_weight)
        call accel_set_kernel_arg(kernel, 4, psib%pack%buffer)
        call accel_set_kernel_arg(kernel, 5, log2(psib%pack%size(1)))
        call accel_set_kernel_arg(kernel, 6, this%buff_density)

        wgsize = accel_kernel_workgroup_size(kernel)
        
        call accel_kernel_run(kernel, (/pad(this%gr%mesh%np, wgsize)/), (/wgsize/))

        call accel_finish()
        
        call accel_release_buffer(buff_weight)
        
      end select

    else if(this%gr%have_fine_mesh) then

      SAFE_ALLOCATE(psi(1:this%gr%mesh%np_part))
      SAFE_ALLOCATE(fpsi(1:this%gr%fine%mesh%np))
      SAFE_ALLOCATE(sqpsi(1:this%gr%fine%mesh%np))

      do ist = 1, psib%nst

        call batch_get_state(psib, ist, this%gr%mesh%np, psi)

        call zmultigrid_coarse2fine(this%gr%fine%tt, this%gr%der, this%gr%fine%mesh, psi, fpsi, order = 2)

        do ip = 1, this%gr%fine%mesh%np
          sqpsi(ip) = abs(fpsi(ip))**2
        end do

        nrm = dmf_integrate(this%gr%fine%mesh, sqpsi)
        
        do ip = 1, this%gr%fine%mesh%np
          this%density(ip, ispin) = this%density(ip, ispin) + weight(ist)*sqpsi(ip)/nrm
        end do
        
      end do


      ! some debugging output that I will keep here for the moment, XA
      !      call dio_function_output(1, "./", "n_fine", this%gr%fine%mesh, frho, unit_one, ierr)
      !      call dio_function_output(1, "./", "n_coarse", this%gr%mesh, crho, unit_one, ierr)
      !        forall(ip = 1:this%gr%fine%mesh%np) this%density(ip, ispin) = this%density(ip, ispin) + frho(ip)

      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(fpsi)
      SAFE_DEALLOCATE_A(sqpsi)
      
    else !SPINORS

      ! in this case wavefunctions are always complex
      ASSERT(.not. this%gr%have_fine_mesh)

      SAFE_ALLOCATE(zpsi(1:this%gr%mesh%np, 1:this%st%d%dim))

      do ist = 1, psib%nst
        call batch_get_state(psib, ist, this%gr%mesh%np, zpsi)
        
        do ip = 1, this%gr%fine%mesh%np

          psi1 = zpsi(ip,(dim_-1)*2+1)
          psi2 = zpsi(ip,(dim_-1)*2+2)

          this%density(ip, 1) = this%density(ip, 1) + weight(ist)*(real(psi1, REAL_PRECISION)**2 + aimag(psi1)**2)
          this%density(ip, 2) = this%density(ip, 2) + weight(ist)*(real(psi2, REAL_PRECISION)**2 + aimag(psi2)**2)

          term = weight(ist)*psi1*conjg(psi2)
          this%density(ip, 3) = this%density(ip, 3) + real(term, REAL_PRECISION)
          this%density(ip, 4) = this%density(ip, 4) + aimag(term)

        end do
      end do

      SAFE_DEALLOCATE_A(zpsi)
      
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
  subroutine density_calc(st, gr, density)
    type(states_t),          intent(inout)  :: st
    type(grid_t),            intent(in)     :: gr
    FLOAT,                   intent(out)    :: density(:, :)

    integer :: ik, ib,ndim, idim, dim_stride
    type(density_calc_t) :: dens_calc

    PUSH_SUB(density_calc)

    ASSERT(ubound(density, dim = 1) == gr%fine%mesh%np .or. ubound(density, dim = 1) == gr%fine%mesh%np_part)

    ndim = 1
    dim_stride = 1
    
    ! In order to handle the cases where st%d%dim encodes more than 
    ! just spinors we contract the additional dimensions
    select case(st%d%ispin)
    case(UNPOLARIZED)
      if(st%d%dim>1) then
        ndim = st%d%dim
      end if

    case(SPIN_POLARIZED)
      if(st%d%dim>1) then
        ndim = st%d%dim
      end if

    case(SPINORS)
      if(st%d%dim>2) then
        dim_stride = 2
        ndim = st%d%dim
      end if

    end select
    
    
    call density_calc_init(dens_calc, st, gr, density)
    
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
         do idim=1,ndim, dim_stride
             call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik),dim=idim)
         end do
      end do
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

    integer :: ist, ik, ib, nblock
    type(states_t) :: staux
    CMPLX, allocatable :: psi(:, :, :)
    type(batch_t)  :: psib
    type(density_calc_t) :: dens_calc
#ifdef HAVE_MPI
    integer :: status(MPI_STATUS_SIZE)
#endif

    PUSH_SUB(states_freeze_orbitals)

    if(n >= st%nst) then
      write(message(1),'(a)') 'Attempting to freeze a number of orbitals which is larger or equal to'
      write(message(2),'(a)') 'the total number. The program has to stop.'
      call messages_fatal(2)
    end if

    ASSERT(.not. st%parallel_in_states)

    if(.not.associated(st%frozen_rho)) then
      SAFE_ALLOCATE(st%frozen_rho(1:gr%mesh%np, 1:st%d%nspin))
    end if

    call density_calc_init(dens_calc, st, gr, st%frozen_rho)

    do ik = st%d%kpt%start, st%d%kpt%end
      if(n < st%st_start .or. n > st%st_end) cycle

      do ib =  st%group%block_start, st%group%block_end
        if(states_block_max(st, ib) <= n) then

          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))
          if(states_block_max(st, ib) == n) exit

        else 

          nblock = n - states_block_min(st, ib) + 1

          SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim, 1:nblock))

          do ist = 1, nblock
            call states_get_state(st, gr%mesh, states_block_min(st, ib) + ist - 1, ik, psi(:, :, ist)) 
          end do

          call batch_init(psib, st%d%dim, states_block_min(st, ib), n, psi)
          call density_calc_accumulate(dens_calc, ik, psib)
          call batch_end(psib)
          SAFE_DEALLOCATE_A(psi)
          
          exit

        end if

      end do
    end do

    call density_calc_end(dens_calc)

    call states_copy(staux, st)

    st%nst = st%nst - n

    call states_deallocate_wfns(st)
    call states_distribute_nodes(st, mc)
    call states_allocate_wfns(st, gr%mesh, TYPE_CMPLX)

    SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim, 1))
    
#if defined(HAVE_MPI) 

    if(staux%parallel_in_states) then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = staux%st_start, staux%st_end
          if(ist <= n) cycle
          if(.not.state_is_local(st, ist-n)) then

            call states_get_state(staux, gr%mesh, ist, ik, psi(:, :, 1))

            ! I think this can cause a deadlock. XA
            call MPI_Send(psi(1, 1, 1), gr%mesh%np_part*st%d%dim, MPI_CMPLX, staux%node(ist), ist, &
              st%mpi_grp%comm, mpi_err)

            call MPI_Recv(psi(1, 1, 1), gr%mesh%np_part*st%d%dim, MPI_CMPLX, st%node(ist - n), &
              ist, st%mpi_grp%comm, status, mpi_err)

            call states_set_state(st, gr%mesh, ist - n, ik, psi(:, :, 1))
            
          else
            call states_get_state(staux, gr%mesh, ist, ik, psi(:, :, 1))
            call states_set_state(st, gr%mesh, ist - n, ik, psi(:, :, 1))
          end if
   
        end do
      end do
   else
     do ik = st%d%kpt%start, st%d%kpt%end
       do ist = staux%st_start, staux%st_end
         if(ist <= n) cycle
         call states_get_state(staux, gr%mesh, ist, ik, psi(:, :, 1))
         call states_set_state(st, gr%mesh, ist - n, ik, psi(:, :, 1))
       end do
     end do
   end if

#else

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        call states_get_state(staux, gr%mesh, ist + n, ik, psi(:, :, 1))
        call states_set_state(st, gr%mesh, ist, ik, psi(:, :, 1))
      end do
    end do

#endif

    SAFE_DEALLOCATE_A(psi)
    
    ! Change the smearing method by fixing the occupations to 
    ! that of the ground-state such that the unfrozen states inherit 
    ! those values.
    st%smear%method = SMEAR_FIXED_OCC
  
    ! Set total charge
    st%qtot = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%qtot = st%qtot + staux%occ(n+ist, ik) * st%d%kweights(ik)
      end do
    end do

    SAFE_DEALLOCATE_P(st%eigenval)
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)

    SAFE_DEALLOCATE_P(st%occ)
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
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
  subroutine states_total_density(st, mesh, total_rho)
    type(states_t),  intent(in)  :: st
    type(mesh_t),    intent(in)  :: mesh
    FLOAT,           intent(out) :: total_rho(:,:)

    integer :: is, ip

    PUSH_SUB(states_total_density)

    forall(ip = 1:mesh%np, is = 1:st%d%nspin)
      total_rho(ip, is) = st%rho(ip, is)
    end forall

    if(associated(st%rho_core)) then
      forall(ip = 1:mesh%np, is = 1:st%d%spin_channels)
        total_rho(ip, is) = total_rho(ip, is) + st%rho_core(ip)/st%d%spin_channels
      end forall
    end if

    ! Add, if it exists, the frozen density from the inner orbitals.
    if(associated(st%frozen_rho)) then
      forall(ip = 1:mesh%np, is = 1:st%d%spin_channels)
        total_rho(ip, is) = total_rho(ip, is) + st%frozen_rho(ip, is)
      end forall
    end if
  
    POP_SUB(states_total_density)
  end subroutine states_total_density

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
