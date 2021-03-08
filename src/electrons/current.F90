!! Copyright (C) 2008-2019 X. Andrade, F. Bonafe, R. Jestaedt, H. Appel
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

module current_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use exchange_operator_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_elec_base_oct_m
  use hamiltonian_elec_oct_m
  use lalg_basic_oct_m
  use lda_u_oct_m
  use math_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use projector_oct_m
  use scissor_oct_m
  use simul_box_oct_m
  use states_elec_dim_oct_m
  use states_elec_oct_m
  use string_oct_m
  use symmetries_oct_m
  use symmetrizer_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m  
  use varinfo_oct_m
  use wfs_elec_oct_m
  use xc_oct_m
  
  implicit none

  private

  type current_t
    private
    integer :: method
  end type current_t
    

  public ::                               &
    current_t,                            &
    current_init,                         &
    current_calculate,                    &
    current_heat_calculate,               &
    current_calculate_mel

  integer, parameter, public ::           &
    CURRENT_GRADIENT           = 1,       &
    CURRENT_GRADIENT_CORR      = 2,       &
    CURRENT_HAMILTONIAN        = 3

contains

  subroutine current_init(this, namespace)
    type(current_t),   intent(out)   :: this
    type(namespace_t), intent(in)    :: namespace

    PUSH_SUB(current_init)

    !%Variable CurrentDensity
    !%Default gradient_corrected
    !%Type integer
    !%Section Hamiltonian
    !%Description
    !% This variable selects the method used to
    !% calculate the current density. For the moment this variable is
    !% for development purposes and users should not need to use
    !% it.
    !%Option gradient 1
    !% The calculation of current is done using the gradient operator. (Experimental)
    !%Option gradient_corrected 2
    !% The calculation of current is done using the gradient operator
    !% with additional corrections for the total current from non-local operators.
    !%Option hamiltonian 3
    !% The current density is obtained from the commutator of the
    !% Hamiltonian with the position operator. (Experimental)
    !%End

    call parse_variable(namespace, 'CurrentDensity', CURRENT_GRADIENT_CORR, this%method)
    if(.not.varinfo_valid_option('CurrentDensity', this%method)) call messages_input_error(namespace, 'CurrentDensity')
    if(this%method /= CURRENT_GRADIENT_CORR) then
      call messages_experimental("CurrentDensity /= gradient_corrected")
    end if
    
    POP_SUB(current_init)
  end subroutine current_init

  ! ---------------------------------------------------------

  subroutine current_batch_accumulate(st, der, ik, ib, psib, gpsib)
    type(states_elec_t), intent(inout) :: st
    type(derivatives_t), intent(inout) :: der
    integer,             intent(in)    :: ik
    integer,             intent(in)    :: ib
    type(wfs_elec_t),    intent(in)    :: psib
    type(wfs_elec_t),    intent(in)    :: gpsib(:)

    integer :: ist, idir, ii, ip, idim, wgsize
    CMPLX, allocatable :: psi(:, :), gpsi(:, :)
    FLOAT, allocatable :: current_tmp(:, :)
    CMPLX :: c_tmp
    FLOAT :: ww
    FLOAT, allocatable :: weight(:)
    type(accel_mem_t) :: buff_weight, buff_current
    type(accel_kernel_t), save :: kernel
        
    SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:der%mesh%np_part, 1:st%d%dim))

    if(st%d%ispin == SPINORS .or. (psib%status() == BATCH_DEVICE_PACKED .and. der%dim /= 3)) then

      do idir = 1, der%dim
        do ist = states_elec_block_min(st, ib), states_elec_block_max(st, ib)

          ww = st%d%kweights(ik)*st%occ(ist, ik)
          if(abs(ww) <= M_EPSILON) cycle

          do idim = 1, st%d%dim
            ii = st%group%psib(ib, ik)%inv_index((/ist, idim/))
            call batch_get_state(psib, ii, der%mesh%np, psi(:, idim))
            call batch_get_state(gpsib(idir), ii, der%mesh%np, gpsi(:, idim))
          end do

          if(st%d%ispin /= SPINORS) then
            !$omp parallel do
            do ip = 1, der%mesh%np
              st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) + ww*aimag(conjg(psi(ip, 1))*gpsi(ip, 1))
            end do
            !$omp end parallel do
          else
            !$omp parallel do private(c_tmp)
            do ip = 1, der%mesh%np
              st%current(ip, idir, 1) = st%current(ip, idir, 1) + ww*aimag(conjg(psi(ip, 1))*gpsi(ip, 1))
              st%current(ip, idir, 2) = st%current(ip, idir, 2) + ww*aimag(conjg(psi(ip, 2))*gpsi(ip, 2))
              c_tmp = conjg(psi(ip, 1))*gpsi(ip, 2) - psi(ip, 2)*conjg(gpsi(ip, 1))
              st%current(ip, idir, 3) = st%current(ip, idir, 3) + ww*TOFLOAT(c_tmp)
              st%current(ip, idir, 4) = st%current(ip, idir, 4) + ww*aimag(c_tmp)
            end do
            !$omp end parallel do
          end if

        end do
      end do

    else if(psib%status() == BATCH_DEVICE_PACKED) then

      ASSERT(der%dim == 3)

      SAFE_ALLOCATE(weight(1:psib%nst))
      do ist = 1, psib%nst
        weight(ist) = st%d%kweights(ik)*st%occ(psib%ist(ist), ik)
      end do


      call accel_create_buffer(buff_weight, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, psib%nst)
      call accel_write_buffer(buff_weight, psib%nst, weight)

      call accel_create_buffer(buff_current, ACCEL_MEM_WRITE_ONLY, TYPE_FLOAT, der%mesh%np*3)
     
      call accel_kernel_start_call(kernel, 'density.cl', 'current_accumulate')
      
      call accel_set_kernel_arg(kernel, 0, psib%nst)
      call accel_set_kernel_arg(kernel, 1, der%mesh%np)
      call accel_set_kernel_arg(kernel, 2, buff_weight)
      call accel_set_kernel_arg(kernel, 3, psib%ff_device)
      call accel_set_kernel_arg(kernel, 4, log2(psib%pack_size(1)))
      call accel_set_kernel_arg(kernel, 5, gpsib(1)%ff_device)
      call accel_set_kernel_arg(kernel, 6, gpsib(2)%ff_device)
      call accel_set_kernel_arg(kernel, 7, gpsib(3)%ff_device)
      call accel_set_kernel_arg(kernel, 8, log2(gpsib(1)%pack_size(1)))
      call accel_set_kernel_arg(kernel, 9, buff_current)
      
      wgsize = accel_kernel_workgroup_size(kernel)
      
      call accel_kernel_run(kernel, (/pad(der%mesh%np, wgsize)/), (/wgsize/))
      
      SAFE_ALLOCATE(current_tmp(1:der%dim, der%mesh%np))

      call accel_finish()

      call accel_read_buffer(buff_current, der%mesh%np*3, current_tmp)

      do ip = 1, der%mesh%np
        do idir = 1, der%dim
          st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) + current_tmp(idir, ip)
        end do
      end do
      
      SAFE_DEALLOCATE_A(current_tmp)
      
      call accel_release_buffer(buff_weight)
      call accel_release_buffer(buff_current)

      SAFE_DEALLOCATE_A(weight)
      
    else

      ASSERT(psib%is_packed() .eqv. gpsib(1)%is_packed())

      do ii = 1, psib%nst
        ist = states_elec_block_min(st, ib) + ii - 1
        ww = st%d%kweights(ik)*st%occ(ist, ik)
        if(abs(ww) <= M_EPSILON) cycle

        if(psib%is_packed()) then
          do idir = 1, der%dim
            !$omp parallel do
            do ip = 1, der%mesh%np
              st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) &
                + ww*aimag(conjg(psib%zff_pack(ii, ip))*gpsib(idir)%zff_pack(ii, ip))
            end do
            !$omp end parallel do
          end do
        else
          do idir = 1, der%dim
            !$omp parallel do
            do ip = 1, der%mesh%np
              st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) &
                + ww*aimag(conjg(psib%zff(ip, 1, ii))*gpsib(idir)%zff(ip, 1, ii))
            end do
            !$omp end parallel do
          end do          
        end if
        
      end do

    end if

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(gpsi)

  end subroutine current_batch_accumulate

  ! ---------------------------------------------------------
  subroutine current_calculate(this, namespace, der, hm, geo, st, symm)
    type(current_t),          intent(in)    :: this
    type(namespace_t),        intent(in)    :: namespace
    type(derivatives_t),      intent(inout) :: der
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(geometry_t),         intent(in)    :: geo
    type(states_elec_t),      intent(inout) :: st
    type(symmetries_t),       intent(in)    :: symm

    integer :: ik, ist, idir, idim, ip, ib, ii, ispin
    CMPLX, allocatable :: gpsi(:, :, :), psi(:, :), hpsi(:, :), rhpsi(:, :), rpsi(:, :), hrpsi(:, :)
    FLOAT, allocatable :: symmcurrent(:, :)
    type(profile_t), save :: prof
    type(symmetrizer_t) :: symmetrizer
    type(wfs_elec_t) :: hpsib, rhpsib, rpsib, hrpsib, epsib
    class(wfs_elec_t), allocatable :: commpsib(:)
    FLOAT :: ww
    CMPLX :: c_tmp

    call profiling_in(prof, "CURRENT")
    PUSH_SUB(current_calculate)

    ! spin not implemented or tested
    ASSERT(all(ubound(st%current) == (/der%mesh%np_part, der%dim, st%d%nspin/)))
    ASSERT(all(ubound(st%current_kpt) == (/der%mesh%np, der%dim, st%d%kpt%end/)))
    ASSERT(all(lbound(st%current_kpt) == (/1, 1, st%d%kpt%start/)))

    SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:der%mesh%np, 1:der%dim, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rhpsi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rpsi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hrpsi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE_TYPE_ARRAY(wfs_elec_t, commpsib, (1:der%dim))

    st%current = M_ZERO
    st%current_kpt = M_ZERO

    select case(this%method)

    case(CURRENT_HAMILTONIAN)

      do ik = st%d%kpt%start, st%d%kpt%end
        ispin = st%d%get_spin_index(ik)
        do ib = st%group%block_start, st%group%block_end

          call st%group%psib(ib, ik)%do_pack(copy = .true.)

          call st%group%psib(ib, ik)%copy_to(hpsib)
          call st%group%psib(ib, ik)%copy_to(rhpsib)
          call st%group%psib(ib, ik)%copy_to(rpsib)
          call st%group%psib(ib, ik)%copy_to(hrpsib)

          call boundaries_set(der%boundaries, st%group%psib(ib, ik))
          call zhamiltonian_elec_apply_batch(hm, namespace, der%mesh, st%group%psib(ib, ik), hpsib, set_bc = .false.)

          do idir = 1, der%dim

            call batch_mul(der%mesh%np, der%mesh%x(:, idir), hpsib, rhpsib)
            call batch_mul(der%mesh%np_part, der%mesh%x(:, idir), st%group%psib(ib, ik), rpsib)

            call zhamiltonian_elec_apply_batch(hm, namespace, der%mesh, rpsib, hrpsib, set_bc = .false.)

            do ist = states_elec_block_min(st, ib), states_elec_block_max(st, ib)
              ww = st%d%kweights(ik)*st%occ(ist, ik)
              if(ww <= M_EPSILON) cycle

              do idim = 1, st%d%dim
                ii = st%group%psib(ib, ik)%inv_index((/ist, idim/))
                call batch_get_state(st%group%psib(ib, ik), ii, der%mesh%np, psi(:, idim))
                call batch_get_state(hrpsib, ii, der%mesh%np, hrpsi(:, idim))
                call batch_get_state(rhpsib, ii, der%mesh%np, rhpsi(:, idim))
              end do

              if(st%d%ispin /= SPINORS) then
                !$omp parallel do
                do ip = 1, der%mesh%np
                  st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) &
                    - ww*aimag(conjg(psi(ip, 1))*hrpsi(ip, 1) - conjg(psi(ip, 1))*rhpsi(ip, 1))
                end do
                !$omp end parallel do
              else
                !$omp parallel do  private(c_tmp)
                do ip = 1, der%mesh%np
                  st%current(ip, idir, 1) = st%current(ip, idir, 1) + &
                    ww*aimag(conjg(psi(ip, 1))*hrpsi(ip, 1) - conjg(psi(ip, 1))*rhpsi(ip, 1))
                  st%current(ip, idir, 2) = st%current(ip, idir, 2) + &
                    ww*aimag(conjg(psi(ip, 2))*hrpsi(ip, 2) - conjg(psi(ip, 2))*rhpsi(ip, 2))
                  c_tmp = conjg(psi(ip, 1))*hrpsi(ip, 2) - conjg(psi(ip, 1))*rhpsi(ip, 2) &
                    -psi(ip, 2)*conjg(hrpsi(ip, 1)) - psi(ip, 2)*conjg(rhpsi(ip, 1))
                  st%current(ip, idir, 3) = st%current(ip, idir, 3) + ww*TOFLOAT(c_tmp)
                  st%current(ip, idir, 4) = st%current(ip, idir, 4) + ww*aimag(c_tmp)
                end do
                !$omp end parallel do
              end if

            end do

          end do

          call st%group%psib(ib, ik)%do_unpack(copy = .false.)

          call hpsib%end()
          call rhpsib%end()
          call rpsib%end()
          call hrpsib%end()

        end do
      end do

    case(CURRENT_GRADIENT, CURRENT_GRADIENT_CORR)

      if(this%method == CURRENT_GRADIENT_CORR .and. .not. family_is_mgga_with_exc(hm%xc) &
        .and. hm%lda_u_level == DFT_U_NONE .and. hm%theory_level /= HARTREE_FOCK &
        .and. hm%theory_level /= RDMFT) then

        ! we can use the packed version
        
        do ik = st%d%kpt%start, st%d%kpt%end
          ispin = st%d%get_spin_index(ik)
          do ib = st%group%block_start, st%group%block_end

            call st%group%psib(ib, ik)%do_pack(copy = .true.)
            call st%group%psib(ib, ik)%copy_to(epsib)
            call boundaries_set(der%boundaries, st%group%psib(ib, ik))

            if (allocated(hm%hm_base%phase)) then
              call hamiltonian_elec_base_phase(hm%hm_base, der%mesh, der%mesh%np_part, &
                conjugate = .false., psib = epsib, src = st%group%psib(ib, ik))
            else
              call st%group%psib(ib, ik)%copy_data_to(der%mesh%np_part, epsib)
            end if

            ! this now takes non-orthogonal axis into account
            do idir = 1, der%dim
              call epsib%copy_to(commpsib(idir))
            end do
            call zderivatives_batch_grad(der, epsib, commpsib, set_bc=.false.)

            call zhamiltonian_elec_base_nlocal_position_commutator(hm%hm_base, der%mesh, st%d, &
                    der%boundaries, epsib, commpsib)


            call current_batch_accumulate(st, der, ik, ib, epsib, commpsib)

            do idir = 1, der%dim
              call commpsib(idir)%end()
            end do

            call epsib%end()
            call st%group%psib(ib, ik)%do_unpack(copy = .false.)

          end do
        end do

      else

        ! use the slow non-packed version
        
        do ik = st%d%kpt%start, st%d%kpt%end
          ispin = st%d%get_spin_index(ik)
          do ist = st%st_start, st%st_end

            ww = st%d%kweights(ik)*st%occ(ist, ik)
            if(abs(ww) <= M_EPSILON) cycle

            call states_elec_get_state(st, der%mesh, ist, ik, psi)

            do idim = 1, st%d%dim
              call boundaries_set(der%boundaries, psi(:, idim))
            end do

            if (allocated(hm%hm_base%phase)) then 
              call states_elec_set_phase(st%d, psi, hm%hm_base%phase(1:der%mesh%np_part, ik), der%mesh%np_part, .false.)
            end if

            do idim = 1, st%d%dim
              call zderivatives_grad(der, psi(:, idim), gpsi(:, :, idim), set_bc = .false.)
            end do

            if(this%method == CURRENT_GRADIENT_CORR) then
              !A nonlocal contribution from the MGGA potential must be included
              !This must be done first, as this is like a position-dependent mass 
              if (family_is_mgga_with_exc(hm%xc)) then
                do idim = 1, st%d%dim
                  do idir = 1, der%dim
                    !$omp parallel do
                    do ip = 1, der%mesh%np
                      gpsi(ip, idir, idim) = (M_ONE+CNST(2.0)*hm%vtau(ip,ispin))*gpsi(ip, idir, idim)
                    end do
                    !$omp end parallel do
                  end do
                end do
              end if

              !A nonlocal contribution from the pseudopotential must be included
              call zprojector_commute_r_allatoms_alldir(hm%ep%proj, geo, der%mesh, st%d%dim, &
                          der%boundaries, ik, psi, gpsi)                 
              !A nonlocal contribution from the scissor must be included
              if(hm%scissor%apply) then
                call scissor_commute_r(hm%scissor, der%mesh, ik, psi, gpsi)
              end if

              if(hm%lda_u_level /= DFT_U_NONE) then
                call zlda_u_commute_r(hm%lda_u, der%mesh, st%d, namespace, ik, psi, gpsi, allocated(hm%hm_base%phase))
              end if

              call zexchange_operator_commute_r(hm%exxop, der%mesh, st%d, ik, psi, gpsi)

            end if

            if(st%d%ispin /= SPINORS) then
              do idir = 1, der%dim
                !$omp parallel do
                do ip = 1, der%mesh%np
                  st%current_kpt(ip, idir, ik) = st%current_kpt(ip, idir, ik) + &
                    ww*aimag(conjg(psi(ip, 1))*gpsi(ip, idir, 1))
                end do
                !$omp end parallel do
              end do
            else
              do idir = 1, der%dim
                !$omp parallel do  private(c_tmp)
                do ip = 1, der%mesh%np
                  st%current(ip, idir, 1) = st%current(ip, idir, 1) + &
                    ww*aimag(conjg(psi(ip, 1))*gpsi(ip, idir, 1))
                  st%current(ip, idir, 2) = st%current(ip, idir, 2) + &
                    ww*aimag(conjg(psi(ip, 2))*gpsi(ip, idir, 2))
                  c_tmp = conjg(psi(ip, 1))*gpsi(ip, idir, 2) - psi(ip, 2)*conjg(gpsi(ip, idir, 1))
                  st%current(ip, idir, 3) = st%current(ip, idir, 3) + ww*TOFLOAT(c_tmp)
                  st%current(ip, idir, 4) = st%current(ip, idir, 4) + ww*aimag(c_tmp)
                end do
                !$omp end parallel do
              end do
            end if

          end do
        end do

      end if

    case default

      ASSERT(.false.)

    end select

    if(st%d%ispin /= SPINORS) then
      !We sum the current over k-points
      do ik = st%d%kpt%start, st%d%kpt%end
        ispin = st%d%get_spin_index(ik)
        call lalg_axpy(der%mesh%np, der%dim, M_ONE, st%current_kpt(:, :, ik), st%current(:, :, ispin))
      end do
    end if

    if(st%parallel_in_states .or. st%d%kpt%parallel) then
      call comm_allreduce(st%st_kpt_mpi_grp, st%current, dim = (/der%mesh%np, der%dim, st%d%nspin/)) 
    end if

    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symmcurrent(1:der%mesh%np, 1:der%dim))
      call symmetrizer_init(symmetrizer, der%mesh, symm)
      do ispin = 1, st%d%nspin
        call dsymmetrizer_apply(symmetrizer, der%mesh, field_vector = st%current(:, :, ispin), &
          symmfield_vector = symmcurrent, suppress_warning = .true.)
        st%current(1:der%mesh%np, 1:der%dim, ispin) = symmcurrent(1:der%mesh%np, 1:der%dim)
      end do
      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmcurrent)
    end if

    SAFE_DEALLOCATE_A(gpsi)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(rhpsi)
    SAFE_DEALLOCATE_A(rpsi)
    SAFE_DEALLOCATE_A(hrpsi)
    SAFE_DEALLOCATE_A(commpsib)

    call profiling_out(prof)
    POP_SUB(current_calculate)

  end subroutine current_calculate

  
  ! ---------------------------------------------------------
  ! Calculate the current matrix element between two states
  ! I_{ij}(t) = <i| J(t) |j>
  ! This is used only in the floquet_observables utility and 
  ! is highly experimental
  
  subroutine current_calculate_mel(der, hm, geo, psi_i, psi_j, ik,  cmel)
    type(derivatives_t),  intent(inout) :: der
    type(hamiltonian_elec_t),  intent(in)    :: hm
    type(geometry_t),     intent(in)    :: geo
    CMPLX,                intent(in)    :: psi_i(:,:)
    CMPLX,                intent(in)    :: psi_j(:,:)
    integer,              intent(in)    :: ik
    CMPLX,                intent(out)   :: cmel(:,:) ! the current vector cmel(1:der%dim, 1:st%d%nspin)

    integer ::  idir, idim, ip, ispin
    CMPLX, allocatable :: gpsi_j(:, :, :), ppsi_j(:,:),  gpsi_i(:, :, :), ppsi_i(:,:)

    PUSH_SUB(current_calculate_mel)

    SAFE_ALLOCATE(gpsi_i(1:der%mesh%np, 1:der%dim, 1:hm%d%dim))
    SAFE_ALLOCATE(ppsi_i(1:der%mesh%np_part,1:hm%d%dim))
    SAFE_ALLOCATE(gpsi_j(1:der%mesh%np, 1:der%dim, 1:hm%d%dim))
    SAFE_ALLOCATE(ppsi_j(1:der%mesh%np_part,1:hm%d%dim))

    cmel = M_z0

    ispin = hm%d%get_spin_index(ik)
    ppsi_i(:,:) = M_z0        
    ppsi_i(1:der%mesh%np,:) = psi_i(1:der%mesh%np,:)    
    ppsi_j(:,:) = M_z0        
    ppsi_j(1:der%mesh%np,:) = psi_j(1:der%mesh%np,:)    

      
    do idim = 1, hm%d%dim
      call boundaries_set(der%boundaries, ppsi_i(:, idim))
      call boundaries_set(der%boundaries, ppsi_j(:, idim))
    end do

    if (allocated(hm%hm_base%phase)) then 
      ! Apply the phase that contains both the k-point and vector-potential terms.
      do idim = 1, hm%d%dim
        !$omp parallel do
        do ip = 1, der%mesh%np_part
          ppsi_i(ip, idim) = hm%hm_base%phase(ip, ik)*ppsi_i(ip, idim)
          ppsi_j(ip, idim) = hm%hm_base%phase(ip, ik)*ppsi_j(ip, idim)
        end do
        !$omp end parallel do
      end do
    end if

    do idim = 1, hm%d%dim
      call zderivatives_grad(der, ppsi_i(:, idim), gpsi_i(:, :, idim), set_bc = .false.)
      call zderivatives_grad(der, ppsi_j(:, idim), gpsi_j(:, :, idim), set_bc = .false.)
    end do
    
    !A nonlocal contribution from the MGGA potential must be included
    !This must be done first, as this is like a position-dependent mass 
    if (family_is_mgga_with_exc(hm%xc)) then
      do idim = 1, hm%d%dim
        do idir = 1, der%dim
          !$omp parallel do
          do ip = 1, der%mesh%np
            gpsi_i(ip, idir, idim) = (M_ONE+CNST(2.0)*hm%vtau(ip,ispin))*gpsi_i(ip, idir, idim)
            gpsi_j(ip, idir, idim) = (M_ONE+CNST(2.0)*hm%vtau(ip,ispin))*gpsi_j(ip, idir, idim)
          end do
          !$omp end parallel do
        end do
      end do 
     
      !A nonlocal contribution from the pseudopotential must be included
      call zprojector_commute_r_allatoms_alldir(hm%ep%proj, geo, der%mesh, hm%d%dim, &
               der%boundaries, ik, ppsi_i, gpsi_i)                 
      call zprojector_commute_r_allatoms_alldir(hm%ep%proj, geo, der%mesh, hm%d%dim, &
               der%boundaries, ik, ppsi_j, gpsi_j)                 
      !A nonlocal contribution from the scissor must be included
      if(hm%scissor%apply) then
        call scissor_commute_r(hm%scissor, der%mesh, ik, ppsi_i, gpsi_i)
        call scissor_commute_r(hm%scissor, der%mesh, ik, ppsi_j, gpsi_j)
      end if

    end if


    do idir = 1, der%dim
      
      do idim = 1, hm%d%dim
          
        cmel(idir,ispin) = M_zI * zmf_dotp(der%mesh, psi_i(:, idim), gpsi_j(:, idir,idim), reduce = .false.)
        cmel(idir,ispin) = cmel(idir,ispin) - M_zI * zmf_dotp(der%mesh, gpsi_i(:, idir, idim), psi_j(:, idim), reduce = .false.)
          
      end do
    end do

    if(der%mesh%parallel_in_domains) call der%mesh%allreduce(cmel)

    

    SAFE_DEALLOCATE_A(gpsi_i)
    SAFE_DEALLOCATE_A(ppsi_i)
    SAFE_DEALLOCATE_A(gpsi_j)
    SAFE_DEALLOCATE_A(ppsi_j)

    POP_SUB(current_calculate_mel)

  end subroutine current_calculate_mel

  ! ---------------------------------------------------------
  subroutine current_heat_calculate(der, hm, st, current)
    type(derivatives_t),  intent(in)    :: der
    type(hamiltonian_elec_t),  intent(in)    :: hm
    type(states_elec_t),  intent(in)    :: st
    FLOAT,                intent(out)   :: current(:, :, :)

    integer :: ik, ist, idir, idim, ip, ispin, ndim
    CMPLX, allocatable :: gpsi(:, :, :), psi(:, :), g2psi(:, :, :, :)
    CMPLX :: tmp

    PUSH_SUB(current_heat_calculate)

    ASSERT(simul_box_is_periodic(der%mesh%sb))
    ASSERT(st%d%dim == 1)

    ndim = der%dim
    
    SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:der%mesh%np_part, 1:ndim, 1:st%d%dim))
    SAFE_ALLOCATE(g2psi(1:der%mesh%np, 1:ndim, 1:ndim, 1:st%d%dim))
    
    do ip = 1, der%mesh%np
      current(ip, 1:ndim, 1:st%d%nspin) = st%current(ip, 1:ndim, 1:st%d%nspin)*hm%ep%vpsl(ip)
    end do
    
    
    do ik = st%d%kpt%start, st%d%kpt%end
      ispin = st%d%get_spin_index(ik)
      do ist = st%st_start, st%st_end

        if(abs(st%d%kweights(ik)*st%occ(ist, ik)) <= M_EPSILON) cycle
        
        call states_elec_get_state(st, der%mesh, ist, ik, psi)
        do idim = 1, st%d%dim
          call boundaries_set(der%boundaries, psi(:, idim))
        end do

        if (allocated(hm%hm_base%phase)) then 
          call states_elec_set_phase(st%d, psi, hm%hm_base%phase(1:der%mesh%np_part, ik), der%mesh%np_part,  conjugate = .false.)
        end if

        do idim = 1, st%d%dim
          call zderivatives_grad(der, psi(:, idim), gpsi(:, :, idim), set_bc = .false.)
        end do
        do idir = 1, ndim
          if (allocated(hm%hm_base%phase)) then 
            call states_elec_set_phase(st%d, gpsi(:, idir, :), hm%hm_base%phase(1:der%mesh%np_part, ik), der%mesh%np, &
              conjugate = .true.)
          end if
            
          !do idim = 1, st%d%dim
          !  call boundaries_set(der%boundaries, psi(:, idim))
          !end do
           
          do idim = 1, st%d%dim
            call boundaries_set(der%boundaries, gpsi(:,idir, idim))
          end do
            
          if (allocated(hm%hm_base%phase)) then 
            call states_elec_set_phase(st%d, gpsi(:, idir, :), hm%hm_base%phase(1:der%mesh%np_part, ik), &
                                  der%mesh%np_part,  conjugate = .false.)
          end if
            
          do idim = 1, st%d%dim
            call zderivatives_grad(der, gpsi(:, idir, idim), g2psi(:, :, idir, idim), set_bc = .false.)
          end do
        end do
        idim = 1
        do ip = 1, der%mesh%np
          do idir = 1, ndim
            !tmp = sum(conjg(g2psi(ip, idir, 1:ndim, idim))*gpsi(ip, idir, idim)) - sum(conjg(gpsi(ip, 1:ndim, idim))*g2psi(ip, idir, 1:ndim, idim))
            tmp = sum(conjg(g2psi(ip, 1:ndim, idir, idim))*gpsi(ip, 1:ndim, idim)) - &
                  sum(conjg(gpsi(ip, 1:ndim, idim))*g2psi(ip, 1:ndim, idir, idim))
            tmp = tmp - conjg(gpsi(ip, idir, idim))*sum(g2psi(ip, 1:ndim, 1:ndim, idim)) + &
                  sum(conjg(g2psi(ip, 1:ndim, 1:ndim, idim)))*gpsi(ip, idir, idim)
            current(ip, idir, ispin) = current(ip, idir, ispin) + st%d%kweights(ik)*st%occ(ist, ik)*aimag(tmp)/CNST(8.0)
          end do
        end do
      end do
    end do


    POP_SUB(current_heat_calculate)
      
  end subroutine current_heat_calculate

end module current_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
