!! Copyright (C) 2015 X. Andrade
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

!---------------------------------------------------------------------------
subroutine X(density_accumulate_grad)(gr, st, psib, grad_psib, grad_rho)
  type(grid_t),        intent(in)    :: gr
  type(states_elec_t), intent(in)    :: st
  type(wfs_elec_t),    intent(in)    :: psib
  type(wfs_elec_t),    intent(in)    :: grad_psib(:)
  FLOAT,               intent(inout) :: grad_rho(:, :)

  integer :: ii, ist, idir, ip
  FLOAT :: ff
  R_TYPE :: psi, gpsi
  integer :: wgsize
  FLOAT, allocatable :: grad_rho_tmp(:, :), weights(:)
  type(accel_mem_t) :: grad_rho_buff, weights_buff
  type(accel_kernel_t), save :: ker_calc_grad_dens

  PUSH_SUB(X(density_accumulate_grad))

  call psib%check_compatibility_with(grad_psib(1))

  select case(psib%status())
  case(BATCH_NOT_PACKED)
    do idir = 1, gr%mesh%sb%dim
      do ii = 1, psib%nst_linear
        ist = psib%linear_to_ist(ii)
      
        ff = st%d%kweights(psib%ik)*st%occ(ist, psib%ik)*M_TWO
        if(abs(ff) <= M_EPSILON) cycle

        do ip = 1, gr%mesh%np
          
          psi = psib%X(ff_linear)(ip, ii)
          gpsi = grad_psib(idir)%X(ff_linear)(ip, ii)
          grad_rho(ip, idir) = grad_rho(ip, idir) + ff*R_REAL(R_CONJ(psi)*gpsi)
          
        end do
      end do
      
    end do

  case(BATCH_PACKED)
    do ii = 1, psib%nst_linear
      ist = psib%linear_to_ist(ii)
      
      ff = st%d%kweights(psib%ik)*st%occ(ist, psib%ik)*M_TWO
      if(abs(ff) <= M_EPSILON) cycle

      do idir = 1, gr%mesh%sb%dim
        do ip = 1, gr%mesh%np
          
          psi = psib%X(ff_pack)(ii, ip)
          gpsi = grad_psib(idir)%X(ff_pack)(ii, ip)
          grad_rho(ip, idir) = grad_rho(ip, idir) + ff*R_REAL(R_CONJ(psi)*gpsi)
          
        end do
      end do
      
    end do
      
  case(BATCH_DEVICE_PACKED)
    call accel_create_buffer(grad_rho_buff, ACCEL_MEM_WRITE_ONLY, TYPE_FLOAT, gr%mesh%np*gr%sb%dim)

    SAFE_ALLOCATE(weights(1:psib%pack_size(1)))

    weights = CNST(0.0)
    do ii = 1, psib%nst_linear
      ist = psib%linear_to_ist(ii)
      weights(ii) = st%d%kweights(psib%ik)*st%occ(ist, psib%ik)*M_TWO
    end do

    call accel_create_buffer(weights_buff, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, psib%pack_size(1))
    call accel_write_buffer(weights_buff, psib%pack_size(1), weights)
   
    SAFE_DEALLOCATE_A(weights)
    
    call accel_kernel_start_call(ker_calc_grad_dens, 'forces.cl', trim(TOSTRING(X(density_gradient))), &
      flags = '-D' + R_TYPE_CL)

    do idir = 1, gr%mesh%sb%dim
      call accel_set_kernel_arg(ker_calc_grad_dens, 0, idir - 1)
      call accel_set_kernel_arg(ker_calc_grad_dens, 1, psib%pack_size(1))
      call accel_set_kernel_arg(ker_calc_grad_dens, 2, gr%mesh%np)
      call accel_set_kernel_arg(ker_calc_grad_dens, 3, weights_buff)
      call accel_set_kernel_arg(ker_calc_grad_dens, 4, grad_psib(idir)%ff_device)
      call accel_set_kernel_arg(ker_calc_grad_dens, 5, psib%ff_device)
      call accel_set_kernel_arg(ker_calc_grad_dens, 6, log2(psib%pack_size(1)))
      call accel_set_kernel_arg(ker_calc_grad_dens, 7, grad_rho_buff)

      wgsize = accel_kernel_workgroup_size(ker_calc_grad_dens)
      call accel_kernel_run(ker_calc_grad_dens, (/pad(gr%mesh%np, wgsize)/), (/wgsize/))

    end do

    call accel_release_buffer(weights_buff)

    SAFE_ALLOCATE(grad_rho_tmp(1:gr%mesh%np, 1:gr%sb%dim))

    call accel_read_buffer(grad_rho_buff, gr%mesh%np*gr%sb%dim, grad_rho_tmp)

    call accel_release_buffer(grad_rho_buff)

    do idir = 1, gr%mesh%sb%dim
      do ip = 1, gr%mesh%np
        grad_rho(ip, idir) = grad_rho(ip, idir) + grad_rho_tmp(ip, idir)
      end do
    end do

    SAFE_DEALLOCATE_A(grad_rho_tmp)
  end select

  POP_SUB(X(density_accumulate_grad))

end subroutine X(density_accumulate_grad)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
