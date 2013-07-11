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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

! ---------------------------------------------------------
subroutine xc_get_fxc(xcs, mesh, rho, ispin, fxc, zfxc)
  type(xc_t), target, intent(in)    :: xcs
  type(mesh_t),       intent(in)    :: mesh
  FLOAT, intent(in)                 :: rho(:, :)
  integer, intent(in)               :: ispin
  FLOAT,              intent(inout) :: fxc(:,:,:)
  CMPLX, intent(inout), optional    :: zfxc(:, :, :, :, :)

  logical :: spinors_fxc
  CMPLX, allocatable :: mmatrix(:, :, :, :), zeigref_(:, :, :)
  FLOAT, allocatable :: l_vdedd(:), vdedd(:, :)
  FLOAT, allocatable :: dens(:,:), dedd(:,:), l_dens(:), l_dedd(:)
  integer :: ip, ixc, spin_channels
  type(xc_functl_t), pointer :: functl(:)

  if(present(zfxc)) then
    ASSERT(ispin  ==  SPINORS)
    spinors_fxc = .true.
    zfxc = M_z0
  else
    ASSERT(ispin /= SPINORS)
    spinors_fxc = .false.
  end if

  if(xcs%kernel_family == XC_FAMILY_NONE) return ! nothing to do

  PUSH_SUB(xc_get_fxc)

  ! is there anything to do? (only LDA by now)
  if(iand(xcs%kernel_family, NOT(XC_FAMILY_LDA)) /= XC_FAMILY_NONE) then
    message(1) = "Only LDA functionals are authorized for now in xc_get_fxc."
    call messages_fatal(1)
  end if

  if(ispin == UNPOLARIZED) then
    functl => xcs%kernel(:, 1)
  else
    functl => xcs%kernel(:, 2)
  end if

  ! This is a bit ugly (why functl(1) and not functl(2)?, but for the moment it works.
  spin_channels = functl(1)%spin_channels
    
  call  lda_init()

  dedd = M_ZERO
  if(spinors_fxc) vdedd = M_ZERO
  space_loop: do ip = 1, mesh%np

    l_dens (:)   = dens (ip, :)

    ! Calculate fxc
    functl_loop: do ixc = 1, 2

      l_dedd = M_ZERO
      if(spinors_fxc) l_vdedd = M_ZERO
      select case(functl(ixc)%family)
      case(XC_FAMILY_LDA)
        call XC_F90(lda_fxc)(functl(ixc)%conf, 1, l_dens(1), l_dedd(1))
        if(spinors_fxc)  call XC_F90(lda_vxc)(functl(ixc)%conf, 1, l_dens(1), l_vdedd(1))
        
      case default
        cycle
      end select

      ! store results
      dedd(ip, :) = dedd(ip, :) + l_dedd(:)
      if(spinors_fxc) vdedd(ip, :) = vdedd(ip, :) + l_vdedd(:)

    end do functl_loop
  end do space_loop

  call  lda_process()
    
  ! clean up allocated memory
  call  lda_end()
  
  POP_SUB(xc_get_fxc)


contains

  ! ---------------------------------------------------------
  ! Takes care of the initialization of the LDA part of the functionals
  !   *) allocates dens(ity) and dedd, and their local variants
  !   *) calculates the density taking into account nlcc and non-collinear spin
  subroutine lda_init()
    integer :: is
    FLOAT   :: d(spin_channels)
    CMPLX :: densitymatrix(2, 2), zeigenval(2)

    PUSH_SUB(xc_get_fxc.lda_init)

    is = 1
    if(ispin == SPIN_POLARIZED .or. ispin == SPINORS) is = 3

    ! allocate some general arrays
    SAFE_ALLOCATE(  dens(1:mesh%np, 1:spin_channels))
    SAFE_ALLOCATE(  dedd(1:mesh%np, 1:is))
    SAFE_ALLOCATE(l_dens(1:spin_channels))
    SAFE_ALLOCATE(l_dedd(1:is))
    dedd = M_ZERO

    if(spinors_fxc) then
      SAFE_ALLOCATE(l_vdedd(1:spin_channels))
      SAFE_ALLOCATE(vdedd(1:mesh%np, 1:spin_channels))
      SAFE_ALLOCATE(zeigref_(1:2, 1:2, 1:mesh%np))
      SAFE_ALLOCATE(mmatrix(1:2, 1:2, 1:2, 1:mesh%np))
    end if

    ! get the density
    do ip = 1, mesh%np
      d(1:spin_channels) = rho(ip, 1:spin_channels)

      select case(ispin)
      case(UNPOLARIZED)
        dens(ip, 1) = max(d(1), M_ZERO)
      case(SPIN_POLARIZED)
        dens(ip, 1) = max(d(1), M_ZERO)
        dens(ip, 2) = max(d(2), M_ZERO)
      case(SPINORS)

        densitymatrix(1, 1) = rho(ip, 1)
        densitymatrix(2, 2) = rho(ip, 2)
        densitymatrix(1, 2) = rho(ip, 3) + M_zI * rho(ip, 4)
        densitymatrix(2, 1) = rho(ip, 3) - M_zI * rho(ip, 4)

        !CHECK OF DERIVATIVES OF EIGENVALUES
        !densitymatrix(1, 1) = ( CNST(0.1), M_ZERO)
        !densitymatrix(2, 2) = (-CNST(0.1), M_ZERO)
        !densitymatrix(1, 2) = ( CNST(0.2), CNST(0.25))
        !densitymatrix(2, 1) = ( CNST(0.2),-CNST(0.25))
        !call lalg_check_zeigenderivatives(2, densitymatrix)
        !stop
        !ENDOFCHECK
        if(maxval(abs(densitymatrix)) < tiny) densitymatrix = M_z0

        call lalg_zeigenderivatives(2, densitymatrix, zeigref_(:, :, ip), zeigenval, mmatrix(:, :, :, ip))
        dens(ip, 1) = max(real(zeigenval(1), REAL_PRECISION), M_ZERO)
        dens(ip, 2) = max(real(zeigenval(2), REAL_PRECISION), M_ZERO)

      end select
    end do

    POP_SUB(xc_get_fxc.lda_init)
  end subroutine lda_init


  ! ---------------------------------------------------------
  ! deallocates variables allocated in lda_init
  subroutine lda_end()
    PUSH_SUB(xc_get_fxc.lda_end)

    SAFE_DEALLOCATE_A(dens)
    SAFE_DEALLOCATE_A(dedd)
    SAFE_DEALLOCATE_A(l_dens)
    SAFE_DEALLOCATE_A(l_dedd)
    if(ispin == SPINORS) then
      SAFE_DEALLOCATE_A(l_vdedd)
      SAFE_DEALLOCATE_A(mmatrix)
    end if

    POP_SUB(xc_get_fxc.lda_end)
  end subroutine lda_end


  ! ---------------------------------------------------------
  ! calculates the LDA part of vxc, taking into account non-collinear spin
  subroutine lda_process()
    integer :: alpha, beta, delta, gamma, ip, i, j
    FLOAT :: localfxc(2, 2)
    PUSH_SUB(xc_get_fxc.lda_process)

    select case(ispin)
    case(UNPOLARIZED)
      fxc(:,1,1) = fxc(:,1,1) + dedd(:,1)
    case(SPIN_POLARIZED)
      fxc(:,1,1) = fxc(:,1,1) + dedd(:,1)
      fxc(:,2,2) = fxc(:,2,2) + dedd(:,3)
      fxc(:,1,2) = fxc(:,1,2) + dedd(:,2)
      fxc(:,2,1) = fxc(:,2,1) + dedd(:,2)
    case(SPINORS)
      do ip = 1, mesh%np
        localfxc(1, 1) = dedd(ip, 1)
        localfxc(1, 2) = dedd(ip, 2)
        localfxc(2, 1) = dedd(ip, 2)
        localfxc(2, 2) = dedd(ip, 3)
        do alpha = 1, 2
          do beta = 1, 2
            do gamma = 1, 2
              do delta = 1, 2
                do i = 1, 2
                  zfxc(alpha, beta, gamma, delta, ip) = zfxc(alpha, beta, gamma, delta, ip) + &
                    lalg_zd2ni(zeigref_(1:2, i, ip), mmatrix(:, :, i, ip), beta, alpha, delta, gamma) * &
                    vdedd(ip, i)
                  do j = 1, 2
                    zfxc(alpha, beta, gamma, delta, ip) = zfxc(alpha, beta, gamma, delta, ip) + &
                       lalg_zdni(zeigref_(1:2, i, ip), beta, alpha) * &
                       lalg_zdni(zeigref_(1:2, j, ip), delta, gamma) * localfxc(j, i)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end select

    POP_SUB(xc_get_fxc.lda_process)
  end subroutine lda_process

end subroutine xc_get_fxc

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
