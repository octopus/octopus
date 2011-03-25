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
!! $Id$


! ---------------------------------------------------------
! this is for the third derivative of Exc
subroutine xc_get_kxc(xcs, mesh, rho, ispin, kxc)
  type(xc_t), target, intent(in)    :: xcs
  type(mesh_t),       intent(in)    :: mesh
  FLOAT,              intent(in)    :: rho(:, :)
  integer,            intent(in)    :: ispin
  FLOAT,              intent(inout) :: kxc(:,:,:,:)

  FLOAT, allocatable :: dens(:,:), dedd(:,:), l_dens(:), l_dedd(:)
  integer :: i, ixc, spin_channels
  type(xc_functl_t), pointer :: functl(:)

  PUSH_SUB(xc_get_kxc)

  if(ispin == UNPOLARIZED) then
    functl => xcs%kernel(:, 1)
  else
    functl => xcs%kernel(:, 2)
  end if

  ! is there anything to do? (only LDA by now)
  if(iand(xcs%kernel_family, NOT(XC_FAMILY_LDA)).ne.XC_FAMILY_NONE) then
    message(1) = "Only LDA functionals are authorized for now in xc_get_kxc."
    call messages_fatal(1)
  end if

  if(xcs%kernel_family == XC_FAMILY_NONE) then
    POP_SUB(xc_get_kxc)
    return
  endif

  ! really start

  ! This is a bit ugly (why functl(1) and not functl(2)?, but for the moment it works.
  spin_channels = functl(1)%spin_channels

  call lda_init()

  space_loop: do i = 1, mesh%np

    ! make a local copy with the correct memory order
    l_dens (:)   = dens (i, :)

    ! Calculate kxc
    functl_loop: do ixc = 1, 2

      select case(functl(ixc)%family)
      case(XC_FAMILY_LDA)
        call XC_F90(lda_kxc)(functl(ixc)%conf, 1, l_dens(1), l_dedd(1))

      case default
        cycle
      end select

      ! store results
      dedd(i,:) = dedd(i,:) + l_dedd(:)

    end do functl_loop
  end do space_loop

  call lda_process()

  ! clean up allocated memory
  call lda_end()

  POP_SUB(xc_get_kxc)

contains

  ! ---------------------------------------------------------
  ! Takes care of the initialization of the LDA part of the functionals
  !   *) allocates dens(ity) and dedd, and their local variants
  !   *) calculates the density taking into account nlcc and non-collinear spin
  subroutine lda_init()
    integer :: i
    FLOAT   :: d(spin_channels)

    PUSH_SUB(xc_get_kxc.lda_init)

    ! allocate some general arrays
    SAFE_ALLOCATE(  dens(1:mesh%np, 1:spin_channels))
    SAFE_ALLOCATE(  dedd(1:mesh%np, 1:spin_channels*spin_channels))
    SAFE_ALLOCATE(l_dens(1:spin_channels))
    SAFE_ALLOCATE(l_dedd(1:spin_channels*spin_channels))
    dedd = M_ZERO

    ! get the density
    do i = 1, mesh%np
      d(:) = rho(i, :)

      select case(ispin)
      case(UNPOLARIZED)
        dens(i, 1) = max(d(1), M_ZERO)
      case(SPIN_POLARIZED)
        dens(i, 1) = max(d(1), M_ZERO)
        dens(i, 2) = max(d(2), M_ZERO)
      case(SPINORS)
        message(1) = 'Do not know how to handle spinors.'
        call messages_fatal(1)
      end select
    end do

    POP_SUB(xc_get_kxc.lda_init)
  end subroutine lda_init


  ! ---------------------------------------------------------
  ! SAFE_DEALLOCATE_Ps variables allocated in lda_init
  subroutine lda_end()
    PUSH_SUB(xc_get_kxc.lda_end)

    SAFE_DEALLOCATE_A(dens)
    SAFE_DEALLOCATE_A(dedd)
    SAFE_DEALLOCATE_A(l_dens)
    SAFE_DEALLOCATE_A(l_dedd)

    POP_SUB(xc_get_kxc.lda_end)
  end subroutine lda_end


  ! ---------------------------------------------------------
  ! calculates the LDA part of vxc
  subroutine lda_process()
    integer :: i

    PUSH_SUB(xc_get_kxc.lda_process)

    forall(i = 1:mesh%np) kxc(i,1,1,1) = kxc(i,1,1,1) + dedd(i,1)

    if(ispin == SPIN_POLARIZED) then
      forall(i = 1:mesh%np)
        kxc(i,1,1,2) = kxc(i,1,1,2) + dedd(i,2)
        kxc(i,1,2,1) = kxc(i,1,2,1) + dedd(i,2)
        kxc(i,1,2,2) = kxc(i,1,2,2) + dedd(i,3)
        kxc(i,2,1,1) = kxc(i,2,1,1) + dedd(i,2)
        kxc(i,2,1,2) = kxc(i,2,1,2) + dedd(i,3)
        kxc(i,2,2,1) = kxc(i,2,2,1) + dedd(i,3)
        kxc(i,2,2,2) = kxc(i,2,2,2) + dedd(i,4)
      end forall
    end if

    POP_SUB(xc_get_kxc.lda_process)
  end subroutine lda_process

end subroutine xc_get_kxc


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
