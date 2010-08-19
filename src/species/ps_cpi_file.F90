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
!! $Id: tm.F90 2307 2006-07-29 00:50:22Z appel $

#include "global.h"

module ps_cpi_file_m
  use global_m
  use messages_m
  use profiling_m
  use ps_in_grid_m

  implicit none

  private

  public ::                &
    ps_cpi_file_t,          &
    ps_cpi_file_read,       &
    ps_cpi_file_end

  ! First, the contents of the file.
  type ps_cpi_file_t
    FLOAT              :: zval          ! valence charge
    integer            :: no_l_channels ! number of pseudo components (lmax+1)
    
    integer            :: nr            ! number of mesh points
    FLOAT              :: a             ! mesh multiplicative increment
    
    FLOAT, pointer     :: rofi(:)       ! radial mesh
    FLOAT, pointer     :: vps(:,:)      ! pseudopotential
    FLOAT, pointer     :: rphi(:,:)     ! r times the pseudowavefunctions
    
    logical            :: core_corrections
    FLOAT, pointer     :: chcore(:)     ! r times the core charge
    FLOAT, pointer     :: d1chcore(:)   ! first  derivative of chcore
    FLOAT, pointer     :: d2chcore(:)   ! second derivative of chcore
  end type ps_cpi_file_t
  
contains

  ! ---------------------------------------------------------
  subroutine ps_cpi_file_read(unit, psf)
    integer,             intent(in)    :: unit
    type(ps_cpi_file_t), intent(inout) :: psf

    integer  :: i, l, ios, idummy
    FLOAT    :: a, b, c, d

    PUSH_SUB(read_file_data)

    read(unit, *) psf%zval, psf%no_l_channels
    ! skip 10 lines
    do i = 1, 10
      read(unit, *)
    end do

    read(unit, *) psf%nr, psf%a

    ! add extra point for the zero
    psf%nr = psf%nr + 1

    SAFE_ALLOCATE(psf%rofi   (1:psf%nr))
    SAFE_ALLOCATE(psf%vps    (1:psf%nr, 1:psf%no_l_channels))
    SAFE_ALLOCATE(psf%rphi   (1:psf%nr, 1:psf%no_l_channels))

    do l = 1, psf%no_l_channels
      if(l.ne.1) read(unit, *)

      do i = 2, psf%nr
        read(unit, *) idummy, psf%rofi(i), psf%rphi(i, l), psf%vps(i, l)
      end do
    end do

    ! read core charge (if present)
    read(unit, *, iostat=ios) a, b, c, d
    if(ios == 0) then
      psf%core_corrections = .true.

      SAFE_ALLOCATE(psf%chcore  (1:psf%nr))
      SAFE_ALLOCATE(psf%d1chcore(1:psf%nr))
      SAFE_ALLOCATE(psf%d2chcore(1:psf%nr))

      psf%  chcore(2) = b
      psf%d1chcore(2) = c
      psf%d2chcore(2) = d

      do i = 3, psf%nr
        read(unit, *) a, psf%chcore(i), psf%d1chcore(i), psf%d2chcore(i)
      end do
    else
      psf%core_corrections = .false.
    end if

    ! add extra point at zero
    psf%rofi(1) = M_ZERO
    do l = 1, psf%no_l_channels
      psf%vps(1,  l) = linear_extrapolate(psf%rofi(1), psf%rofi(2), psf%rofi(3), &
        psf%vps(2, l), psf%vps(3, l))

      psf%rphi(1, l) = M_ZERO
    end do

    if(psf%core_corrections) then
      ! At this point, we use the normalization of the siesta format, where
      ! psf%chcore(:) = 4*pi*\tilde{rho} r**2. As in the Fritz-Haber file we
      ! have written 4*pi*\tilde{rho}, we multiply by r**2
      psf%chcore(:) = psf%chcore(:) * psf%rofi(:)**2

      psf%chcore(1) = linear_extrapolate(psf%rofi(1), psf%rofi(2), psf%rofi(3), &
        psf%chcore(2), psf%chcore(3))

      psf%d1chcore(1) = linear_extrapolate(psf%rofi(1), psf%rofi(2), psf%rofi(3), &
        psf%d1chcore(2), psf%d1chcore(3))
      

      psf%d2chcore(1) = linear_extrapolate(psf%rofi(1), psf%rofi(2), psf%rofi(3), &
        psf%d2chcore(2), psf%d2chcore(3))
    end if

    ! WARNING: This should go away
    psf%vps(:,:) = psf%vps(:,:)*M_TWO ! convert to Rydbergs
    POP_SUB(read_file_data)
  end subroutine ps_cpi_file_read


  ! ---------------------------------------------------------
  subroutine ps_cpi_file_end(psf)
    type(ps_cpi_file_t), intent(inout) :: psf

    SAFE_DEALLOCATE_P(psf%rofi)
    SAFE_DEALLOCATE_P(psf%vps)
    SAFE_DEALLOCATE_P(psf%rphi)

    if(psf%core_corrections) then
      SAFE_DEALLOCATE_P(psf%chcore)
      SAFE_DEALLOCATE_P(psf%d1chcore)
      SAFE_DEALLOCATE_P(psf%d2chcore)
    end if

  end subroutine ps_cpi_file_end

end module ps_cpi_file_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
