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
!! $Id: psf.F90 2307 2006-07-29 00:50:22Z appel $

#include "global.h"

module ps_psf_file_m
  use global_m
  use messages_m
  use profiling_m
  use ps_in_grid_m

  implicit none

  private

  public ::                &
    ps_psf_file_t,          &
    ps_psf_file_read,       &
    ps_psf_file_end

  ! First, the contents of the file.
  type ps_psf_file_t
    character(len=2)   :: namatm
    character(len=2)   :: icorr
    character(len=3)   :: irel
    character(len=4)   :: icore
    character(len=10)  :: method(6)
    character(len=70)  :: title

    integer            :: npotd         ! l = 0 .. npotd-1
    integer            :: npotu         ! l = 1 .. npotu
    integer            :: nr
    FLOAT              :: a, b
    FLOAT              :: zval          ! valence charge

    FLOAT, pointer     :: rofi(:)
    FLOAT, pointer     :: vps(:,:)
    FLOAT, pointer     :: chcore(:)
    FLOAT, pointer     :: rho_val(:)
    FLOAT, pointer     :: vso(:,:)
  end type ps_psf_file_t
  
contains

  ! ---------------------------------------------------------
  subroutine ps_psf_file_read(unit, ascii, psf)
    integer,             intent(in)    :: unit
    logical,             intent(in)    :: ascii
    type(ps_psf_file_t), intent(inout) :: psf

    integer  :: ndown, nup, i, l
    character(len=70) :: aux_s

    PUSH_SUB(ps_psf_file_read)

    ! formats used in this routine
8000 format(1x,i2)
9000 format(1x,a2,1x,a2,1x,a3,1x,a4)
9010 format(1x,6a10,/,1x,a70)
9015 format(1x,2i3,i5,3f20.10)
9030 format(4(g20.12))
9040 format(1x,a)

    ! Reads the header line of the file, with general info about the ps.
    if(ascii) then
      read(unit, 9000) psf%namatm, psf%icorr, psf%irel, psf%icore
      read(unit, 9010) (psf%method(l),l=1,6), psf%title
      read(unit, 9015) psf%npotd, psf%npotu, psf%nr, psf%b, psf%a, psf%zval
    else
      read(unit) psf%namatm, psf%icorr, psf%irel, psf%icore,     &
        (psf%method(l), l=1, 6), psf%title, psf%npotd, psf%npotu,  &
        psf%nr, psf%b, psf%a, psf%zval
    end if

    ! add extra point for the zero
    psf%nr = psf%nr + 1

    ! Allocates the variables to psf%nr:  ! Reads the pseudo-valence charge density, in bohr^(-3)
    !   rho_val(1:nrval) : pseudo-valence charge distribution
    SAFE_ALLOCATE(psf%rofi   (1:psf%nr))
    SAFE_ALLOCATE(psf%vps    (1:psf%nr, 1:psf%npotd))
    SAFE_ALLOCATE(psf%chcore (1:psf%nr))
    SAFE_ALLOCATE(psf%rho_val(1:psf%nr))
    SAFE_ALLOCATE(psf%vso    (1:psf%nr, 1:psf%npotu))

    ! Reads the radial values, in bohrs
    !   rofi(1:nr) : radial values ( rofi(i) = b*( exp(a*(i-1)) - 1 ) ) [bohr]
    if(ascii) then
      read(unit, 9040) aux_s
      read(unit, 9030) (psf%rofi(i), i=2, psf%nr)
    else
      read(unit) (psf%rofi(i), i=2, psf%nr)
    end if
    psf%rofi(1) = M_ZERO

    ! Reads the pseudoptential functions, times r, in Rydberg*bohr.
    ! Inmediately afterwards, it is divided by r, so that its final units are Rydbergs
    do ndown = 1, psf%npotd
      if(ascii) then
        read(unit, 9040) aux_s
        read(unit, 8000) l
        read(unit, 9030) (psf%vps(i, ndown), i=2, psf%nr)
      else
        read(unit) l, (psf%vps(i, ndown), i=2, psf%nr)
      end if

      if(l /= ndown-1) then
        message(1) = 'Unexpected angular momentum'
        message(2) = 'Pseudopotential should be ordered by increasing l'
        call messages_warning(2)
      end if
      
      psf%vps(2:, ndown) = psf%vps(2:, ndown) / psf%rofi(2:)
      psf%vps(1,  ndown) = linear_extrapolate(psf%rofi(1), psf%rofi(2), psf%rofi(3), &
        psf%vps(2, ndown), psf%vps(3, ndown))
    end do

    ! Reads --or skips-- the "down" pseudopotentials.
    do nup = 1, psf%npotu
      if(ascii) then
        read(unit, 9040) aux_s
        read(unit, 8000) l
        read(unit, 9030) (psf%vso(i, nup), i=2, psf%nr)
      else
        read(unit) l, (psf%vso(i, nup), i=2, psf%nr)
      end if

      if( (l /= nup) .and. (psf%irel.eq.'rel') ) then
        message(1) = 'Unexpected angular momentum'
        message(2) = 'Pseudopotential should be ordered by increasing l'
        call messages_warning(2)
      end if

      psf%vso(2:, nup) = psf%vso(2:, nup) / psf%rofi(2:)
      psf%vso(1,  nup) = linear_extrapolate(psf%rofi(1), psf%rofi(2), psf%rofi(3), &
        psf%vso(2, nup), psf%vso(3, nup))

    end do
    if(psf%irel.ne.'rel') then
      psf%vso(:,:) = M_ZERO
    end if

    ! Reads the core correcction charge density, in bohr^(-3)
    !   chcore(1:nrval) : core-correction charge distribution
    if(ascii) then
      read(unit, 9040) aux_s
      read(unit, 9030) (psf%chcore(i), i=2, psf%nr)
    else
      read(unit) (psf%chcore(i), i=2, psf%nr)
    end if

    psf%chcore(1) = linear_extrapolate(psf%rofi(1), psf%rofi(2), psf%rofi(3), &
      psf%chcore(2), psf%chcore(3))

    ! Reads the pseudo-valence charge density, in bohr^(-3)
    !   rho_val(1:nrval) : pseudo-valence charge distribution
    !   rho_val(1:nrval) : pseudo-valence charge distribution
    if(ascii) then
      read(unit, 9040) aux_s
      read(unit, 9030) (psf%rho_val(i), i=2, psf%nr)
    else
      read(unit) (psf%rho_val(i), i=2, psf%nr)
    end if

    psf%rho_val(1) = linear_extrapolate(psf%rofi(1), psf%rofi(2), psf%rofi(3), &
      psf%rho_val(2), psf%rho_val(3))

    POP_SUB(ps_psf_file_read)
  end subroutine ps_psf_file_read


  ! ---------------------------------------------------------
  subroutine ps_psf_file_end(psf)
    type(ps_psf_file_t), intent(inout) :: psf

    SAFE_DEALLOCATE_P(psf%rofi)
    SAFE_DEALLOCATE_P(psf%vps)
    SAFE_DEALLOCATE_P(psf%chcore)
    SAFE_DEALLOCATE_P(psf%rho_val)
    SAFE_DEALLOCATE_P(psf%vso)
  end subroutine ps_psf_file_end


end module ps_psf_file_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
