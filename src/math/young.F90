!! Copyright (C) 2009
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

#include "global.h"

module young_m
  use global_m
  use loct_m
  use math_m
  use messages_m
  use profiling_m

  implicit none

  private

  public :: young_init,           &
            young_write_allspins, &
            young_write,          &
            young_write_one,      &
            young_copy,           &
            young_nullify,        &
            young_end,            &
            young_t

  type young_t 
   integer :: nyoung, nup, ndown, iyoung
   integer, pointer :: young_up(:,:)
   integer, pointer :: young_down(:,:)
  end type young_t 

  contains

  !------------------------------------------------------------
  subroutine young_init (this, nup, ndown)
    integer, intent(in) :: nup, ndown
    type(young_t), intent(inout) :: this

    ! local vars
    integer :: ipart

    PUSH_SUB(young_init)

    if (ndown > nup) then
      write (message(1),'(a)') 'Error: we only make 2-row Young diagrams with nup >= ndown'
      call messages_fatal(1)
    end if

    this%nup = nup
    this%ndown = ndown

    this%nyoung = factorial(nup+ndown)
    do ipart = 1, ndown
      ! hook factor for down spin
      this%nyoung = this%nyoung / (ndown-ipart+1)
      ! hook factor for up spins which are paired to a down spin
      this%nyoung = this%nyoung / (nup  -ipart+2)
    end do
    do ipart = ndown+1, nup
      ! hook factor for unpaired up spins
      this%nyoung = this%nyoung / (nup  -ipart+1)
    end do

    SAFE_ALLOCATE (this%young_up  (1:nup,  1:this%nyoung))
    SAFE_ALLOCATE (this%young_down(1:ndown,1:this%nyoung))

    this%young_up(:, :)   = -999
    this%young_down(:, :) = -999

    this%iyoung = 1
    call young_fill (this, nup+ndown)

    POP_SUB(young_init)
  end subroutine young_init


  !------------------------------------------------------------
  recursive subroutine young_fill (this, nn)
    integer, intent(in) :: nn
    type(young_t), intent(inout) :: this

    ! local vars
    integer :: idown, iup
    PUSH_SUB(young_fill)

    if (this%iyoung > this%nyoung) then
    POP_SUB(young_fill)
      return
    end if

    ! find next lower right hand corner, in the down spins
    do idown = this%ndown, 1, -1
      if (this%young_down(idown, this%iyoung) == -999) then
        this%young_down(idown, this%iyoung) = nn
        ! call again with smaller diagram
        if (nn > 1) then
          call young_fill (this, nn-1)
        else
          this%iyoung = this%iyoung+1
          if (this%iyoung <= this%nyoung) then
            this%young_up(:,this%iyoung)   = this%young_up(:,this%iyoung-1)
            this%young_down(:,this%iyoung) = this%young_down(:,this%iyoung-1)
            call young_reset_1val (this, 0)
          end if
        end if
        exit
      end if
    end do

    if (this%iyoung > this%nyoung) then
    POP_SUB(young_fill)
      return
    end if

    ! find next lower right hand corner, in the up spins
    do iup = this%nup, 1, -1
      if (this%young_up(iup, this%iyoung) == -999) then
        ! either in the unpaired spins
        if (iup > this%ndown) then
          this%young_up(iup, this%iyoung) = nn
          ! call again with smaller diagram
          if (nn > 1) then
            call young_fill (this, nn-1)
          else
            this%iyoung = this%iyoung+1
            if (this%iyoung <= this%nyoung) then
              this%young_up(:,this%iyoung)   = this%young_up(:,this%iyoung-1)
              this%young_down(:,this%iyoung) = this%young_down(:,this%iyoung-1)
              call young_reset_1val (this, 0)
            end if
          end if
        ! or in the paired spins, provided the box below has been filled
        else if (this%young_down(iup, this%iyoung) /= -999) then
          this%young_up(iup, this%iyoung) = nn
          ! call again with smaller diagram
          if (nn > 1) then
            call young_fill (this, nn-1)
          else
            this%iyoung = this%iyoung+1
            if (this%iyoung <= this%nyoung) then
              this%young_up(:,this%iyoung)   = this%young_up(:,this%iyoung-1)
              this%young_down(:,this%iyoung) = this%young_down(:,this%iyoung-1)
              call young_reset_1val (this, 0)
            end if
          end if
        end if
        exit
      end if
    end do

    if (this%iyoung > this%nyoung) then
    POP_SUB(young_fill)
      return
    end if
    
    call young_reset_1val (this, nn)

!
!    if (all(this%young_up(:,this%iyoung) /= nn) .and. all(this%young_down(:,this%iyoung) /= nn)) then
!      write (message(1),'(a,I7,a)') 'Error : nn = ', nn, &
!            ' was not attributed to any box- this should not happen!'
!      call messages_fatal(1)
!    end if

    POP_SUB(young_fill)
  end subroutine


  !------------------------------------------------------------
  subroutine young_reset_1val (this, nn)
    type(young_t), intent(inout) :: this
    integer, intent(in) :: nn

    integer :: iup, idown

    PUSH_SUB(young_reset_1val)

    ! remove last entry in new diagram
    do iup = 1, this%nup
      if(this%young_up(iup,this%iyoung) == nn+1) this%young_up(iup,this%iyoung) = -999
    end do
    do idown = 1, this%ndown
      if(this%young_down(idown,this%iyoung) == nn+1) this%young_down(idown,this%iyoung) = -999
    end do

    POP_SUB(young_reset_1val)
  end subroutine young_reset_1val


  !------------------------------------------------------------
  subroutine young_write (iunit, this)
    integer, intent(in) :: iunit
    type(young_t), intent(inout) :: this

    integer :: iyoung

    PUSH_SUB(young_write)
    
    write (iunit, '(a,I4,a,I4,a)') ' Young diagrams for ', this%nup, ' spins up, and ', this%ndown, ' down '
    do iyoung = 1, this%nyoung
      call young_write_one (iunit, this, iyoung)
    end do
    
    POP_SUB(young_write)
  end subroutine young_write


  !------------------------------------------------------------
  subroutine young_write_one (iunit, this, iyoung)
    integer, intent(in) :: iunit
    type(young_t), intent(inout) :: this

    integer, intent(in) :: iyoung

    PUSH_SUB(young_write_one)
    
    write (iunit,'(a,I7)') ' Young diagram ', iyoung
    write (iunit,'(10I7)') this%young_up(:, iyoung)
    write (iunit,'(10I7)') this%young_down(:, iyoung)
    
    POP_SUB(young_write_one)
  end subroutine young_write_one


  !------------------------------------------------------------
  ! routine gets all Young diagrams for all distributions of spins
  subroutine young_write_allspins (iunit, nparticles)
    integer, intent(in) :: iunit, nparticles
    integer :: nup, ndown
    type(young_t) :: this

    PUSH_SUB(young_write_allspins)

    call young_nullify (this)
    do ndown = 0, floor(nparticles * 0.5)
      nup = nparticles - ndown
      call young_init (this, nup, ndown)
      call young_write (iunit, this)
      call young_end (this)
    end do
    POP_SUB(young_write_allspins)

  end subroutine young_write_allspins 


  !------------------------------------------------------------
  subroutine young_nullify (this)
    type(young_t), intent(inout) :: this

    PUSH_SUB(young_nullify)

    nullify(this%young_up)
    nullify(this%young_down)

    POP_SUB(young_nullify)
  end subroutine young_nullify


  !------------------------------------------------------------
  subroutine young_copy (young_in, young_out)
    type(young_t), intent(inout) :: young_in, young_out

    PUSH_SUB(young_copy)

    young_out%nup = young_in%nup
    young_out%ndown = young_in%ndown
    young_out%nyoung = young_in%nyoung

    call loct_pointer_copy(young_out%young_up,young_in%young_up)
    call loct_pointer_copy(young_out%young_down,young_in%young_down)

    POP_SUB(young_copy)
  end subroutine young_copy


  !------------------------------------------------------------
  subroutine young_end (this)
    type(young_t), intent(inout) :: this

    PUSH_SUB(young_end)

    SAFE_DEALLOCATE_P (this%young_up)
    SAFE_DEALLOCATE_P (this%young_down)

    POP_SUB(young_end)
  end subroutine young_end

end module young_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
