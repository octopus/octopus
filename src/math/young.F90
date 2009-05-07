


#include "global.h"

module young_m

  use global_m
  use loct_m
  use math_m
  use messages_m
  use profiling_m

  implicit none

  private

  public :: young_init, &
            young_write, &
            young_copy, &
            young_nullify, &
            young_end, &
            young_t

  type young_t 
   integer :: nyoung, nup, ndown, iyoung
   integer, pointer :: young_up(:,:)
   integer, pointer :: young_down(:,:)
  end type young_t 



  contains

  subroutine young_init (this, nup, ndown)

    implicit none

    integer, intent(in) :: nup, ndown
    type(young_t), intent(inout) :: this

    ! local vars
    integer :: ipart

    call push_sub('math.young_init')

    if (ndown > nup) then
      write (message(1),'(a)') 'Error: we only make 2 row Young diagrams with nup >= ndown'
      call write_fatal(1)
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

    write (message(1), '(a,I7)') 'young_init found nyoung = ', this%nyoung
    call write_info(1)

    SAFE_ALLOCATE (this%young_up  (1:nup,  1:this%nyoung))
    SAFE_ALLOCATE (this%young_down(1:ndown,1:this%nyoung))

    this%young_up(:, :)   = -999
    this%young_down(:, :) = -999

    this%iyoung = 1
    call young_fill (this, nup+ndown)

    call pop_sub()
  end subroutine young_init



  recursive subroutine young_fill (this, nn)
    implicit none

    integer, intent(in) :: nn
    type(young_t), intent(inout) :: this

    ! local vars
    integer :: idown, iup
    call push_sub('math.young_fill')

    if (this%iyoung > this%nyoung) return

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

    if (this%iyoung > this%nyoung) return

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

    if (this%iyoung > this%nyoung) return
    
    call young_reset_1val (this, nn)

!
!    if (all(this%young_up(:,this%iyoung) /= nn) .and. all(this%young_down(:,this%iyoung) /= nn)) then
!      write (message(1),'(a,I7,a)') 'Error : nn = ', nn, &
!            ' was not attributed to any box- this should not happen!'
!      call write_fatal(1)
!    end if

    call pop_sub()
  end subroutine

  subroutine young_reset_1val (this, nn)
    implicit none

    type(young_t), intent(inout) :: this
    integer, intent(in) :: nn

    integer :: iup, idown

    ! remove last entry in new diagram
    do iup = 1, this%nup
      if(this%young_up(iup,this%iyoung) == nn+1) this%young_up(iup,this%iyoung) = -999
    end do
    do idown = 1, this%ndown
      if(this%young_down(idown,this%iyoung) == nn+1) this%young_down(idown,this%iyoung) = -999
    end do

  end subroutine young_reset_1val

  subroutine young_write (this)
    implicit none

    type(young_t), intent(inout) :: this
    ! local vars
    integer :: iyoung
    call push_sub('math.young_write')
    
    do iyoung = 1, this%nyoung
      write (message(1),'(a,I7)') 'Filled Young diagram ', iyoung
      write (message(2),'(10I7)') this%young_up(:, iyoung)
      write (message(3),'(10I7)') this%young_down(:, iyoung)
      call write_info(3)
    end do
    
    call pop_sub()
  end subroutine young_write

  subroutine young_nullify (this)
    implicit none

    type(young_t), intent(inout) :: this
    call push_sub('math.young_nullify')
    nullify(this%young_up)
    nullify(this%young_down)
    call pop_sub()
  end subroutine young_nullify

  subroutine young_copy (young_in, young_out)
    implicit none

    type(young_t), intent(inout) :: young_in, young_out

    call push_sub('math.young_copy')

    young_out%nup = young_in%nup
    young_out%ndown = young_in%ndown
    young_out%nyoung = young_in%nyoung

    call loct_pointer_copy(young_out%young_up,young_in%young_up)
    call loct_pointer_copy(young_out%young_down,young_in%young_down)
    call pop_sub()
  end subroutine young_copy

  subroutine young_end (this)
    implicit none

    type(young_t), intent(inout) :: this

    call push_sub('math.young_end')
    SAFE_DEALLOCATE_P (this%young_up)
    SAFE_DEALLOCATE_P (this%young_down)
    call pop_sub()
  end subroutine young_end

end module young_m
