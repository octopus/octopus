!! Copyright (C) 2008 X. Andrade
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
!! $Id: multicomm.F90 4396 2008-07-21 16:14:17Z xavier $

#include "global.h"

module distributed_m
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use parser_m
  use messages_m
  use mpi_m
  use multicomm_m
  use profiling_m
  use utils_m
  use varinfo_m

  implicit none

  private

  public                            &
       distributed_t,               &
       distributed_nullify,         &
       distributed_init,            &
       distributed_copy,            &
       distributed_end
  

  type distributed_t
    integer :: start
    integer :: end
    integer :: nlocal
    integer :: nglobal
    logical :: parallel
    integer, pointer :: node(:)
    integer, pointer :: range(:, :)
    integer, pointer :: num(:)
    type(mpi_grp_t)  :: mpi_grp
  end type distributed_t
  
contains

  subroutine distributed_nullify(this, total)
    type(distributed_t), intent(out) :: this
    integer,             intent(in)  :: total

    call push_sub('distributed.distributed_nullify')

    this%start           = 1
    this%end             = total
    this%nlocal          = total
    this%nglobal         = total
    this%parallel        = .false.
    nullify(this%node, this%range, this%num)
    call mpi_grp_init(this%mpi_grp, -1)

    call pop_sub()
  end subroutine distributed_nullify
  
  subroutine distributed_init(this, total, mc, strategy, tag)
    type(distributed_t), intent(out) :: this
    integer,             intent(in)  :: total
    type(multicomm_t),   intent(in)  :: mc
    integer,             intent(in)  :: strategy
    character(len=*),    intent(in)  :: tag

#ifdef HAVE_MPI
    integer :: kk, size, rank
#endif

    call push_sub('distributed.distributed_init')
    
    this%nglobal = total

#ifdef HAVE_MPI
    if(  (.not. multicomm_strategy_is_parallel(mc, strategy)) .or. &
         (this%nglobal.eq.1) ) then
#endif      
      
      SAFE_ALLOCATE(this%node(1:total))
      ! Defaults.
      this%node(1:total)   = 0
      this%start           = 1
      this%end             = total
      this%nlocal          = total
      this%nglobal         = total
      this%parallel        = .false.
      nullify(this%range, this%num)
      call mpi_grp_init(this%mpi_grp, -1)
      
#ifdef HAVE_MPI
    else

      this%parallel = .true.

      call mpi_grp_init(this%mpi_grp, mc%group_comm(strategy))

      size = mc%group_sizes(strategy)
      rank = mc%who_am_i(strategy)

      SAFE_ALLOCATE(this%range(1:2, 0:size - 1))
      SAFE_ALLOCATE(this%num(0:size - 1))
      SAFE_ALLOCATE(this%node(1:this%nglobal))

      call multicomm_divide_range(this%nglobal, size, this%range(1, :), this%range(2, :), this%num)

      message(1) = 'Info: Parallelization in ' // trim(tag)
      call write_info(1)

      do kk = 1, size
        write(message(1),'(a,i4,a,i6,a,i6)') 'Info: Node in group ', kk - 1, &
             ' will manage '//trim(tag), this%range(1, kk - 1), " - ", this%range(2, kk - 1)
        call write_info(1)
        if(rank .eq. kk - 1) then
          this%start  = this%range(1, kk - 1)
          this%end    = this%range(2, kk - 1)
          this%nlocal = this%num(kk - 1)
        endif
        
        this%node(this%range(1, kk - 1):this%range(2, kk - 1)) = kk - 1
        
      end do
      
    end if
#endif
    
    call pop_sub()
  end subroutine distributed_init

  subroutine distributed_copy(in, out)
    type(distributed_t), intent(in)  :: in
    type(distributed_t), intent(out) :: out

    integer :: size

    call push_sub('distributed.distributed_copy')

    out%start    = in%start
    out%end      = in%end
    out%nlocal   = in%nlocal
    out%nglobal  = in%nglobal
    out%parallel = in%parallel

    size = in%mpi_grp%size

    call mpi_grp_init(out%mpi_grp, in%mpi_grp%comm)
    
    nullify(out%node, out%range, out%num)

    if(associated(in%node)) then
      SAFE_ALLOCATE(out%node(1:in%nglobal))
      out%node(1:in%nglobal) = in%node(1:in%nglobal)
    end if

    if(associated(in%range)) then
      SAFE_ALLOCATE(out%range(1:2, 0:size - 1))
      out%range(1:2, 0:size - 1) = in%range(1:2, 0:size - 1)
    end if

    if(associated(in%num)) then
      SAFE_ALLOCATE(out%num(0:size - 1))
      out%num(0:size - 1) = in%num(0:size - 1)
    end if

    call pop_sub()
  end subroutine distributed_copy

  subroutine distributed_end(this)
    type(distributed_t), intent(inout) :: this
    
    call push_sub('distributed.distributed_end')

    SAFE_DEALLOCATE_P(this%node)
    SAFE_DEALLOCATE_P(this%range)
    SAFE_DEALLOCATE_P(this%num)

    nullify(this%node, this%range, this%num)

    call pop_sub()
  end subroutine distributed_end

end module distributed_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
