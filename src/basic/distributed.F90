!! Copyright (C) 2008-2016 X. Andrade
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

module distributed_oct_m
  use global_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use profiling_oct_m

  implicit none

  private

  public                            &
       distributed_t,               &
       distributed_nullify,         &
       distributed_init,            &
       distributed_copy,            &
       distributed_end,             &
       distributed_allgather
  

  type distributed_t
    integer              :: start
    integer              :: end
    integer              :: nlocal
    integer              :: nglobal
    logical              :: parallel
    integer, allocatable :: node(:)
    integer, allocatable :: range(:, :)
    integer, allocatable :: num(:)
    type(mpi_grp_t)      :: mpi_grp
  end type distributed_t
  
contains

  ! ---------------------------------------------------------
  subroutine distributed_nullify(this, total)
    type(distributed_t), intent(out) :: this
    integer, optional,   intent(in)  :: total

    PUSH_SUB(distributed_nullify)

    this%start           = 1

    if(present(total)) then
      this%end             = total
      this%nlocal          = total
      this%nglobal         = total
    end if

    this%parallel        = .false.

    call mpi_grp_init(this%mpi_grp, -1)

    POP_SUB(distributed_nullify)
  end subroutine distributed_nullify
  

  ! ---------------------------------------------------------
  subroutine distributed_init(this, total, comm, tag, scalapack_compat)
    type(distributed_t),        intent(out) :: this
    integer,                    intent(in)  :: total
    integer,                    intent(in)  :: comm
    character(len=*), optional, intent(in)  :: tag
    logical,          optional, intent(in)  :: scalapack_compat

    integer :: kk

    PUSH_SUB(distributed_init)
    
    this%nglobal = total

    call mpi_grp_init(this%mpi_grp, comm)
    if(this%mpi_grp%size == 1 .or. this%nglobal == 1) then
      
      SAFE_ALLOCATE(this%node(1:total))
      SAFE_ALLOCATE(this%range(1:2, 0:this%mpi_grp%size - 1))
      SAFE_ALLOCATE(this%num(0:this%mpi_grp%size - 1))
      
      ! Defaults.
      this%node(1:total)   = 0
      this%start           = 1
      this%end             = total
      this%nlocal          = total
      this%nglobal         = total
      this%parallel        = .false.
      this%range(1, 0)     = 1
      this%range(2, 0)     = total
      this%num(0)          = total
      
      call mpi_grp_init(this%mpi_grp, -1)
      
    else

      this%parallel = .true.

      SAFE_ALLOCATE(this%range(1:2, 0:this%mpi_grp%size - 1))
      SAFE_ALLOCATE(this%num(0:this%mpi_grp%size - 1))
      SAFE_ALLOCATE(this%node(1:this%nglobal))

      call multicomm_divide_range(this%nglobal, this%mpi_grp%size, this%range(1, :), this%range(2, :), &
        lsize = this%num, scalapack_compat = scalapack_compat)

      if(present(tag)) then
        message(1) = 'Info: Parallelization in ' // trim(tag)
        call messages_info(1)
      end if

      do kk = 1, this%mpi_grp%size

        if(present(tag)) then
          write(message(1),'(a,i4,a,i6,a)') 'Info: Node in group ', kk - 1, &
            ' will manage ', this%num(kk - 1), ' '//trim(tag)
          if(this%num(kk - 1) > 0) then
            write(message(1),'(a,a,i6,a,i6)') trim(message(1)), ':', this%range(1, kk - 1), " - ", this%range(2, kk - 1)
          end if
          call messages_info(1)
        end if
        
        if(this%mpi_grp%rank  ==  kk - 1) then
          this%start  = this%range(1, kk - 1)
          this%end    = this%range(2, kk - 1)
          this%nlocal = this%num(kk - 1)
        end if
        
        this%node(this%range(1, kk - 1):this%range(2, kk - 1)) = kk - 1
        
      end do
      
    end if
    
    POP_SUB(distributed_init)
  end subroutine distributed_init


  ! ---------------------------------------------------------
  subroutine distributed_copy(in, out)
    type(distributed_t), intent(in)  :: in
    type(distributed_t), intent(out) :: out

    integer :: size

    PUSH_SUB(distributed_copy)

    out%start    = in%start
    out%end      = in%end
    out%nlocal   = in%nlocal
    out%nglobal  = in%nglobal
    out%parallel = in%parallel

    size = in%mpi_grp%size

    call mpi_grp_init(out%mpi_grp, in%mpi_grp%comm)
    
    if(allocated(in%node)) then
      SAFE_ALLOCATE(out%node(1:in%nglobal))
      out%node(1:in%nglobal) = in%node(1:in%nglobal)
    end if

    if(allocated(in%range)) then
      SAFE_ALLOCATE(out%range(1:2, 0:size - 1))
      out%range(1:2, 0:size - 1) = in%range(1:2, 0:size - 1)
    end if

    if(allocated(in%num)) then
      SAFE_ALLOCATE(out%num(0:size - 1))
      out%num(0:size - 1) = in%num(0:size - 1)
    end if

    POP_SUB(distributed_copy)
  end subroutine distributed_copy


  ! ---------------------------------------------------------
  subroutine distributed_end(this)
    type(distributed_t), intent(inout) :: this
    
    PUSH_SUB(distributed_end)

    SAFE_DEALLOCATE_A(this%node)
    SAFE_DEALLOCATE_A(this%range)
    SAFE_DEALLOCATE_A(this%num)

    POP_SUB(distributed_end)
  end subroutine distributed_end

  ! --------------------------------------------------------

  subroutine distributed_allgather(this, aa)
    type(distributed_t), intent(in)    :: this
    FLOAT,               intent(inout) :: aa(:)

    integer, allocatable :: displs(:)

    if(.not. this%parallel) return
    
    PUSH_SUB(distributed_allgather)

    SAFE_ALLOCATE(displs(0:this%mpi_grp%size - 1))

    displs(0:this%mpi_grp%size - 1) = this%range(1, 0:this%mpi_grp%size - 1) - 1

#ifdef HAVE_MPI    
    call MPI_Allgatherv(MPI_IN_PLACE, this%nlocal, MPI_FLOAT, &
      aa(1), this%num(0), displs(0), MPI_FLOAT, this%mpi_grp%comm, mpi_err)
#endif
    
    SAFE_DEALLOCATE_A(displs)
    
    POP_SUB(distributed_allgather)
  end subroutine distributed_allgather

end module distributed_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
