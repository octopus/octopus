!! Copyright (C) 2005-2006 Heiko Appel
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

#include "global.h"

module datasets_m
  use parser_m

  implicit none
  private

  public ::        &
    datasets_init, &
    datasets_end,  &
    datasets_check

  ! variables to treat multi datasets
  character(len=32), public, allocatable :: dataset_label(:)
  integer,           public, allocatable :: dataset_runmode(:), dataset_run_order(:)
  integer,           public              :: current_dataset
  integer,           public              :: no_datasets
  character(len=32), public              :: current_label
  character(len=32), public              :: tmpdir, inputdir, outputdir

contains

  ! ---------------------------------------------------------
  ! first we read the required information from the input file
  ! and prompt the user for possible errors in the input
  subroutine datasets_init(calc_mode, blk)
    integer,                 intent(in) :: calc_mode
    type(block_t), optional, intent(in) :: blk

    integer :: ids, n_rows, n_cols
    integer, allocatable :: order(:)
#ifdef HAVE_MPI
    integer :: mpi_err
#endif

    ! The block should be of the form
    ! %CalculationMode
    !   gs | unocc | td
    !   "run1" | "run2"
    !   1 | 2 | 3
    ! %
    ! where the second (the labels) and the third (the order) line are optional

    if(present(blk)) then
      no_datasets = parse_block_cols(blk, 0)
    else
      no_datasets = 1
    end if

    allocate(dataset_label(no_datasets))
    allocate(dataset_runmode(no_datasets))
    allocate(dataset_run_order(no_datasets))

    ! assign some reasonable defaults
    dataset_runmode(1) = calc_mode

    if(.not.present(blk)) then
      dataset_label(1)     = ''
      dataset_run_order(1) = 1
    else
      do ids = 1, no_datasets
        write(dataset_label(ids), '(a,i2.2,a)') "ds", ids, "_"
        dataset_run_order(ids) = ids
      end do

      ! get run modes
      do ids = 1, no_datasets
        call parse_block_integer(blk, 0, ids - 1, dataset_runmode(ids))
      end do

      ! the rest of the information is optional
      n_rows = parse_block_n(blk)

      ! first we get the labels
      if(n_rows >= 2) then
        n_cols = parse_block_cols(blk, 1)
        do ids = 1, n_cols
          call parse_block_string(blk, 1, ids - 1, dataset_label(ids))
        end do
      end if

      ! and now the order
      if(n_rows >= 3) then
        n_cols = parse_block_cols(blk, 1)
        do ids = 1, n_cols
          call parse_block_integer(blk, 2, ids - 1, dataset_run_order(ids))
        end do
      end if

      call parse_block_end(blk)
    end if

    ! consistency check: does the order match the run modes?
    allocate(order(no_datasets))
    order(:) = 1
    do ids = 1, no_datasets
      if(dataset_run_order(ids) >= 1 .and. dataset_run_order(ids) <= no_datasets) then
        order(dataset_run_order(ids)) = 0
      end if
    end do
    if(.not. all(order(:) == 0)) then
      write(0,'(a)') '*** Fatal Error (description follows)'
      write(0,'(a)') 'Ordering of the datasets seems to be wrong in block "%CalculationMode"'
#ifdef HAVE_MPI
      call MPI_Finalize(mpi_err)
#endif
      stop
    end if
    deallocate(order)

    current_label   = dataset_label(dataset_run_order(1))
    current_dataset = dataset_run_order(1)
  end subroutine datasets_init


  ! ---------------------------------------------------------
  subroutine datasets_end()

    deallocate(dataset_label, dataset_runmode, dataset_run_order)

  end subroutine datasets_end


  ! ---------------------------------------------------------
  character(len=64) function datasets_check(variable) result(var_name)
    character(len = * ), intent(in)  :: variable
    character(len = 64)              :: composite_name

    composite_name = trim(current_label)//trim(variable)

    if(parse_isdef(composite_name) .ne. 0) then
      ! composite name has been defined in the input file
      var_name = composite_name
    else
      ! could not find composite name in the input;
      ! will use bare variable name
      var_name = variable
    end if

  end function datasets_check

end module datasets_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
