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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module datasets_m
  use lib_oct_parser_m

  implicit none
  private

  public ::        &
    datasets_init, &
    datasets_end,  &
    check_inp

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
    integer, intent(in) :: calc_mode
    C_POINTER, optional, intent(in) :: blk
    integer :: i, n_rows, n_cols
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
      no_datasets = loct_parse_block_cols(blk, 0)
    else
      no_datasets = 1
    end if

    allocate(dataset_label(no_datasets), dataset_runmode(no_datasets))
    allocate(dataset_run_order(no_datasets))

    ! assign some reasonable defaults
    dataset_runmode(1) = calc_mode

    if(.not.present(blk)) then
      dataset_label(1)     = ''
      dataset_run_order(1) = 1
    else
      do i = 1, no_datasets
        write(dataset_label(i), '(a,i2.2,a)') "ds", i, "_"
        dataset_run_order(i) = i
      end do

      ! get run modes
      do i = 1, no_datasets
        call loct_parse_block_int(blk, 0, i-1, dataset_runmode(i))
      end do

      ! the rest of the information is optional
      n_rows = loct_parse_block_n(blk)

      ! first we get the labels
      if(n_rows >= 2) then
        n_cols = loct_parse_block_cols(blk, 1)
        do i = 1, n_cols
          call loct_parse_block_string(blk, 1, i-1, dataset_label(i))
        end do
      end if

      ! and now the order
      if(n_rows >= 3) then
        n_cols = loct_parse_block_cols(blk, 1)
        do i = 1, n_cols
          call loct_parse_block_int(blk, 2, i-1, dataset_run_order(i))
        end do
      end if

      call loct_parse_block_end(blk)
    end if

    ! consistency check: does the order match the run modes?
    allocate(order(no_datasets))
    order(:) = 1
    do i = 1, no_datasets
      if(dataset_run_order(i) >= 1 .and. dataset_run_order(i) <= no_datasets) then
        order(dataset_run_order(i)) = 0
      end if
    end do
    if(.not.all(order(:) == 0)) then
      write(6,'(a)') '*** Fatal Error (description follows)'
      write(6,'(a)') 'Ordering of the datasets seems to be wrong in block "%CalculationMode"'
#ifdef HAVE_MPI
      call MPI_Finalize(mpi_err)
#endif
      stop
    end if
    deallocate(order)

    current_label = trim(dataset_label(dataset_run_order(1)))
    current_dataset = dataset_run_order(1)
  end subroutine datasets_init


  ! ---------------------------------------------------------
  subroutine datasets_end()

    deallocate(dataset_label, dataset_runmode, dataset_run_order)

  end subroutine datasets_end


  ! ---------------------------------------------------------
  character(len=64) function check_inp(variable) result(var_name)
    character(len = * ), intent(in)  :: variable
    character(len = 64)              :: composite_name

    composite_name = trim(dataset_label(current_dataset))//trim(variable)

    if(loct_parse_isdef(composite_name).ne.0) then
      ! composite name has been defined in the input file
      var_name = composite_name
    else
      ! could not find composite name in the input;
      ! will use bare variable name
      var_name = variable
    end if

  end function check_inp

end module datasets_m
