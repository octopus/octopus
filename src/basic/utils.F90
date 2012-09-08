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
!! $Id: util.F90 3013 2007-06-21 15:53:17Z xavier $

#include "global.h"

!> This module is intended to contain simple general-purpose utility functions
!! and procedures.

module utils_m
  use global_m
  use io_m
  use loct_m
  use messages_m
  use mpi_m
  use unit_m
  use unit_system_m
  use string_m

  implicit none

  private
  public ::                       &
    get_divisors,                 &
    index2axis,                   &
    output_tensor,                &
    output_dipole,                &
    print_header,                 &
    leading_dimension_is_known,   &
    lead_dim

  interface leading_dimension_is_known
    module procedure dleading_dimension_is_known, zleading_dimension_is_known
  end interface leading_dimension_is_known

  interface lead_dim
    module procedure dlead_dim, zlead_dim
  end interface lead_dim

contains

  ! ---------------------------------------------------------
  subroutine get_divisors(nn, n_divisors, divisors)
    integer, intent(in)    :: nn
    integer, intent(inout) :: n_divisors
    integer, intent(out)   :: divisors(:)

    integer :: ii, max_d

    PUSH_SUB(get_divisors)

    ASSERT(n_divisors > 1)
    max_d = n_divisors

    n_divisors = 1
    divisors(n_divisors) = 1
    do ii = 2, nn / 2
      if(mod(nn, ii)==0) then
        n_divisors = n_divisors + 1

        if(n_divisors > max_d - 1) then
          message(1) = "Internal error in get_divisors. Please increase n_divisors"
          call messages_fatal(1)
        end if

        divisors(n_divisors) = ii
      end if
    end do
    n_divisors = n_divisors + 1
    divisors(n_divisors) = nn

    POP_SUB(get_divisors)
  end subroutine get_divisors


  ! ---------------------------------------------------------
  character function index2axis(idir) result(ch)
    integer, intent(in) :: idir
    
    PUSH_SUB(index2axis)

    select case(idir)
      case(1)
        ch = 'x'
      case(2)
        ch = 'y'
      case(3)
        ch = 'z'
      case(4)
        ch = 'w'
      case default
        write(ch,'(i1)') idir
    end select

    POP_SUB(index2axis)
  end function index2axis


  ! ---------------------------------------------------------
  subroutine output_tensor(iunit, tensor, ndim, unit, write_average)
    integer,           intent(in) :: iunit
    FLOAT,             intent(in) :: tensor(:,:)
    integer,           intent(in) :: ndim
    type(unit_t),      intent(in) :: unit
    logical, optional, intent(in) :: write_average
    
    FLOAT :: trace
    integer :: jj, kk
    logical :: write_average_

    PUSH_SUB(output_tensor)

    write_average_ = .true.
    if(present(write_average)) write_average_ = write_average

    trace = M_z0
    do jj = 1, ndim
      write(iunit, '(3f20.6)') (units_from_atomic(unit, tensor(jj, kk)), kk=1,ndim)
      trace = trace + tensor(jj, jj)
    end do

    trace = units_from_atomic(unit, trace/TOFLOAT(ndim))

    if(write_average_) write(iunit, '(a, f20.6)')  'Isotropic average', trace

    POP_SUB(output_tensor)
  end subroutine output_tensor


  ! ---------------------------------------------------------
  subroutine output_dipole(iunit, dipole, ndim) 
    integer, intent(in) :: iunit
    FLOAT,   intent(in) :: dipole(:)
    integer, intent(in) :: ndim
    
    integer :: idir

    PUSH_SUB(output_dipole)

    write(iunit, '(a,a20,a17)') 'Dipole:', '[' // trim(units_abbrev(units_out%length)) // ']', &
          '[' // trim(units_abbrev(unit_debye)) // ']'
    do idir = 1, ndim
      write(iunit, '(6x,3a,es14.5,3x,2es14.5)') '<', index2axis(idir), '> = ', &
        units_from_atomic(units_out%length, dipole(idir)), units_from_atomic(unit_debye, dipole(idir))
    end do

    POP_SUB(output_dipole)
  end subroutine output_dipole

  !> This subroutine prints the logo followed by information about 
  !! the compilation and the system. It also prints the start time 
  !! of the execution
  ! ---------------------------------------------------------
  subroutine print_header()
    
    character(len=256) :: sys_name
    
    ! Let us print our logo
    if(mpi_grp_is_root(mpi_world)) then
      call io_dump_file(stdout, trim(trim(conf%share) // '/logo'))
    end if

    ! Let us print the version
    message(1) = ""
    message(2) = str_center("Running octopus", 70)
    message(3) = ""
    call messages_info(3)

    message(1) = &
         "Version                : " // trim(conf%version)
    message(2) = &
         "Revision               : "// trim(conf%latest_svn)
    message(3) = &
         "Build time             : "// trim(conf%build_time)
    call messages_info(3)

    write(message(1), '(a, i1)') &
         'Configuration options  : max-dim=', MAX_DIM
!!$
#ifdef HAVE_OPENMP
    message(1) = trim(message(1))//' openmp'
#endif
#ifdef HAVE_MPI
    message(1) = trim(message(1))//' mpi'
#endif
#ifdef HAVE_OPENCL
    message(1) = trim(message(1))//' opencl'
#endif
#ifdef HAVE_CLAMDFFT
    message(1) = trim(message(1))//' clamdfft'
#endif
#ifdef HAVE_M128D
    message(1) = trim(message(1))//' sse2'
#endif
#ifdef HAVE_M256D
    message(1) = trim(message(1))//' avx'
#endif
#ifdef HAVE_BLUE_GENE
    message(1) = trim(message(1))//' bluegene'
#endif

    message(2) = &
         'Optional libraries     :'
#ifdef HAVE_MPI2
    message(2) = trim(message(2))//' mpi2'
#endif
#ifdef HAVE_NETCDF
    message(2) = trim(message(2))//' netcdf'
#endif
#ifdef HAVE_METIS
    message(2) = trim(message(2))//' metis'
#endif
#ifdef HAVE_GDLIB
    message(2) = trim(message(2))//' gdlib'
#endif
#ifdef HAVE_PAPI
    message(2) = trim(message(2))//' papi'
#endif
#ifdef HAVE_SPARSKIT
    message(2) = trim(message(2))//' sparskit'
#endif
#ifdef HAVE_ETSF_IO
    message(2) = trim(message(2))//' etsf_io'
#endif
#ifdef HAVE_BERKELEYGW
    message(2) = trim(message(2))//' berkeleygw'
#endif
#ifdef HAVE_PFFT
    message(2) = trim(message(2))//' pfft'
#endif
#ifdef HAVE_NFFT
    message(2) = trim(message(2))//' nfft'
#endif
#ifdef HAVE_SCALAPACK
    message(2) = trim(message(2))//' scalapack'
#endif
#ifdef HAVE_ARPACK
    message(2) = trim(message(2))//' arpack'
#endif
#ifdef HAVE_LIBFM
    message(2) = trim(message(2))//' libfm'
#endif

    message(3) = &
         'Architecture           : '// TOSTRING(OCT_ARCH)
    call messages_info(3)

    message(1) = &
         "C compiler             : "//trim(conf%cc)
    message(2) = &
         "C compiler flags       : "//trim(conf%cflags)
    message(3) = &
         "Fortran compiler       : "//trim(conf%fc)
    message(4) = &
         "Fortran compiler flags : "//trim(conf%fcflags)
    call messages_info(4)

    message(1) = ""
    call messages_info(1)

    ! Let us print where we are running
    call loct_sysname(sys_name)
    write(message(1), '(a)') str_center("The octopus is swimming in " // trim(sys_name), 70)
    message(2) = ""
    call messages_info(2)

#if defined(HAVE_MPI)
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

    call print_date("Calculation started on ")
  end subroutine print_header

  ! ---------------------------------------------------------
  
  logical function dleading_dimension_is_known(array) result(known)
    FLOAT, intent(in) :: array(:, :)
    
    known = .true.

#if defined(HAVE_FORTRAN_LOC) && defined(HAVE_FC_SIZEOF)
    if(ubound(array, dim = 2) > 1) then
      known = ubound(array, dim = 1) == (loc(array(1, 2)) - loc(array(1, 1)))/sizeof(array(1, 1))
    end if
#endif

  end function dleading_dimension_is_known


  ! ---------------------------------------------------------
  
  logical function zleading_dimension_is_known(array) result(known)
    CMPLX, intent(in) :: array(:, :)
    
    known = .true.

#if defined(HAVE_FORTRAN_LOC) && defined(HAVE_FC_SIZEOF)
    if(ubound(array, dim = 2) > 1) then
      known = ubound(array, dim = 1) == (loc(array(1, 2)) - loc(array(1, 1)))/sizeof(array(1, 1))
    end if
#endif
    
  end function zleading_dimension_is_known

  ! ---------------------------------------------------------

  integer function dlead_dim(array) result(lead_dim)
    FLOAT, intent(in) :: array(:, :)
    
    ASSERT(leading_dimension_is_known(array))
    
    lead_dim = ubound(array, dim = 1)
  end function dlead_dim

  ! ---------------------------------------------------------

  integer function zlead_dim(array) result(lead_dim)
    CMPLX, intent(in) :: array(:, :)
    
    ASSERT(leading_dimension_is_known(array))
    
    lead_dim = ubound(array, dim = 1)
  end function zlead_dim

end module utils_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
