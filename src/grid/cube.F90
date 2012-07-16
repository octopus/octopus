!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira, J. Alberdi
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

module cube_m
  use datasets_m
  use global_m
  use index_m
  use io_m
  use messages_m
  use mpi_m
  use fft_m
  use opencl_m
  use pfft_m
  use parser_m
  use profiling_m
  use simul_box_m

  implicit none
  private
  public ::             &
    cube_t,             &
    cube_init,          &
    cube_partition,     &
    cube_global2local,  &
    cube_end

  type cube_t
    logical :: parallel_in_domains !< will the cube be divided in domains?
    type(mpi_grp_t) :: mpi_grp     !< the mpi group describing parallelization in domains

    integer :: rs_n_global(1:3) !< the dimensions of the cube in real space
    integer :: fs_n_global(1:3) !< the dimensions of the cube in fourier space
    integer :: rs_n(1:3)        !< the dimensions of the local portion of the cube in real space
    integer :: fs_n(1:3)        !< the dimensions of the local portion of the cube in fourier space
    integer :: rs_istart(1:3)   !< where does the local portion of the cube start in real space
    integer :: fs_istart(1:3)   !< where does the local portion of the cube start in fourier space
    integer :: center(1:3)      !< the coordinates of the center of the cube

    integer, pointer :: np_local(:) !< Number of points in each partition
    integer, pointer :: xlocal(:)   !< where does each process start when gathering a function
    integer, pointer :: local(:,:)  !< local to global map used when gathering a function

    type(fft_t), pointer :: fft !< the fft object
  end type cube_t

contains

  ! ---------------------------------------------------------
  subroutine cube_init(cube, nn, sb, fft_type, fft_library, dont_optimize, nn_out, verbose)
    type(cube_t),      intent(out) :: cube
    integer,           intent(in)  :: nn(3)
    type(simul_box_t), intent(in)  :: sb
    integer, optional, intent(in)  :: fft_type  !< Is the cube going to be used to perform FFTs?
    integer, optional, intent(in)  :: fft_library !< What fft library to use
    logical, optional, intent(in)  :: dont_optimize !< if true, do not optimize grid for FFT
    integer, optional, intent(out) :: nn_out(3) !< What are the FFT dims?
                                                !! If optimized, may be different from input nn.
    logical, optional, intent(in)  :: verbose   !< Print info to the screen.

    integer :: mpi_comm, tmp_n(3), fft_type_, optimize_parity(3), default_lib, fft_library_
    integer :: effdim_fft
    logical :: optimize(3)

    PUSH_SUB(cube_init)

    ASSERT(all(nn(1:3) > 0))

    fft_type_ = optional_default(fft_type, FFT_NONE)

    effdim_fft = min (3, sb%dim)

    nullify(cube%fft)

    if (fft_type_ /= FFT_NONE) then

      if (present(fft_library)) then
        fft_library_ = fft_library
      else
        !%Variable FFTLibrary
        !%Type logical
        !%Section Mesh::FFTs
        !%Default fftw 
        !%Description
        !% (experimental) You can select the FFT library to use.
        !%Option fftw 1
        !% Uses FFTW3 library.
        !%Option pfft 2
        !% (experimental) Uses PFFT library, which has to be linked.
        !%Option clfft 3
        !% (experimental) Uses clAmdFft (GPU) library, which has to be linked.
        !%End
        default_lib = FFTLIB_FFTW
#ifdef HAVE_CLAMDFFT
        ! disabled by default since there are some problems for dim != 3
        ! if(opencl_is_enabled() .and. sb%dim == 3) default_lib = FFTLIB_CLAMD
#endif
        call parse_integer(datasets_check('FFTLibrary'), default_lib, fft_library_)
        if(optional_default(verbose, .false.)) call messages_print_var_option(stdout, 'FFTLibrary', fft_library_)
      end if
#ifndef HAVE_PFFT
      if (fft_library_ == FFTLIB_PFFT) then
        write(message(1),'(a)')'You have selected the PFFT for FFT, but it was not linked.'
        call messages_fatal(1)
      end if
#endif

      if (fft_library_ == FFTLIB_CLAMD) then
#ifndef HAVE_CLAMDFFT
        call messages_write('You have selected the OpenCL FFT, but Octopus was compiled', new_line = .true.)
        call messages_write('without clAmdFft (or OpenCL) support.')
        call messages_fatal()
#endif
        if(.not. opencl_is_enabled()) then
          call messages_write('You have selected the OpenCL FFT, but OpenCL is disabled.')
          call messages_fatal()
        end if
      end if

    else
      fft_library_ = FFTLIB_NONE
    end if

    cube%parallel_in_domains = fft_library_ == FFTLIB_PFFT
    if (fft_library_ == FFTLIB_NONE) then
      cube%rs_n_global = nn
      cube%fs_n_global = nn
      cube%rs_n = cube%rs_n_global
      cube%fs_n = cube%fs_n_global
      cube%rs_istart = 1
      cube%fs_istart = 1
      mpi_comm = -1
      if(present(nn_out)) nn_out(1:3) = nn(1:3)
    else
      SAFE_ALLOCATE(cube%fft)
      tmp_n = nn

      optimize(1:3) = .false.
      optimize_parity(1:3) = 0
      optimize(sb%periodic_dim+1:effdim_fft) = .true.
      optimize_parity(sb%periodic_dim+1:effdim_fft) = 1

      if(present(dont_optimize)) then
        if(dont_optimize) optimize = .false.
      endif
      call fft_init(cube%fft, tmp_n, sb%dim, fft_type_, fft_library_, optimize, optimize_parity, &
           mpi_comm=mpi_comm)
      if(present(nn_out)) nn_out(1:3) = tmp_n(1:3)

      call fft_get_dims(cube%fft, cube%rs_n_global, cube%fs_n_global, cube%rs_n, cube%fs_n, &
           cube%rs_istart, cube%fs_istart)
    end if
    cube%center(1:3) = cube%rs_n_global(1:3)/2 + 1

    call mpi_grp_init(cube%mpi_grp, mpi_comm)

    call cube_do_mapping(cube)

    if (cube%parallel_in_domains) call cube_partition_messages_debug(cube)

    POP_SUB(cube_init)
  end subroutine cube_init

  ! ---------------------------------------------------------
  subroutine cube_end(cube)
    type(cube_t), intent(inout) :: cube
    
    PUSH_SUB(cube_end)

    if(associated(cube%fft)) then
      call fft_end(cube%fft)
      SAFE_DEALLOCATE_P(cube%fft)
    end if

    SAFE_DEALLOCATE_P(cube%np_local)
    SAFE_DEALLOCATE_P(cube%xlocal)
    SAFE_DEALLOCATE_P(cube%local)

    POP_SUB(cube_end)
  end subroutine cube_end

  ! ---------------------------------------------------------
  !> True if global coordinates belong to this process. On output
  !! lxyz contains the local coordinates
  logical function cube_global2local(cube, ixyz, lxyz) result(is_here)
    type(cube_t), intent(in)  :: cube
    integer,      intent(in)  :: ixyz(3) !< global coordinates
    integer,      intent(out) :: lxyz(3) !< local coordinates

    lxyz(1) = ixyz(1) - cube%rs_istart(1) + 1
    lxyz(2) = ixyz(2) - cube%rs_istart(2) + 1
    lxyz(3) = ixyz(3) - cube%rs_istart(3) + 1
    is_here = lxyz(1) >= 1 .and. lxyz(1) <= cube%rs_n(1) .and. &
              lxyz(2) >= 1 .and. lxyz(2) <= cube%rs_n(2) .and. &
              lxyz(3) >= 1 .and. lxyz(3) <= cube%rs_n(3)
    
  end function cube_global2local

  ! ---------------------------------------------------------
  !> do the mapping between global and local points of the cube
  subroutine cube_do_mapping(cube)
    type(cube_t), intent(inout) :: cube
    
    integer :: tmp_local(6), position, process, ix, iy, iz, index
    integer, allocatable :: local_sizes(:)
    type(profile_t), save ::  prof_gt, prof_map

    PUSH_SUB(cube_do_mapping)

    !!BEGIN:gather the local information into a unique vector.
    !!do a gather in 3d of all the box, into a loop
    tmp_local(1) = cube%rs_istart(1)
    tmp_local(2) = cube%rs_istart(2)
    tmp_local(3) = cube%rs_istart(3)
    tmp_local(4) = cube%rs_n(1)
    tmp_local(5) = cube%rs_n(2)
    tmp_local(6) = cube%rs_n(3)

    if (cube%parallel_in_domains) then
      SAFE_ALLOCATE(local_sizes(6*cube%mpi_grp%size))
      call profiling_in(prof_gt,"CUBE_GAT")
#ifdef HAVE_MPI
      call MPI_Allgather(tmp_local, 6, MPI_INTEGER, local_sizes, 6, MPI_INTEGER,&
                         cube%mpi_grp%comm, mpi_err)
#endif
      call profiling_out(prof_gt)
    else
      SAFE_ALLOCATE(local_sizes(6))
      local_sizes = tmp_local
    end if

    call profiling_in(prof_map,"CUBE_MAP")

    SAFE_ALLOCATE(cube%xlocal(cube%mpi_grp%size))
    SAFE_ALLOCATE(cube%np_local(cube%mpi_grp%size))
    SAFE_ALLOCATE(cube%local(cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3), 3))

    index = 1
    do process = 1, cube%mpi_grp%size
      position = ((process-1)*6)+1
      if (position == 1) then
        cube%xlocal(1) = 1
        cube%np_local(1) = local_sizes(4)*local_sizes(5)*local_sizes(6)
      else
        ! calculate the begin index and size of each process
        cube%xlocal(process) = cube%xlocal(process-1) + cube%np_local(process-1)
        cube%np_local(process) = local_sizes(position+3)*local_sizes(position+4)*local_sizes(position+5)
      end if

      ! save the mapping between the global x,y,z and the global index
      ! and determine which partition the point belongs to
      do iz = local_sizes(position+2), local_sizes(position+2)+local_sizes(position+5)-1
        do iy = local_sizes(position+1), local_sizes(position+1)+local_sizes(position+4)-1
          do ix = local_sizes(position), local_sizes(position)+local_sizes(position+3)-1
            cube%local(index, 1) = ix
            cube%local(index, 2) = iy
            cube%local(index, 3) = iz
            index = index + 1
          end do
        end do
      end do
    end do

    call profiling_out(prof_map)

    SAFE_DEALLOCATE_A(local_sizes)

    POP_SUB(cube_do_mapping)
  end subroutine cube_do_mapping

  ! ---------------------------------------------------------
  subroutine cube_partition(cube, part)
    type(cube_t), intent(in)  :: cube
    integer,      intent(out) :: part(:,:,:)

    integer :: tmp_local(6), position, process, ix, iy, iz
    integer, allocatable :: local_sizes(:)

    PUSH_SUB(cube_partition)

    !!gather the local information into a unique vector.
    tmp_local(1) = cube%rs_istart(1)
    tmp_local(2) = cube%rs_istart(2)
    tmp_local(3) = cube%rs_istart(3)
    tmp_local(4) = cube%rs_n(1)
    tmp_local(5) = cube%rs_n(2)
    tmp_local(6) = cube%rs_n(3)

    if (cube%parallel_in_domains) then
      SAFE_ALLOCATE(local_sizes(6*cube%mpi_grp%size))
#ifdef HAVE_MPI
      call MPI_Allgather(tmp_local, 6, MPI_INTEGER, local_sizes, 6, MPI_INTEGER,&
                         cube%mpi_grp%comm, mpi_err)
#endif
    else
      SAFE_ALLOCATE(local_sizes(6))
      local_sizes = tmp_local
    end if

    do process = 1, cube%mpi_grp%size
      position = ((process-1)*6)+1

      do iz = local_sizes(position+2), local_sizes(position+2)+local_sizes(position+5)-1
        do iy = local_sizes(position+1), local_sizes(position+1)+local_sizes(position+4)-1
          do ix = local_sizes(position), local_sizes(position)+local_sizes(position+3)-1
            part(ix, iy, iz) = process
          end do
        end do
      end do
    end do

    POP_SUB(cube_partition)
  end subroutine cube_partition

  ! ---------------------------------------------------------
  subroutine cube_partition_messages_debug(cube)
    type(cube_t), intent(in) :: cube

    integer          :: nn, ii, jj, kk ! Counters.
    integer          :: npart
    integer          :: iunit          ! For debug output to files.
    character(len=3) :: filenum
    integer, allocatable :: part(:,:,:)

    PUSH_SUB(cube_partition_messages_debug)

    if(in_debug_mode) then
      SAFE_ALLOCATE(part(cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3)))
      call cube_partition(cube, part)
  
      if(mpi_grp_is_root(mpi_world)) then
        call io_mkdir('debug/cube_partition')
        npart = cube%mpi_grp%size
      
        ! Debug output. Write points of each partition in a different file.
        do nn = 1, npart

          write(filenum, '(i3.3)') nn

          iunit = io_open('debug/cube_partition/cube_partition.'//filenum, &
               action='write')
          do kk = 1, cube%rs_n_global(3)
            do jj = 1, cube%rs_n_global(2)
              do ii = 1, cube%rs_n_global(1)
                if(part(ii, jj, kk) .eq. nn) write(iunit, '(3i8)') ii, jj, kk
              end do
            end do
          end do
          call io_close(iunit)
        end do


      end if

      SAFE_DEALLOCATE_A(part)
    end if

#ifdef HAVE_MPI
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

    POP_SUB(cube_partition_messages_debug)
  end subroutine cube_partition_messages_debug

end module cube_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
