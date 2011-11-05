!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!! $Id: cube_function.F90 8159 2011-08-09 01:45:03Z dstrubbe $

#include "global.h"

module cube_m
  use datasets_m
  use global_m
  use index_m
  use io_m
  use messages_m
  use mpi_m
  use fft_m
  use pfft_m
  use parser_m
  use par_vec_m
  use profiling_m
  use simul_box_m

  implicit none
  private
  public ::               &
    cube_t,               &
    cube_init,            &
#ifdef HAVE_PFFT
    cube_get_pfft_index,  &
#endif
    cube_end


  type cube_t
    integer :: n(1:3)      !< the global dimensions of the cube
    integer :: nx          ! = n(1)/2 + 1, first dimension of the FS array

    logical :: parallel_in_domains !< will the cube be divided in domains?
    type(mpi_grp_t) :: mpi_grp     !< the mpi group describing parallelization in domains

    integer :: rs_n(1:3)      !< the dimensions of the local portion of the cube in real space
    integer :: fs_n(1:3)      !< the dimensions of the local portion of the cube in fourier space
    integer :: rs_istart(1:3) !< where does the local portion of the cube start in real space
    integer :: fs_istart(1:3) !< where does the local portion of the cube start in fourier space
    integer :: np             !< the number of points in the local portion of the cube

    integer, pointer :: part(:,:,:)            !< point -> partition
    integer, pointer :: get_local_index(:,:,:) !< Mapping vector between local and global pfft indexing
    integer, pointer :: begin_indexes(:)       !< where does each process start
    integer, pointer :: block_sizes(:)         !< size of the block that is going to be used in the gatherv

    integer :: fft_library !< which FFT library to use. Options: NONE=0, FFTW=1, PFFT=2
    type(fft_t), pointer :: dfftw
    type(fft_t), pointer :: zfftw
#ifdef HAVE_PFFT
    type(pfft_t), pointer :: pfft
#endif
  end type cube_t


  integer, public, parameter :: &
       FFTLIB_NONE = 0, &
       FFTLIB_FFTW = 1, &
       FFTLIB_PFFT = 2

contains

  ! ---------------------------------------------------------
  subroutine cube_init(cube, n, sb, fft)
    type(cube_t),      intent(out) :: cube
    integer,           intent(in)  :: n(3)
    type(simul_box_t), intent(in)  :: sb
    logical, optional, intent(in)  :: fft !< Is the cube going to be used to perform FFTs?

    integer :: i1(1:3), i2(1:3)
    integer :: step, ip, mpi_comm
    logical :: fft_

    PUSH_SUB(cube_init)

    ASSERT(all(n(1:3) > 0))

    fft_ = .false.
    if (present(fft)) fft_ = fft

    cube%n = n

#ifdef HAVE_PFFT
    nullify(cube%pfft)
#endif
    nullify(cube%dfftw)
    nullify(cube%zfftw)

    if (fft_) then
      !%Variable FFTLibrary
      !%Type logical
      !%Section Hamiltonian::Poisson
      !%Default fftw 
      !%Description
      !% (experimental) You can select the FFT library to use.
      !%Option fftw 1
      !% Uses FFTW3 library.
      !%Option pfft 2
      !% (experimental) Uses PFFT library, which has to be linked.
      !%End
      call parse_integer(datasets_check('FFTLibrary'), FFTLIB_FFTW, cube%fft_library)
#ifndef HAVE_PFFT
      if (cube%fft_library == FFTLIB_PFFT) then
        write(message(1),'(a)')'You have selected the PFFT for FFT, but it is not compiled.'
        call messages_fatal(1)
      end if
#endif
    else
      cube%fft_library = FFTLIB_NONE
    end if


    select case (cube%fft_library)
    case (FFTLIB_NONE)
      cube%parallel_in_domains = .false.

    case (FFTLIB_FFTW)
      cube%parallel_in_domains = .false.

      SAFE_ALLOCATE(cube%dfftw)
      SAFE_ALLOCATE(cube%zfftw)

      call fft_init(cube%n, sb%dim, 0, cube%dfftw, optimize = .not.simul_box_is_periodic(sb))
      call fft_init(cube%n, sb%dim, 1, cube%zfftw, optimize = .not.simul_box_is_periodic(sb))
      cube%nx = cube%n(1)/2 + 1

      cube%np = cube%n(1)*cube%n(2)*cube%n(3)
      cube%rs_n = n
      cube%fs_n(1) = cube%nx
      cube%fs_n(2:3) = n(2:3)
      cube%rs_istart = 1
      cube%fs_istart = 1

    case (FFTLIB_PFFT)
#ifdef HAVE_PFFT
      cube%parallel_in_domains = .true.

      SAFE_ALLOCATE(cube%pfft)
      call pfft_init(cube%n, sb%dim, 1, cube%np, cube%rs_istart, cube%fs_istart, cube%rs_n, &
           cube%fs_n, mpi_comm, cube%pfft, optimize = .not.simul_box_is_periodic(sb))
      cube%nx = cube%n(1)/2 + 1
      call mpi_grp_init(cube%mpi_grp, mpi_comm)
#endif
    end select

    if (cube%parallel_in_domains) then
      call cube_do_mapping(cube)
      call cube_partition_messages_debug(cube)
    end if

  POP_SUB(cube_init)
end subroutine  cube_init

  ! ---------------------------------------------------------
  subroutine cube_end(cube)
    type(cube_t), intent(inout) :: cube
    
    PUSH_SUB(cube_end)

#ifdef HAVE_PFFT
    if(associated(cube%pfft)) then
      call pfft_end(cube%pfft)
      SAFE_DEALLOCATE_P(cube%pfft)
    end if
#endif  

    if (cube%parallel_in_domains) then
      SAFE_DEALLOCATE_P(cube%part)
      SAFE_DEALLOCATE_P(cube%begin_indexes)
      SAFE_DEALLOCATE_P(cube%block_sizes)
      SAFE_DEALLOCATE_P(cube%get_local_index)
    end if

    if(associated(cube%dfftw)) then
      call fft_end(cube%dfftw)
      SAFE_DEALLOCATE_P(cube%dfftw)
    end if
    if(associated(cube%zfftw)) then
      call fft_end(cube%zfftw)
      SAFE_DEALLOCATE_P(cube%zfftw)
    end if

    POP_SUB(cube_end)
  end subroutine  cube_end

  !> returns the local index for the PFFT library using the global x, y and z
  integer function cube_get_pfft_index(cube, ix, iy, iz) result(index)
    type(cube_t), intent(in) :: cube
    integer,      intent(in) :: ix !< x index
    integer,      intent(in) :: iy !< y index
    integer,      intent(in) :: iz !< z index
    
    integer :: normalized(1:3)

    !normalize to the min values
    normalized(1) = ix - cube%rs_istart(1) + 1 
    normalized(2) = iy - cube%rs_istart(2) + 1 
    normalized(3) = iz - cube%rs_istart(3) + 1 
      
    index = (normalized(1)-1) + &                              ! x component
         ((normalized(2)-1) * cube%rs_n(1)) + &                ! y component
         ((normalized(3)-1) * cube%rs_n(1) * cube%rs_n(2)) + 1 ! z component and the sum of 1
    
  end function cube_get_pfft_index

  ! ---------------------------------------------------------
  !> do the mapping between global and local points of PFFT arrays
  subroutine cube_do_mapping(cube)
    type(cube_t), intent(inout) :: cube
    
    integer :: tmp_local(6), position, process,ii, jj, kk, index
    integer, allocatable ::local_sizes(:)
    type(profile_t), save ::  prof_gt, prof_a

    PUSH_SUB(cube_do_mapping)
    call profiling_in(prof_gt,"PFFT_GAT")

    !!BEGIN:gather the local information into a unique vector.
    !!do a gather in 3d of all the box, into a loop
    tmp_local(1) = cube%rs_istart(1)
    tmp_local(2) = cube%rs_istart(2)
    tmp_local(3) = cube%rs_istart(3)
    tmp_local(4) = cube%rs_n(1)
    tmp_local(5) = cube%rs_n(2)
    tmp_local(6) = cube%rs_n(3) 

    SAFE_ALLOCATE(local_sizes(6*cube%mpi_grp%size))
#ifdef HAVE_MPI
    call MPI_Allgather(tmp_local, 6, MPI_INTEGER, local_sizes, 6, MPI_INTEGER,&
                       cube%mpi_grp%comm, mpi_err)
#endif
    call profiling_out(prof_gt)

    call profiling_in(prof_a,"PFFT_ALLOC")
    SAFE_ALLOCATE(cube%part(cube%n(1), cube%n(2), cube%n(3)))
    SAFE_ALLOCATE(cube%begin_indexes(cube%mpi_grp%size))
    SAFE_ALLOCATE(cube%block_sizes(cube%mpi_grp%size))
    SAFE_ALLOCATE(cube%get_local_index(cube%n(1),cube%n(2),cube%n(3)))

    do process = 1, cube%mpi_grp%size
      position = ((process-1)*6)+1
      if (position == 1) then
        cube%begin_indexes(1) = 1
        cube%block_sizes(1)  = local_sizes(4)*local_sizes(5)*local_sizes(6)
      else
        ! calculate the begin index and size of each process
        cube%begin_indexes(process) =  cube%begin_indexes(process-1) +  cube%block_sizes(process-1)
        cube%block_sizes(process) = local_sizes(position+3)*local_sizes(position+4)*local_sizes(position+5)
      end if

      !save the mapping between the global x,y,z and the local index and
      !get the partition
      index = 0
      do kk = local_sizes(position+2), local_sizes(position+2)+local_sizes(position+5)-1
        do jj = local_sizes(position+1), local_sizes(position+1)+local_sizes(position+4)-1
          do ii = local_sizes(position), local_sizes(position)+local_sizes(position+3)-1
            cube%part(ii,jj,kk) = process
            cube%get_local_index(ii,jj,kk) = index + cube%begin_indexes(process)
            index = index + 1
          end do
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(local_sizes)

    call profiling_out(prof_a)
    POP_SUB(cube_do_mapping)
  end subroutine cube_do_mapping

  subroutine cube_partition_messages_debug(cube)
    type(cube_t), intent(in)    :: cube

    integer          :: nn, ii, jj, kk ! Counters.
    integer          :: npart
    integer          :: iunit          ! For debug output to files.
    character(len=3) :: filenum

    PUSH_SUB(cube_partition_messages_debug)

    if(in_debug_mode .and. mpi_grp_is_root(mpi_world)) then

      call io_mkdir('debug/cube_partition')

      npart = cube%mpi_grp%size

      ! Debug output. Write points of each partition in a different file.
      do nn = 1, npart

        write(filenum, '(i3.3)') nn

        iunit = io_open('debug/cube_partition/cube_partition.'//filenum, &
             action='write')
        do kk = 1, cube%n(3)
          do jj = 1, cube%n(2)
            do ii = 1, cube%n(1)
              if(cube%part(ii, jj, kk) .eq. nn) write(iunit, '(3i8)') ii, jj, kk
            end do
          end do
        end do
        call io_close(iunit)
      end do
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
