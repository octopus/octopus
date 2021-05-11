!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio,
!! G. Bertsch, M. Oliveira, J. Alberdi-Rodriguez
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

module cube_oct_m
  use accel_oct_m
  use fft_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use nfft_oct_m
  use parser_oct_m
  use pfft_oct_m
  use pnfft_oct_m
  use profiling_oct_m
  use space_oct_m

  implicit none
  private
  public ::             &
    cube_t,             &
    dimensions_t,       &
    cube_init,          &
    cube_point_to_process, &
    cube_partition,     &
    cube_global2local,  &
    cube_getFFTLibrary, &
    cube_end

  type cube_t
    ! Components are public by default
    logical :: parallel_in_domains !< will the cube be divided in domains?
    type(mpi_grp_t) :: mpi_grp     !< the mpi group describing parallelization in domains

    integer :: rs_n_global(1:3) !< the dimensions of the cube in real space
    integer :: fs_n_global(1:3) !< the dimensions of the cube in fourier space
    integer :: rs_n(1:3)        !< the dimensions of the local portion of the cube in real space
    integer :: fs_n(1:3)        !< the dimensions of the local portion of the cube in fourier space
    integer :: rs_istart(1:3)   !< where does the local portion of the cube start in real space
    integer :: fs_istart(1:3)   !< where does the local portion of the cube start in fourier space
    integer :: center(1:3)      !< the coordinates of the center of the cube

    FLOAT, allocatable :: Lrs(:,:)  !< The real space coordinates vector: Lrs(i,{1,2,3})={x,y,z}(i)
    FLOAT, allocatable :: Lfs(:,:)  !< The fourier space coordinates vector: Lfs(i,{1,2,3})={kx,ky,kz}(i)

    integer, allocatable :: np_local(:)    !< Number of points in each partition
    integer, allocatable :: xlocal(:)      !< where does each process start when gathering a function
    integer, allocatable :: local(:,:)     !< local to global map used when gathering a function
    integer, allocatable :: np_local_fs(:) !< Number of points in each fs partition
    integer, allocatable :: xlocal_fs(:)   !< where does each process start when gathering a fs function
    integer, allocatable :: local_fs(:,:)  !< local to global map used when gathering a fs function


    type(fft_t), allocatable :: fft !< the fft object
    logical, private :: has_cube_mapping !< Saves if a mapping with the cube is needed.
                                         !! Until now, is needed with par_states (without par_domains) and PES.
  end type cube_t

  !> It is intended to be used within a vector.
  !!
  !! Each index of the vector corresponds to an MPI process.  A
  !! mapping between x,y,z index and process is saved, in a compact
  !! way.
  type dimensions_t
    private
    integer :: start_xyz(1:3) !< First index X, Y, Z, which this process has
    integer :: end_xyz(1:3)   !< Last  index X, Y, Z, which this process has
  end type dimensions_t

contains

  ! ---------------------------------------------------------
  subroutine cube_init(cube, nn, namespace, space, fft_type, fft_library, dont_optimize, nn_out, &
                       mpi_grp, need_partition, spacing, tp_enlarge, blocksize)
    type(cube_t),      intent(out) :: cube
    integer,           intent(in)  :: nn(3)
    type(namespace_t), intent(in)  :: namespace
    type(space_t),     intent(in)  :: space
    integer, optional, intent(in)  :: fft_type  !< Is the cube going to be used to perform FFTs?
    integer, optional, intent(in)  :: fft_library !< What fft library to use
    logical, optional, intent(in)  :: dont_optimize !< if true, do not optimize grid for FFT
    integer, optional, intent(out) :: nn_out(3) !< What are the FFT dims?
                                                !! If optimized, may be different from input nn.
    type(mpi_grp_t), optional, intent(in) :: mpi_grp !< The mpi group to be use for cube parallelization
    logical, optional, intent(in)  :: need_partition !< Should we calculate and store the cube partition?
    FLOAT, optional, intent(in)    :: spacing(3)
    FLOAT, optional, intent(in)    :: tp_enlarge(3)  !< Two point enlargement factor. Can be used with (p)nfft
                                                     !! enlarge the box by moving outward the cube first and last points
                                                     !! in each direction. The resulting boundaries are rescaled by tp_enlarge.
    integer, optional, intent(in)  :: blocksize !< just use a fixed block decomposition without caring about FFT library.
                                                !! See description for cube_set_blocksize.

    integer :: mpi_comm, tmp_n(3), fft_type_, optimize_parity(3), fft_library_
    integer :: effdim_fft
    logical :: optimize(3)
    type(mpi_grp_t) :: mpi_grp_
    FLOAT   :: tp_enlarge_(3)

    PUSH_SUB(cube_init)

    ASSERT(all(nn(1:3) > 0))

    fft_type_ = optional_default(fft_type, FFT_NONE)
    tp_enlarge_(:) = (/M_ONE, M_ONE, M_ONE/)
    if(present(tp_enlarge)) tp_enlarge_(:)=tp_enlarge(:)

    if (present(tp_enlarge)) then
      ASSERT(present(spacing))
    end if

    effdim_fft = min (3, space%dim)

    mpi_grp_ = mpi_world
    if (present(mpi_grp)) mpi_grp_ = mpi_grp

    if (fft_type_ /= FFT_NONE) then

      if (present(fft_library)) then
        fft_library_ = fft_library
      else
        fft_library_ = fft_default_lib
      end if

#ifndef HAVE_PFFT
      if (fft_library_ == FFTLIB_PFFT) then
        write(message(1),'(a)')'You have selected the PFFT for FFT, but it was not linked.'
        call messages_fatal(1)
      end if
#endif

    else
      fft_library_ = FFTLIB_NONE
    end if

    ! Note: later we set parallel_in_domains if blocksize is given, too
    cube%parallel_in_domains = (fft_library_ == FFTLIB_PFFT .or. fft_library_ == FFTLIB_PNFFT)
    if (present(blocksize)) then
      ASSERT(present(need_partition).and.need_partition)
      ASSERT(fft_library_ == FFTLIB_NONE)
      ! For all the different FFT libraries there are strange (?)
      ! rules about how the decomposition is chosen.  What we want
      ! (for libvdwxc) is a cube parallelized according to the simple
      ! but contrary rule "just do what I say".  Hence the blocksize
      ! parameter.  (Later to be expanded to allow 2D distributions.)
      cube%rs_n_global = nn
      cube%fs_n_global = nn ! not to be used
      cube%fs_n = cube%fs_n_global ! not to be used
      cube%fs_istart = 1 ! not to be used

      mpi_comm = mpi_grp_%comm
      cube%parallel_in_domains = (mpi_grp_%size > 1) ! XXX whether comm size > 1
      call cube_set_blocksize(cube%rs_n_global, blocksize, mpi_grp_%rank, cube%rs_n, cube%rs_istart)
    else if (fft_library_ == FFTLIB_NONE) then
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
      optimize(space%periodic_dim + 1:effdim_fft) = .true.
      optimize_parity(space%periodic_dim + 1:effdim_fft) = 1

      if(present(dont_optimize)) then
        if(dont_optimize) optimize = .false.
      end if

      if(present(tp_enlarge)) call cube_tp_fft_defaults(cube, fft_library_)

      call fft_init(cube%fft, tmp_n, space%dim, fft_type_, fft_library_, optimize, optimize_parity, &
        mpi_comm=mpi_comm, mpi_grp = mpi_grp_, use_aligned=.true.)
      if(present(nn_out)) nn_out(1:3) = tmp_n(1:3)

      call fft_get_dims(cube%fft, cube%rs_n_global, cube%fs_n_global, cube%rs_n, cube%fs_n, &
        cube%rs_istart, cube%fs_istart)

      if(present(tp_enlarge) .or. present(spacing)) then
        call cube_init_coords(cube, tp_enlarge_, spacing, fft_library_)
      end if

      if(fft_library_ == FFTLIB_NFFT .or. fft_library_ == FFTLIB_PNFFT) then
        call fft_init_stage1(cube%fft, namespace, cube%Lrs, cube%rs_n_global)
        !set local dimensions after stage1 - needed for PNFFT
        call fft_get_dims(cube%fft, cube%rs_n_global, cube%fs_n_global, cube%rs_n, cube%fs_n, &
             cube%rs_istart, cube%fs_istart)
      end if

    end if

    if(present(spacing) .and. .not. allocated(cube%Lrs)) then
      call cube_init_coords(cube, tp_enlarge_, spacing, fft_library_)
    end if


    
    cube%center(1:3) = cube%rs_n_global(1:3)/2 + 1
    

    call mpi_grp_init(cube%mpi_grp, mpi_comm)

    ! Initialize mapping only if needed
    if (present(need_partition) .and. cube%parallel_in_domains) then
      cube%has_cube_mapping = need_partition
    else
      cube%has_cube_mapping = .false.
    end if
    if (cube%has_cube_mapping) then
      call cube_do_mapping(cube, fs = fft_library_ == FFTLIB_PNFFT)
    end if

    if (cube%parallel_in_domains) call cube_partition_messages_debug(cube, namespace)

    POP_SUB(cube_init)
  end subroutine cube_init

  ! ---------------------------------------------------------
  subroutine cube_end(cube)
    type(cube_t), intent(inout) :: cube
    
    PUSH_SUB(cube_end)

    if (allocated(cube%fft)) then
      call fft_end(cube%fft)
      SAFE_DEALLOCATE_A(cube%fft)
    end if

    if (cube%has_cube_mapping) then
      SAFE_DEALLOCATE_A(cube%np_local)
      SAFE_DEALLOCATE_A(cube%xlocal)
      SAFE_DEALLOCATE_A(cube%local)

      SAFE_DEALLOCATE_A(cube%np_local_fs)
      SAFE_DEALLOCATE_A(cube%xlocal_fs)
      SAFE_DEALLOCATE_A(cube%local_fs)
    end if
    
    SAFE_DEALLOCATE_A(cube%Lrs)
    SAFE_DEALLOCATE_A(cube%Lfs)
    
    POP_SUB(cube_end)
  end subroutine cube_end


  ! ---------------------------------------------------------
  subroutine cube_tp_fft_defaults(cube, fft_library)
    type(cube_t), intent(inout) :: cube
    integer,      intent(in)    :: fft_library
    
    PUSH_SUB(cube_tp_fft_defaults)
    select case (fft_library)
      case (FFTLIB_NFFT)
        !Set NFFT defaults to values that gives good performance for two-point enlargement
        !These values are overridden by the NFFT options in the input file 
        cube%fft%nfft%set_defaults = .true.
        cube%fft%nfft%guru = .true.
        cube%fft%nfft%mm = 2 
        cube%fft%nfft%sigma = CNST(1.1)
        cube%fft%nfft%precompute = NFFT_PRE_PSI

      case (FFTLIB_PNFFT)
        cube%fft%pnfft%set_defaults = .true.
        cube%fft%pnfft%m = 2 
        cube%fft%pnfft%sigma = CNST(1.1)

      case default 
      !do nothing  
    end select
    
    POP_SUB(cube_tp_fft_defaults)
  end subroutine cube_tp_fft_defaults


  ! ---------------------------------------------------------
  subroutine cube_init_coords(cube, tp_enlarge, spacing, fft_library)
    type(cube_t), intent(inout) :: cube
    FLOAT,        intent(in) :: tp_enlarge(3)
    FLOAT,        intent(in) :: spacing(3)
    integer,      intent(in) :: fft_library
    
    FLOAT   :: temp
    integer :: ii, nn(3), maxn, idim

    PUSH_SUB(cube_init_coords)


    nn(1:3) = cube%fs_n_global(1:3) 

    maxn = maxval(nn) 
    SAFE_ALLOCATE(cube%Lrs(1:maxn, 1:3))
    cube%Lrs(:,:) = M_ZERO
    
    !! Real space coordinates 
    do idim = 1,3
      if (tp_enlarge(idim) > M_ONE ) then
        do ii = 2, nn(idim) - 1 
          cube%Lrs(ii, idim) = (ii - int(nn(idim)/2) -1) * spacing(idim)
        end do
        cube%Lrs(1,        idim) = (-int(nn(idim)/2)) * spacing(idim) * tp_enlarge(idim)
        cube%Lrs(nn(idim), idim) = (int(nn(idim)/2))  * spacing(idim) * tp_enlarge(idim)
      else
        do ii = 1, nn(idim) 
          cube%Lrs(ii, idim) = (ii - int(nn(idim)/2) -1) * spacing(idim)
        end do
      end if
    end do


    !! Fourier space coordinates 
    if(fft_library /= FFTLIB_NONE) then

      SAFE_ALLOCATE(cube%Lfs(1:maxn, 1:3))
      cube%Lfs(:,:) = M_ZERO

      do idim = 1,3
        temp = M_TWO * M_PI / (nn(idim) * spacing(idim))
!temp = M_PI / (nn * spacing(1))
        do ii = 1, nn(idim)
          if (fft_library == FFTLIB_NFFT .or. fft_library == FFTLIB_PNFFT ) then
            !The Fourier space is shrunk by the tp_enlarge factor
            !cube%Lfs(ii, 1:3) = (ii - nn/2 - 1)*temp/tp_enlarge
!HH NOTE:
!not sure this is the right general factor
            cube%Lfs(ii, idim) = (ii - nn(idim)/2 - 1)*temp/tp_enlarge(idim)
          else
            cube%Lfs(ii, idim) = pad_feq(ii,nn(idim), .true.) * temp
          end if
        end do
      end do
    end if
    
    POP_SUB(cube_init_coords)
  end subroutine cube_init_coords


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
  !> Returns the FFT library of the cube.
  !! Possible values are FFTLIB_NONE, FFTLIB_FFTW, FFTLIB_PFFT 
  !! FFTLIB_ACCEL, FFTLIB_NFFT and FFTLIB_PNFFT 
  !! (defined in fft.F90)
  integer function cube_getFFTLibrary(cube) result(fft_library)
    type(cube_t), intent(in)  :: cube

    if (allocated(cube%fft)) then
       fft_library = cube%fft%library
    else
       fft_library = FFTLIB_NONE
    end if
  end function cube_getFFTLibrary

  ! ---------------------------------------------------------
  !> do the mapping between global and local points of the cube
  subroutine cube_do_mapping(cube, fs)
    type(cube_t), intent(inout) :: cube
    logical,      intent(in)    :: fs !< fill the fs mapping too
    
    integer :: tmp_local(6), position, process, ix, iy, iz, index
    integer, allocatable :: local_sizes(:)
    type(profile_t), save ::  prof_gt, prof_map
    integer(8) :: number_points

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
      SAFE_ALLOCATE(local_sizes(1:6*cube%mpi_grp%size))
      call profiling_in(prof_gt,"CUBE_GAT")
#ifdef HAVE_MPI
      call MPI_Allgather(tmp_local, 6, MPI_INTEGER, local_sizes, 6, MPI_INTEGER,&
                         cube%mpi_grp%comm, mpi_err)
#endif
      call profiling_out(prof_gt)
    else
      SAFE_ALLOCATE(local_sizes(1:6))
      local_sizes = tmp_local
    end if

    call profiling_in(prof_map,"CUBE_MAP")

    SAFE_ALLOCATE(cube%xlocal(1:cube%mpi_grp%size))
    SAFE_ALLOCATE(cube%np_local(1:cube%mpi_grp%size))
    ! make sure we do not run into integer overflow here
    number_points = cube%rs_n_global(1) * cube%rs_n_global(2)
    number_points = number_points * cube%rs_n_global(3)
    if (number_points >= HUGE(0)) then
      message(1) = "Error: too many points for the normal cube. Please try to use a distributed FFT."
      call messages_fatal(1)
    end if
    SAFE_ALLOCATE(cube%local(1:cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3), 3))

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
    
    if(optional_default(fs,.false.)) then

      tmp_local(1) = cube%fs_istart(1)
      tmp_local(2) = cube%fs_istart(2)
      tmp_local(3) = cube%fs_istart(3)
      tmp_local(4) = cube%fs_n(1)
      tmp_local(5) = cube%fs_n(2)
      tmp_local(6) = cube%fs_n(3)

      local_sizes = 0
      if (cube%parallel_in_domains) then
        call profiling_in(prof_gt,"CUBE_GAT_FS")
  #ifdef HAVE_MPI
        call MPI_Allgather(tmp_local, 6, MPI_INTEGER, local_sizes, 6, MPI_INTEGER,&
                           cube%mpi_grp%comm, mpi_err)
  #endif
        call profiling_out(prof_gt)
      else
        local_sizes = tmp_local
      end if

      call profiling_in(prof_map,"CUBE_MAP_FS")

      SAFE_ALLOCATE(cube%xlocal_fs(1:cube%mpi_grp%size))
      SAFE_ALLOCATE(cube%np_local_fs(1:cube%mpi_grp%size))
      ! make sure we do not run into integer overflow here
      number_points = cube%fs_n_global(1) * cube%fs_n_global(2)
      number_points = number_points * cube%fs_n_global(3)
      if (number_points >= HUGE(0)) then
        message(1) = "Error: too many points for the normal cube. Please try to use a distributed FFT."
        call messages_fatal(1)
      end if
      SAFE_ALLOCATE(cube%local_fs(1:cube%fs_n_global(1)*cube%fs_n_global(2)*cube%fs_n_global(3), 3))

      index = 1
      do process = 1, cube%mpi_grp%size
        position = ((process-1)*6)+1
        if (position == 1) then
          cube%xlocal_fs(1) = 1
          cube%np_local_fs(1) = local_sizes(4)*local_sizes(5)*local_sizes(6)
        else
          ! calculate the begin index and size of each process
          cube%xlocal_fs(process) = cube%xlocal_fs(process-1) + cube%np_local_fs(process-1)
          cube%np_local_fs(process) = local_sizes(position+3)*local_sizes(position+4)*local_sizes(position+5)
        end if

        ! save the mapping between the global x,y,z and the global index
        ! and determine which partition the point belongs to
        do iz = local_sizes(position+2), local_sizes(position+2)+local_sizes(position+5)-1
          do iy = local_sizes(position+1), local_sizes(position+1)+local_sizes(position+4)-1
            do ix = local_sizes(position), local_sizes(position)+local_sizes(position+3)-1
              cube%local_fs(index, 1) = ix
              cube%local_fs(index, 2) = iy
              cube%local_fs(index, 3) = iz
              index = index + 1
            end do
          end do
        end do
      end do
      
      call profiling_out(prof_map)

    end if

    

    SAFE_DEALLOCATE_A(local_sizes)

    POP_SUB(cube_do_mapping)
  end subroutine cube_do_mapping

  !!> Given a x, y, z point of the cube, it returns the corresponding process
  !!
  !! last_found is used to speed-up the search
  integer pure function cube_point_to_process(xyz, part) result(process)
    integer, intent(in)   :: xyz(1:3)
    type(dimensions_t), intent(in) :: part(:)
    
    integer :: proc
    logical :: found

    ! No PUSH/POP because it is a PURE function

    found = .false.
    do proc = 1, mpi_world%size
      !Compare XYZ index
      if ( all(xyz >= part(proc)%start_xyz) .and. all(xyz < part(proc)%end_xyz) ) then
        process = proc
        found = .true.
        exit
      end if
    end do
      
    ! An error message should be raised, if this point is reached
    if (.not. found) then
      process = -1
    end if

  end function cube_point_to_process

  ! Sets a 1D decomposition with fixed-size blocks over the last (least-contiguous) axis.
  ! Each core will have <blocksize> slices except the last one which will typically have
  ! less.  (In some cases, there can be multiple trailing cores without any slices.)
  subroutine cube_set_blocksize(rs_n_global, blocksize, rank, rs_n, rs_istart)
    integer,  intent(in) :: rs_n_global(1:3)
    integer,  intent(in) :: blocksize
    integer,  intent(in) :: rank
    integer, intent(out) :: rs_n(1:3)
    integer, intent(out) :: rs_istart(1:3)

    integer :: imin, imax

    rs_n = rs_n_global
    rs_istart = 1

    imin = min(blocksize * rank, rs_n_global(3))
    imax = min(imin + blocksize, rs_n_global(3))
    rs_istart(3) = 1 + imin
    rs_n(3) = imax - imin
  end subroutine cube_set_blocksize

  ! ---------------------------------------------------------
  subroutine cube_partition(cube, part)
    type(cube_t),       intent(in)  :: cube
    type(dimensions_t), intent(out) :: part(:)

    integer :: tmp_local(6), position, process
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
      SAFE_ALLOCATE(local_sizes(1:6*cube%mpi_grp%size))
#ifdef HAVE_MPI
      call MPI_Allgather(tmp_local, 6, MPI_INTEGER, local_sizes, 6, MPI_INTEGER,&
                         cube%mpi_grp%comm, mpi_err)
#endif
    else
      SAFE_ALLOCATE(local_sizes(1:6))
      local_sizes = tmp_local
    end if

    do process = 1, cube%mpi_grp%size
      position = ((process-1)*6)+1
      
      part(process)%start_xyz(1) = local_sizes(position)
      part(process)%start_xyz(2) = local_sizes(position+1) 
      part(process)%start_xyz(3) = local_sizes(position+2) 
      part(process)%end_xyz(1)   = local_sizes(position)+local_sizes(position+3)
      part(process)%end_xyz(2)   = local_sizes(position+1)+local_sizes(position+4)
      part(process)%end_xyz(3)   = local_sizes(position+2)+local_sizes(position+5)
      
    end do

    POP_SUB(cube_partition)
  end subroutine cube_partition

  ! ---------------------------------------------------------
  subroutine cube_partition_messages_debug(cube, namespace)
    type(cube_t),      intent(in) :: cube
    type(namespace_t), intent(in) :: namespace

    integer          :: nn, ii, jj, kk ! Counters.
    integer          :: ixyz(3)        ! Current value of xyz 
    integer          :: npart
    integer          :: iunit          ! For debug output to files.
    character(len=3) :: filenum
    type(dimensions_t), allocatable :: part(:)
    
    PUSH_SUB(cube_partition_messages_debug)

    if(debug%info) then
      SAFE_ALLOCATE(part(1:cube%mpi_grp%size))
      call cube_partition(cube, part)
  
      if(mpi_grp_is_root(mpi_world)) then
        call io_mkdir('debug/cube_partition', namespace)
        npart = cube%mpi_grp%size
      
        ! Debug output. Write points of each partition in a different file.
        do nn = 1, npart

          write(filenum, '(i3.3)') nn

          iunit = io_open('debug/cube_partition/cube_partition.'//filenum, &
            namespace, action='write')
          do kk = 1, cube%rs_n_global(3)
            do jj = 1, cube%rs_n_global(2)
              do ii = 1, cube%rs_n_global(1)
                ixyz(1) = ii; ixyz(2) = jj; ixyz(3) = kk
                if(cube_point_to_process(ixyz, part) == nn) then
                  write(iunit, '(3i8)') ii, jj, kk
                end if
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

end module cube_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
