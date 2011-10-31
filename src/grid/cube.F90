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
  use messages_m
  use mpi_m
  use fft_m
  use pfft_m
  use parser_m
  use par_vec_m
  use pfft_m
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
    integer :: n_part(1:3) !< the dimensions of the local portion of the cube
    integer :: nx          ! = n(1)/2 + 1, first dimension of the FS array

    integer :: istart(1:3) !< where does the local portion of the cube start

    logical :: parallel_in_domains !< will the cube be divided in domains?
    type(mpi_grp_t) :: mpi_grp     !< the mpi group describing parallelization in domains

    integer :: offset(1:3) !< 

    integer :: fft_library          !< which FFT library has to be used. Options FFTW3=0, PFFT=1
    type(fft_t), pointer :: dfftw
    type(fft_t), pointer :: zfftw
#ifdef HAVE_PFFT
    type(pfft_t), pointer :: pfft
#endif
  end type cube_t


  integer, public, parameter :: &
       FFTW3_LIB  = 0, &
       PFFT_LIB   = 1

  integer, public, parameter :: &
       MCM_POINT = 1, &
       MCM_COUNT = 2


contains

  ! ---------------------------------------------------------
  subroutine cube_init(cube, n, sb, idx, np_global)
    type(cube_t),      intent(out) :: cube
    integer,           intent(in)  :: n(3)
    type(simul_box_t), intent(in)  :: sb
    type(index_t),     intent(in)  :: idx
    integer,           intent(in)  :: np_global

    integer :: i1(1:3), i2(1:3)
    integer :: step, ip

    PUSH_SUB(cube_init)

    ASSERT(all(n(1:3) > 0))

    cube%n = n
    cube%offset = 0
    
#ifdef HAVE_PFFT
    nullify(cube%pfft)
#endif
    nullify(cube%dfftw)
    nullify(cube%zfftw)

    !%Variable FFTLibrary
    !%Type logical
    !%Section Hamiltonian::Poisson
    !%Default fftw 
    !%Description
    !% (experimental) You can select the FFT library to use.
    !%Option fftw 0
    !% Uses FFTW3 library.
    !%Option pfft 1
    !% (experimental) Uses PFFT library, which has to be linked.
    !%End
    call parse_integer(datasets_check('FFTLibrary'), FFTW3_LIB, cube%fft_library)
#ifndef HAVE_PFFT
    if (cube%fft_library == PFFT_LIB) then
      write(message(1),'(a)')'You have selected the PFFT for FFT, but it is not compiled.'
      call messages_fatal(1)
    end if
#endif

    if (cube%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
      SAFE_ALLOCATE(cube%pfft)
      call pfft_init(cube%n, sb%dim, 1, cube%pfft, optimize = .not.simul_box_is_periodic(sb))
      cube%nx = cube%n(1)/2 + 1
#endif
    else
      SAFE_ALLOCATE(cube%dfftw)
      SAFE_ALLOCATE(cube%zfftw)

      call fft_init(cube%n, sb%dim, 0, cube%dfftw, optimize = .not.simul_box_is_periodic(sb))
      call fft_init(cube%n, sb%dim, 1, cube%zfftw, optimize = .not.simul_box_is_periodic(sb))
      cube%nx = cube%n(1)/2 + 1
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
  

#ifdef HAVE_PFFT
  !> returns the local index for the PFFT library using the global x, y and z
  integer function cube_get_pfft_index(pfft, ix, iy, iz) result(index)
    type(pfft_t), intent(in) :: pfft
    integer,      intent(in) :: ix !< x index
    integer,      intent(in) :: iy !< y index
    integer,      intent(in) :: iz !< z index
    
    integer :: normalized(1:3)

    !normalize to the min values
    normalized(1) = ix - pfft%local_i_start(1) + 1 
    normalized(2) = iy - pfft%local_i_start(2) + 1 
    normalized(3) = iz - pfft%local_i_start(3) + 1 
      
    index = (normalized(1)-1) + &                                      ! x component
         ((normalized(2)-1) * pfft%local_ni(1)) + &                    ! y component
         ((normalized(3)-1) * pfft%local_ni(1) * pfft%local_ni(2)) + 1 ! z component and the sum of 1
    
  end function cube_get_pfft_index
#endif

end module cube_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
