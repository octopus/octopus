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
!! $Id$

#include "global.h"

module cube_function_m
  use datasets_m
  use global_m
  use index_m
  use mesh_m
  use mesh_cube_map_m
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
  public ::                        &
    cube_function_t,               &
    cube_function_init,            &
    cube_function_init_from,       &
    cube_function_end,             &
    dcube_function_alloc_RS,       &
    zcube_function_alloc_RS,       &
    dcube_function_free_RS,        &
    zcube_function_free_RS,        &
    cube_function_surface_average, &
    cube_function_phase_factor,    &
#ifdef HAVE_PFFT
    dmesh_to_cube_parallel,        &
    zmesh_to_cube_parallel,        &
    dcube_to_mesh_parallel,        &
    zcube_to_mesh_parallel,        &
    cube_get_pfft_index,           &
#endif
    dmesh_to_cube,                 &
    zmesh_to_cube,                 &
    dcube_to_mesh,                 &
    zcube_to_mesh

  type cube_function_t
    integer :: n(1:3)               !< the linear dimensions of the cube
    integer :: nalloc(1:3)          !< local memory 
    integer :: nglobal(1:3)         !< global dimensions
    integer :: nprocs(1:3)          !< number of processors
    integer :: iprocs(1:3)          !< 
    integer :: offset(1:3)          !< 
    FLOAT, pointer :: dRS(:, :, :)  !< real-space grid
    CMPLX, pointer :: zRS(:, :, :)  !< real-space grid, complex numbers
    CMPLX, pointer :: FS(:, :, :)   !< Fourier-space grid
    integer :: fft_library          !< which FFT library has to be used. Options FFTW3=0, PFFT=1
    integer :: nx     ! = n(1)/2 + 1, first dimension of the FS array
    type(fft_t), pointer :: fft
#ifdef HAVE_PFFT
    type(pfft_t), pointer :: pfft
#endif

    type(mpi_grp_t) :: mpi_grp

  end type cube_function_t

  type(profile_t), save :: prof_m2c, prof_c2m
  
  integer, public, parameter :: &
       FFTW3_LIB = 0, &
       PFFT_LIB  = 1

contains

  ! ---------------------------------------------------------
  ! This function calculates the surface average of any function.
  ! \warning: Some more careful testing should be done on this.
  FLOAT function cube_function_surface_average(cf) result(x)
    type(cube_function_t), intent(in)       :: cf

    integer ix, iy, iz, npoints
    x = M_ZERO

    PUSH_SUB(cube_function_surface_average)

    do iy = 2, cf%n(2) - 1
      do iz = 2, cf%n(3) - 1
        x = x + (cf%dRS(1, iy, iz) + cf%dRS(cf%n(1), iy, iz))
      end do
    end do

    do ix = 2, cf%n(1) - 1
      do iz = 2, cf%n(3) - 1
        x = x + (cf%dRS(ix, 1, iz) + cf%dRS(ix, cf%n(2), iz))
      end do
    end do

    do ix = 2, cf%n(1) - 1
      do iy = 2, cf%n(2) - 1
        x = x + (cf%dRS(ix, iy, 1) + cf%dRS(ix, iy, cf%n(3)))
      end do
    end do

    do iz = 2, cf%n(3) - 1
      x = x + cf%dRS(1, 1, iz) + cf%dRS(cf%n(1), 1, iz) + &
        cf%dRS(1, cf%n(2), iz) + cf%dRS(cf%n(1), cf%n(2), 1)
    end do

    do iy = 2, cf%n(2) - 1
      x = x + cf%dRS(1, iy, 1) + cf%dRS(cf%n(1), iy, 1) + &
        cf%dRS(1, iy, cf%n(3)) + cf%dRS(cf%n(1), iy, cf%n(3))
    end do

    do ix = 2, cf%n(1) - 1
      x = x + cf%dRS(ix, 1, 1) + cf%dRS(ix, cf%n(2), 1) + &
        cf%dRS(ix, 1, cf%n(3)) + cf%dRS(ix, cf%n(2), cf%n(3))
    end do

    x = x + cf%dRS(1, 1, 1)             + cf%dRS(cf%n(1), 1, 1) + &
      cf%dRS(1, cf%n(2), 1)       + cf%dRS(cf%n(1), cf%n(2), 1) + &
      cf%dRS(1, 1, cf%n(3))       + cf%dRS(cf%n(1), 1, cf%n(3)) + &
      cf%dRS(1, cf%n(2), cf%n(3)) + cf%dRS(cf%n(1), cf%n(2), cf%n(3))

    npoints = 2*(cf%n(1)-2)**2 + 4*(cf%n(1)-2) + &
              2*(cf%n(2)-2)**2 + 4*(cf%n(2)-2) + &
              2*(cf%n(3)-2)**2 + 4*(cf%n(3)-2) + 8
    x = x/npoints

    POP_SUB(cube_function_surface_average)
  end function cube_function_surface_average

  ! ---------------------------------------------------------
  ! this routine computes
  ! cube_function_o = cf_o + exp(-k vec) cf_i
  subroutine cube_function_phase_factor(mesh, vec, cf_i, cf_o)
    type(mesh_t),       intent(in)    :: mesh
    FLOAT,              intent(in)    :: vec(3)
    type(cube_function_t), intent(in)    :: cf_i
    type(cube_function_t), intent(inout) :: cf_o

    CMPLX   :: k(3)
    integer :: n(3), ix, iy, iz, ixx, iyy, izz

    PUSH_SUB(cube_function_phase_factor)

    ASSERT(all(cf_i%n == cf_o%n))
    ASSERT(associated(cf_i%FS).and.associated(cf_o%FS))

    k = M_z0
    k(1:mesh%sb%dim) = M_zI * ((M_TWO*M_Pi)/(cf_i%n(1:mesh%sb%dim)*mesh%spacing(1:mesh%sb%dim)))

    n  = cf_i%n
    do iz = 1, n(3)
      izz = pad_feq(iz, n(3), .true.)
      do iy = 1, n(2)
        iyy = pad_feq(iy, n(2), .true.)
        do ix = 1, cf_i%nx
          ixx = pad_feq(ix, n(1), .true.)

          cf_o%FS(ix, iy, iz) = cf_o%FS(ix, iy, iz) + &
            exp( -(k(1)*vec(1)*ixx + k(2)*vec(2)*iyy + k(3)*vec(3)*izz) ) * cf_i%FS(ix, iy, iz)
        end do
      end do
    end do

    POP_SUB(cube_function_phase_factor)
  end subroutine cube_function_phase_factor
  
  ! ---------------------------------------------------------
  ! The following routines handle creation/destruction of the cube
  subroutine cube_function_init(cf, n)
    type(cube_function_t), intent(out) :: cf
    integer,               intent(in)  :: n(3)
    
    integer :: default_fft_library
    
    PUSH_SUB(cube_function_init)

    ASSERT(all(n(1:3) > 0))    

    nullify(cf%zRS)
    nullify(cf%dRS)
    nullify(cf%FS)
    cf%n = n
    cf%offset = 0
    
    nullify(cf%fft)

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
    
    default_fft_library = FFTW3_LIB

    call parse_integer(datasets_check('FFTLibrary'), default_fft_library, cf%fft_library)
#ifndef HAVE_PFFT
    if (cf%fft_library == PFFT_LIB) then
      write(message(1),'(a)')'You have selected the PFFT for FFT, but it is not compiled.'
      call messages_fatal(1)
    end if
#else
    nullify(cf%pfft)
#endif
    
    POP_SUB(cube_function_init) 
  end subroutine cube_function_init
  
  ! ---------------------------------------------------------

  subroutine cube_function_init_from(cf, cf_i)
    type(cube_function_t), intent(out) :: cf
    type(cube_function_t), intent(in)  :: cf_i

    PUSH_SUB(cube_function_init_from)
    ASSERT(all(cf_i%n>0))

    nullify(cf%dRS)
    nullify(cf%zRS)
    nullify(cf%FS)

    cf%n = cf_i%n

    if(associated(cf_i%fft)) then
      SAFE_ALLOCATE(cf%fft)
      call fft_copy(cf_i%fft, cf%fft)
      cf%nx = cf_i%nx
    else
      nullify(cf%fft)
    end if

    POP_SUB(cube_function_init_from)
  end subroutine cube_function_init_from

  ! ---------------------------------------------------------
  subroutine cube_function_end(cf)
    type(cube_function_t), intent(inout) :: cf
    
    PUSH_SUB(cube_function_end)
    
    SAFE_DEALLOCATE_P(cf%dRS)
    SAFE_DEALLOCATE_P(cf%zRS)
    SAFE_DEALLOCATE_P(cf%FS)
    
#ifdef HAVE_PFFT
    if(associated(cf%pfft)) then
      call pfft_end(cf%pfft)
      SAFE_DEALLOCATE_P(cf%pfft)
    end if
#else      

    if(associated(cf%fft)) then
      call fft_end(cf%fft)
      SAFE_DEALLOCATE_P(cf%fft)
    end if
#endif 
    
    POP_SUB(cube_function_end)
  end subroutine cube_function_end

#ifdef HAVE_PFFT
  !> returns the local index for the PFFT library using the global x, y and z
  integer function cube_get_pfft_index(pfft, ix, iy, iz) result(index)
    type(pfft_t), intent(in) :: pfft
    integer, intent(in) :: ix !< x index
    integer, intent(in) :: iy !< y index
    integer, intent(in) :: iz !< z index
    
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

#include "undef.F90"
#include "real.F90"
#include "cube_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "cube_function_inc.F90"

end module cube_function_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
