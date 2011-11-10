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
!! $Id$

#include "global.h"

module cube_function_m
  use cube_m
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
    cube_function_null,            &
    cube_function_end,             &
    dcube_function_alloc_RS,       &
    zcube_function_alloc_RS,       &
    dcube_function_free_RS,        &
    zcube_function_free_RS,        &
    cube_function_surface_average, &
    cube_function_surface_average_parallel, &
    cube_function_phase_factor,    &
    dmesh_to_cube_parallel,        &
    zmesh_to_cube_parallel,        &
    dcube_to_mesh_parallel,        &
    zcube_to_mesh_parallel,        &
    dmesh_to_cube,                 &
    zmesh_to_cube,                 &
    dcube_to_mesh,                 &
    zcube_to_mesh

  type cube_function_t
    FLOAT, pointer :: dRS(:, :, :)  !< real-space grid
    CMPLX, pointer :: zRS(:, :, :)  !< real-space grid, complex numbers
    CMPLX, pointer :: FS(:, :, :)   !< Fourier-space grid

    CMPLX, pointer :: pRS(:) !< real-space grid, complex number, in parallel
    CMPLX, pointer :: pFS(:) !< fourier-space grid, in parallel
  end type cube_function_t

  type(profile_t), save :: prof_m2c, prof_c2m
  
contains

  ! ---------------------------------------------------------
  ! This function calculates the surface average of any function.
  ! \warning: Some more careful testing should be done on this.
  FLOAT function cube_function_surface_average(cube, cf) result(x)
    type(cube_t),          intent(in) :: cube
    type(cube_function_t), intent(in) :: cf

    integer ix, iy, iz, npoints
    x = M_ZERO

    PUSH_SUB(cube_function_surface_average)

    do iy = 2, cube%n(2) - 1
      do iz = 2, cube%n(3) - 1
        x = x + (cf%dRS(1, iy, iz) + cf%dRS(cube%n(1), iy, iz))
      end do
    end do

    do ix = 2, cube%n(1) - 1
      do iz = 2, cube%n(3) - 1
        x = x + (cf%dRS(ix, 1, iz) + cf%dRS(ix, cube%n(2), iz))
      end do
    end do

    do ix = 2, cube%n(1) - 1
      do iy = 2, cube%n(2) - 1
        x = x + (cf%dRS(ix, iy, 1) + cf%dRS(ix, iy, cube%n(3)))
      end do
    end do

    do iz = 2, cube%n(3) - 1
      x = x + cf%dRS(1, 1, iz) + cf%dRS(cube%n(1), 1, iz) + &
        cf%dRS(1, cube%n(2), iz) + cf%dRS(cube%n(1), cube%n(2), 1)
    end do

    do iy = 2, cube%n(2) - 1
      x = x + cf%dRS(1, iy, 1) + cf%dRS(cube%n(1), iy, 1) + &
        cf%dRS(1, iy, cube%n(3)) + cf%dRS(cube%n(1), iy, cube%n(3))
    end do

    do ix = 2, cube%n(1) - 1
      x = x + cf%dRS(ix, 1, 1) + cf%dRS(ix, cube%n(2), 1) + &
        cf%dRS(ix, 1, cube%n(3)) + cf%dRS(ix, cube%n(2), cube%n(3))
    end do

    x = x + cf%dRS(1, 1, 1)           + cf%dRS(cube%n(1), 1, 1) + &
      cf%dRS(1, cube%n(2), 1)         + cf%dRS(cube%n(1), cube%n(2), 1) + &
      cf%dRS(1, 1, cube%n(3))         + cf%dRS(cube%n(1), 1, cube%n(3)) + &
      cf%dRS(1, cube%n(2), cube%n(3)) + cf%dRS(cube%n(1), cube%n(2), cube%n(3))

    npoints = 2*(cube%n(1)-2)**2 + 4*(cube%n(1)-2) + &
              2*(cube%n(2)-2)**2 + 4*(cube%n(2)-2) + &
              2*(cube%n(3)-2)**2 + 4*(cube%n(3)-2) + 8
    x = x/npoints

    POP_SUB(cube_function_surface_average)
  end function cube_function_surface_average

  FLOAT function cube_function_surface_average_parallel(cube, cf) result(x)
    type(cube_t),          intent(in) :: cube
    type(cube_function_t), intent(in) :: cf

    integer ii, jj, kk, ix, iy, iz, npoints
    FLOAT :: tmp_x

    PUSH_SUB(cube_function_surface_average_parallel)

    npoints = 0
    tmp_x = M_ZERO
    do ii = 1, cube%rs_n(1)
      do jj = 1, cube%rs_n(2)
        do kk = 1, cube%rs_n(3)
          ix = ii + cube%rs_istart(1) - 1
          iy = jj + cube%rs_istart(2) - 1
          iz = kk + cube%rs_istart(3) - 1
          if ( (ix == 1 .or. ix == cube%n(1)                                          ) .or. &
             ( (iy == 1 .or. iy == cube%n(2)) .and. (ix /= 1 .and. ix /= cube%n(1))   ) .or. &
             ( (iz == 1 .or. iz == cube%n(3)) .and. (ix /= 1 .and. ix /= cube%n(1) .and. iy /= 1 .and. iy /= cube%n(2))) ) then
            tmp_x = tmp_x + real(cf%pRS(cube_global2local(cube, ix, iy, iz)))
          end if
        end do
      end do
    end do

#ifdef HAVE_MPI
    call MPI_Allreduce(tmp_x, x, 1, MPI_FLOAT, MPI_SUM, cube%mpi_grp%comm, mpi_err)
#endif

    npoints = 2*(cube%n(1)-2)**2 + 4*(cube%n(1)-2) + &
              2*(cube%n(2)-2)**2 + 4*(cube%n(2)-2) + &
              2*(cube%n(3)-2)**2 + 4*(cube%n(3)-2) + 8
    x = x/npoints

    POP_SUB(cube_function_surface_average_parallel)
  end function cube_function_surface_average_parallel

  ! ---------------------------------------------------------
  ! this routine computes
  ! cube_function_o = cf_o + exp(-k vec) cf_i
  subroutine cube_function_phase_factor(mesh, vec, cube, cf_i, cf_o)
    type(mesh_t),          intent(in)    :: mesh
    FLOAT,                 intent(in)    :: vec(3)
    type(cube_t),          intent(in)    :: cube
    type(cube_function_t), intent(in)    :: cf_i
    type(cube_function_t), intent(inout) :: cf_o

    CMPLX   :: k(3)
    integer :: n(3), ix, iy, iz, ixx, iyy, izz

    PUSH_SUB(cube_function_phase_factor)

    ASSERT(associated(cf_i%FS).and.associated(cf_o%FS))

    k = M_z0
    k(1:mesh%sb%dim) = M_zI * ((M_TWO*M_Pi)/(cube%n(1:mesh%sb%dim)*mesh%spacing(1:mesh%sb%dim)))

    n  = cube%n
    do iz = 1, n(3)
      izz = pad_feq(iz, n(3), .true.)
      do iy = 1, n(2)
        iyy = pad_feq(iy, n(2), .true.)
        do ix = 1, cube%nx
          ixx = pad_feq(ix, n(1), .true.)

          cf_o%FS(ix, iy, iz) = cf_o%FS(ix, iy, iz) + &
            exp( -(k(1)*vec(1)*ixx + k(2)*vec(2)*iyy + k(3)*vec(3)*izz) ) * cf_i%FS(ix, iy, iz)
        end do
      end do
    end do

    POP_SUB(cube_function_phase_factor)
  end subroutine cube_function_phase_factor
  
  ! ---------------------------------------------------------
  subroutine cube_function_null(cf)
    type(cube_function_t), intent(out) :: cf
    
    PUSH_SUB(cube_function_null)

    nullify(cf%zRS)
    nullify(cf%dRS)
    nullify(cf%FS)

    nullify(cf%pRS)
    nullify(cf%pFS)

    POP_SUB(cube_function_null) 
  end subroutine cube_function_null

  ! ---------------------------------------------------------
  subroutine cube_function_end(cf)
    type(cube_function_t), intent(inout) :: cf
    
    PUSH_SUB(cube_function_end)
    
    SAFE_DEALLOCATE_P(cf%dRS)
    SAFE_DEALLOCATE_P(cf%zRS)
    SAFE_DEALLOCATE_P(cf%FS)

    nullify(cf%pRS)
    nullify(cf%pFS)

    POP_SUB(cube_function_end)
  end subroutine cube_function_end


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
