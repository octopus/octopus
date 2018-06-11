!! Copyright (C) 2018 N. Tancogne-Dejean
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

module kerker_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use fft_oct_m
  use fourier_space_oct_m
  use global_oct_m
  use mesh_oct_m
  use mesh_cube_parallel_map_oct_m
  use messages_oct_m
  use parser_oct_m
  use poisson_oct_m
  use poisson_fft_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                        &
    kerker_t,                      &
    kerker_init,                   &
    kerker_end,                    &
    kerker_filtering

 type kerker_t
   private
   FLOAT :: q0   ! Value of the Kerker filtering, in Bohr-1
   FLOAT :: q1   ! Metric parameter, in Bohr-1

   type(cube_t), pointer :: cube
   type(mesh_cube_parallel_map_t), pointer :: mesh_cube_map
   type(fourier_space_op_t) :: kerker_op  !< object for Fourier space operations
   type(fourier_space_op_t) :: metric_op
 end type kerker_t


contains

  ! ---------------------------------------------------------
  subroutine kerker_init(this, mesh, poisson)
    type(kerker_t),          intent(inout) :: this
    type(mesh_t),            intent(in)    :: mesh
    type(poisson_t), target, intent(in)    :: poisson

    PUSH_SUB(kerker_init)

    !%Variable KerkerParameter
    !%Type float
    !%Default 0.0
    !%Section SCF::Mixing
    !%Description
    !% (Experimental) Specifies the values of the Kerker preconditioner, which is designed to avoid
    !% charge sloshing, see for instance Kresse and Furthmueller, PRB 54 11169 (1996).
    !% Not used if set to zero.
    !%End
    call parse_variable('KerkerParameter', M_ZERO, this%q0, unit_one/units_inp%length)
    if(abs(this%q0) > CNST(1.0e-6)) call messages_experimental('KerkerParameter')

    !%Variable MetricParameter
    !%Type float
    !%Default 0.0
    !%Section SCF::Mixing
    !%Description
    !% (Experimental) Specifies the values of the metric preconditioner, which is designed to avoid
    !% charge sloshing, see Kresse and Furthmueller, PRB 54 11169 (1996).
    !% Not used if set to zero.
    !%End
    call parse_variable('MetricParameter', M_ZERO, this%q1, unit_one/units_inp%length)
    if(abs(this%q1) > CNST(1.0e-6)) call messages_experimental('MetricParameter')


    nullify(this%cube)
    nullify(this%mesh_cube_map)

    if(abs(this%q0) > CNST(1.0e-6) .or. (abs(this%q1) > CNST(1.0e-6) ) then

      ASSERT(associated(poisson%cube%fft))
      if(poisson%cube%fft%library /= FFTLIB_FFTW) then
        call messages_not_implemented('Kerker preconditioning without FFTW')
      end if

      ASSERT(poisson%fft_solver%kernel /= POISSON_FFT_KERNEL_CORRECTED)

      call build_kerker(poisson%cube)

      this%cube => poisson%cube
      this%mesh_cube_map => poisson%mesh_cube_map
    end if

    POP_SUB(kerker_init)

    contains

    subroutine build_kerker(cube)
      type(cube_t),  intent(in) :: cube
 
      integer :: ix, iy, iz, ixx(3), db(3)
      integer :: ix0; iy0, iz0
      FLOAT :: temp(3), modg2, mdog2_min
      FLOAT :: gg(3)
      FLOAT, allocatable :: kerker_FS(:,:,:), metric_FC(:,:,:)

      PUSH_SUB(kerker_init.build_kerker)

      db(1:3) = cube%rs_n_global(1:3)

      ! Kerker preconditioned in Fourier space
      SAFE_ALLOCATE(kerker_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
      kerker_FS = M_ZERO
      SAFE_ALLOCATE(metric_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
      metric_FS = M_ZERO

      temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))

      mdog2_min = M_HUGE

      do ix = 1, cube%fs_n_global(1)
        ixx(1) = pad_feq(ix, db(1), .true.)
        do iy = 1, cube%fs_n_global(2)
          ixx(2) = pad_feq(iy, db(2), .true.)
          do iz = 1, cube%fs_n_global(3)
            ixx(3) = pad_feq(iz, db(3), .true.)

            gg(1:3) =  ixx(1:3) * temp(1:3)
            gg(1:3) = matmul(mesh%sb%klattice_primitive(1:3,1:3),gg(1:3))
            modg2 = sum(gg(1:3)**2)

            kerker_FS(ix, iy, iz) = modg2/(modg2 + this%q0**2)

            metric_FS(ix, iy, iz) = (modg2 + this%q1**2)/modg2
            if(modg2 > M_EPSILON .and. modg2 < mdog2_min) mdog2_min = modg2
            if(modg2 <= M_EPSILON) then
              ix0 = ix; iy0 = iy; iz0 = iz
            end if
          end do
        end do
      end do

      metric_FS(ix0, iy0, iz0) = (modg2_min + this%q1**2)/modg2_min
      

      call dfourier_space_op_init(this%kerker_op, cube, kerker_FS)
      call dfourier_space_op_init(this%metric_op, cube, kerker_FS)

      SAFE_DEALLOCATE_A(kerker_FS)

      POP_SUB(kerker_init.build_kerker)
    end subroutine build_kerker
    
  end subroutine kerker_init

  ! ---------------------------------------------------------
  subroutine kerker_end(this)
    type(kerker_t),  intent(inout) :: this

    PUSH_SUB(kerker_end)

    if(associated(this%cube)) then
      nullify(this%cube)
      nullify(this%mesh_cube_map)
      call fourier_space_op_end(this%kerker_op)
    end if

    POP_SUB(kerker_end)
  end subroutine kerker_end

  ! --------------------------------------------------------

  subroutine kerker_filtering(this, mesh, rho, rho_out)
    type(kerker_t),intent(in) :: this
    type(mesh_t),  intent(in) :: mesh
    FLOAT,      intent(inout) :: rho(:)
    FLOAT,        intent(out) :: rho_out(:)

    type(cube_function_t) :: cf

    PUSH_SUB(kerker_filtering)

    if(abs(this%q0) <= CNST(1.0e-6)) then
      rho_out = rho
      POP_SUB(kerker_filtering)
      return
    end if

    call cube_function_null(cf)
    call dcube_function_alloc_RS(this%cube, cf, in_device = .true.)

    ! put the density in the cube
    if (this%cube%parallel_in_domains) then
      call dmesh_to_cube_parallel(mesh, rho, this%cube, cf, this%mesh_cube_map)
    else
      if(mesh%parallel_in_domains) then
        call dmesh_to_cube(mesh, rho, this%cube, cf, local = .true.)
      else
        call dmesh_to_cube(mesh, rho, this%cube, cf)
      end if
    end if

    ! apply the Kerker preconditioner in Fourier space
    call dfourier_space_op_apply(this%kerker_op, this%cube, cf)

    ! move the potential back to the mesh
    if (this%cube%parallel_in_domains) then
      call dcube_to_mesh_parallel(this%cube, cf, mesh, rho_out, this%mesh_cube_map)
    else
      if(mesh%parallel_in_domains) then
        call dcube_to_mesh(this%cube, cf, mesh, rho_out, local=.true.)
      else
        call dcube_to_mesh(this%cube, cf, mesh, rho_out)
      end if
    end if

    call dcube_function_free_RS(this%cube, cf)

    POP_SUB(kerker_filtering)

  end subroutine kerker_filtering

end module kerker_oct_m
  
