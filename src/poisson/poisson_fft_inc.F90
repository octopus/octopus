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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

  !-----------------------------------------------------------------

  subroutine X(poisson_fft_solve)(this, mesh, cube, pot, rho, mesh_cube_map, average_to_zero)
    type(poisson_fft_t),            intent(in)    :: this
    type(mesh_t),                   intent(in)    :: mesh
    type(cube_t),                   intent(in)    :: cube
    R_TYPE,                          intent(out)   :: pot(:)
    R_TYPE,                          intent(in)    :: rho(:)
    type(mesh_cube_parallel_map_t), intent(in)    :: mesh_cube_map
    logical,              optional, intent(in)    :: average_to_zero !< default is false

    logical :: average_to_zero_
    R_TYPE :: average
    type(cube_function_t) :: cf

    PUSH_SUB(X(poisson_fft_solve))

    average_to_zero_ = .false.
    if (present(average_to_zero)) average_to_zero_ = average_to_zero
    average = M_z0 !this avoids a non-initialized warning

#ifdef R_TCOMPLEX
    !If we perform complex ffts, the full cube need to be allocated
    ASSERT(cube%rs_n_global(1)==cube%fs_n_global(1))
#endif

    call cube_function_null(cf)
    call X(cube_function_alloc_RS)(cube, cf, in_device = (this%kernel /= POISSON_FFT_KERNEL_CORRECTED))

    ! put the density in the cube
    if (cube%parallel_in_domains) then
      call X(mesh_to_cube_parallel)(mesh, rho, cube, cf, mesh_cube_map)
    else
      if(mesh%parallel_in_domains) then
        call X(mesh_to_cube)(mesh, rho, cube, cf, local = .true.)
      else
        call X(mesh_to_cube)(mesh, rho, cube, cf)
      end if
    end if

    ! apply the Couloumb term in Fourier space
    call X(fourier_space_op_apply)(this%coulb, cube, cf)

    !now the cube has the potential
    if(average_to_zero_) average = X(cube_function_surface_average)(cube, cf)

    ! move the potential back to the mesh
    if (cube%parallel_in_domains) then
      call X(cube_to_mesh_parallel)(cube, cf, mesh, pot, mesh_cube_map)
    else
      if(mesh%parallel_in_domains) then
        call X(cube_to_mesh)(cube, cf, mesh, pot, local=.true.)
      else
        call X(cube_to_mesh)(cube, cf, mesh, pot)
      end if
    end if

    if(average_to_zero_) pot(1:mesh%np) = pot(1:mesh%np) - average

    call X(cube_function_free_RS)(cube, cf) ! memory is no longer needed

    POP_SUB(X(poisson_fft_solve))
  end subroutine X(poisson_fft_solve)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
