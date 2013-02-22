!! Copyright (C) 2013 J. Alberdi-Rodriguez
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

module poisson_libisf_m
  use cube_function_m
  use cube_m
  use datasets_m
  use fft_m
  use fourier_space_m
  use geometry_m
  use global_m
  use lalg_basic_m
  use loct_m
  use loct_math_m
  use math_m
  use mesh_cube_parallel_map_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use par_vec_m
  use parser_m
  use profiling_m
  use simul_box_m
  use splines_m
  use unit_m
  use unit_system_m

#ifdef HAVE_ISF
  !! From BigDFT
  use poisson_solver
#endif

  implicit none

  private
  public ::                        &
    poisson_libisf_t,              &
    poisson_libisf_init,           &
    poisson_libisf_end,            &
    poisson_libisf_global_solve,   &
    poisson_libisf_parallel_solve, &
    poisson_libisf_get_dims

  type poisson_libisf_t
    type(fourier_space_op_t) :: coulb  !< object for Fourier space operations
    FLOAT, pointer           :: kernel(:) !< choice of kernel, one of options above
    !> Indicates the boundary conditions (BC) of the problem:
    !!            'F' free BC, isolated systems.
    !!                The program calculates the solution as if the given density is
    !!                "alone" in R^3 space.
    !!            'S' surface BC, isolated in y direction, periodic in xz plane                
    !!                The given density is supposed to be periodic in the xz plane,
    !!                so the dimensions in these direction mus be compatible with the FFT
    !!                Beware of the fact that the isolated direction is y!
    !!            'P' periodic BC.
    !!                The density is supposed to be periodic in all the three directions,
    !!                then all the dimensions must be compatible with the FFT.
    !!                No need for setting up the kernel.
    character(len = 1) :: geocode  = "F" !< 'F' free boundary contition
    !> Indicates the distribution of the data of the input/output array:
    !!            'G' global data. Each process has the whole array of the density 
    !!                which will be overwritten with the whole array of the potential
    !!            'D' distributed data. Each process has only the needed part of the density
    !!                and of the potential. The data distribution is such that each processor
    !!                has the xy planes needed for the calculation AND for the evaluation of the 
    !!                gradient, needed for XC part, and for the White-Bird correction, which
    !!                may lead up to 8 planes more on each side. Due to this fact, the information
    !!                between the processors may overlap.
    character(len = 1) :: datacode = "G" 

    integer :: isf_order !< order of the interpolating scaling functions used in the decomposition 
    integer :: localnscatterarr(5)
    
    integer     :: rs_n_global(3) !< total size of the fft in each direction in real space
    integer     :: rs_istart(1:3) !< where does the local portion of the function start in real space

  end type poisson_libisf_t
contains

  subroutine poisson_libisf_init(this, mesh, cube, soft_coulb_param)
    type(poisson_libisf_t), intent(out)   :: this
    type(mesh_t),        intent(inout) :: mesh
    type(cube_t),        intent(inout) :: cube
    FLOAT, optional,     intent(in)    :: soft_coulb_param

#ifdef HAVE_LIBISF
    logical data_is_parallel
 
    PUSH_SUB(poisson_libisf_init)
    
    ! Free BC
    this%geocode = "F" 
    this%isf_order = 16

    !%Variable PoissonSolverISFParallelData
    !%Type float
    !%Section Hamiltonian::Poisson
    !%Description
    !% Indicates the distribution of the data of the input/output array
    !%
    !% Default
    !% true
    !%Option true
    !% Z axis of the rho vector is split among the MPI processes
    !%Option false
    !% Entire vector is saved in all the MPI processes
    !%End
    call parse_logical('PoissonSolverISFParallelData', .false., data_is_parallel)
    if (data_is_parallel) then
      this%datacode = "D"
    else 
      this%datacode = "G"
    end if

    call createKernel(mpi_world%rank, mpi_world%size, this%geocode, &
         cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), &
         mesh%spacing(1), mesh%spacing(2), mesh%spacing(3), &
         this%isf_order, this%kernel, .true.)

    POP_SUB(poisson_libisf_init)
#endif
  end subroutine poisson_libisf_init

  
  !-----------------------------------------------------------------
  subroutine poisson_libisf_end(this)
    type(poisson_libisf_t), intent(inout) :: this

    PUSH_SUB(poisson_libisf_end)

    write(message(1),*) "memory has to be deallocated in LIBISF"
    call messages_warning(1)

    POP_SUB(poisson_libisf_end)
  end subroutine poisson_libisf_end

  subroutine poisson_libisf_parallel_solve(this, mesh, cube, pot, rho,  mesh_cube_map)
    type(poisson_libisf_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(in)    :: cube
    FLOAT,               intent(out)   :: pot(:)
    FLOAT,               intent(in)    :: rho(:)
    type(mesh_cube_parallel_map_t), intent(in)    :: mesh_cube_map

#ifdef HAVE_LIBISF
    type(profile_t), save :: prof
    type(cube_function_t) :: cf   
    double precision :: hartree_energy !<  Hartree energy
    double precision :: offset !< offset  Total integral on the supercell of the final potential on output
    !!                To be used only in the periodic case, ignored for other boundary conditions.

     
    double precision, allocatable :: pot_ion(:,:,:) !< pot_ion additional  external potential
    !that is added to the output
    !!                when the XC parameter ixc/=0 and sumpion=.true.
    !!                When sumpion=.true., it is always provided in the distributed form,
    !!                clearly without the overlapping terms which are needed only for the XC part

    double precision :: strten(6)
    PUSH_SUB(poisson_libisf_parallel_solve)

    call cube_function_null(cf)
    call dcube_function_alloc_RS(cube, cf)

    call dmesh_to_cube_parallel(mesh, rho, cube, cf, mesh_cube_map)
  
    call profiling_in(prof,"ISF_LIBRARY")
    call H_potential(this%geocode, this%datacode, mpi_world%rank, mpi_world%size, &
         cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), &
         mesh%spacing(1), mesh%spacing(2), mesh%spacing(3), &
         cf%dRS, this%kernel, pot_ion, hartree_energy, offset, .false., &
         quiet = "NO ", stress_tensor = strten) !optional argument
    call profiling_out(prof)

    call dcube_to_mesh_parallel(cube, cf, mesh, pot, mesh_cube_map)

    call dcube_function_free_RS(cube, cf)

    !! @TODO: All the other deallocations
    write(message(1),*) "memory has to be deallocated in LIBISF"
    call messages_warning(1)
    
    POP_SUB(poisson_libisf_parallel_solve)
#endif
  end subroutine poisson_libisf_parallel_solve

  !-----------------------------------------------------------------
  subroutine poisson_libisf_global_solve(this, mesh, cube, pot, rho, all_nodes)
    type(poisson_libisf_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(in)    :: cube
    FLOAT,               intent(out)   :: pot(:)
    FLOAT,               intent(in)    :: rho(:)
    logical,             intent(in)    :: all_nodes

#ifdef HAVE_LIBISF
    type(cube_function_t) :: cf
    double precision :: hartree_energy !<  Hartree energy
    double precision :: offset !< offset  Total integral on the supercell of the final potential on output
    !!                To be used only in the periodic case, ignored for other boundary conditions.

     
    double precision, allocatable ::  pot_ion(:,:,:) !< pot_ion additional  external potential
    !!                that is added to the output
    !!                when the XC parameter ixc/=0 and sumpion=.true.
    !!                When sumpion=.true., it is always provided in the distributed form,
    !!                clearly without the overlapping terms which are needed only for the XC part

    double precision :: strten(6)

    PUSH_SUB(poisson_libisf_global_solve)
    
    call cube_function_null(cf)
    call dcube_function_alloc_RS(cube, cf)

    if(mesh%parallel_in_domains) then
      call dmesh_to_cube(mesh, rho, cube, cf, local=.true.)
    else
      call dmesh_to_cube(mesh, rho, cube, cf)
    end if

    call H_potential(this%geocode, this%datacode, mpi_world%rank, mpi_world%size, &
         cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), &
         mesh%spacing(1), mesh%spacing(2), mesh%spacing(3), &
         cf%dRS, this%kernel, pot_ion, hartree_energy, offset, .false., &
         quiet = "NO ", stress_tensor = strten) !optional argument

    if(mesh%parallel_in_domains) then
      call dcube_to_mesh(cube, cf, mesh, pot, local=.true.)
    else
      call dcube_to_mesh(cube, cf, mesh, pot)
    end if

    call dcube_function_free_RS(cube, cf)

    !! @TODO: All the other deallocations
    write(message(1),*) "memory has to be deallocated in LIBISF"
    call messages_warning(1)

    POP_SUB(poisson_libisf_global_solve)
#endif
  end subroutine poisson_libisf_global_solve

  subroutine poisson_libisf_get_dims(this, cube) 
    type(poisson_libisf_t), intent(inout) :: this
    type(cube_t), intent(inout) :: cube

#ifdef HAVE_LIBISF
    !>    ixc         eXchange-Correlation code. Indicates the XC functional to be used 
    !!                for calculating XC energies and potential. 
    !!                ixc=0 indicates that no XC terms are computed. The XC functional codes follow
    !!                the ABINIT convention.   
    integer :: ixc = 1
    !>    n3d         third dimension of the density. For distributed data, it takes into account 
    !!                the enlarging needed for calculating the XC functionals.
    !!                For global data it is simply equal to n03. 
    !!                When there are too many processes and there is no room for the density n3d=0
    integer :: n3d 
    !>    n3p         third dimension for the potential. The same as n3d, but without 
    !!                taking into account the enlargment for the XC part. For non-GGA XC, n3p=n3d.
    integer :: n3p
    !>    n3pi        Dimension of the pot_ion array, always with distributed data. 
    !!                For distributed data n3pi=n3p
    integer :: n3pi
    !>!    i3xcsh     Shift of the density that must be performed to enter in the 
    !!                non-overlapping region. Useful for recovering the values of the potential
    !!                when using GGA XC functionals. If the density starts from rhopot(1,1,1),
    !!                the potential starts from rhopot(1,1,i3xcsh+1). 
    !!                For non-GGA XCs and for global distribution data i3xcsh=0
    integer :: i3xcsh
    !>    i3s         Starting point of the density effectively treated by each processor 
    !!                in the third direction.
    !!                It takes into account also the XC enlarging. The array rhopot will correspond
    !!                To the planes of third coordinate from i3s to i3s+n3d-1. 
    !!                The potential to the planes from i3s+i3xcsh to i3s+i3xcsh+n3p-1
    !!                The array pot_ion to the planes from i3s+i3xcsh to i3s+i3xcsh+n3pi-1
    !!                For global disposition i3s is equal to distributed case with i3xcsh=0.
    integer :: i3s

    PUSH_SUB(libisf_get_dims)

    !! Get the dimensions of the cube
    call PS_dim4allocation(this%geocode, this%datacode, mpi_world%rank, mpi_world%size, &
         cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), &
         ixc, n3d, n3p, n3pi, i3xcsh, i3s)
    this%localnscatterarr(:) = (/ n3d, n3p, n3pi, i3xcsh, i3s /)

    cube%rs_n(1:2)      = cube%rs_n_global(1:2)
    cube%rs_n(3)        = this%localnscatterarr(1)
    cube%rs_istart(1:2) = 1
    cube%rs_istart(3)   = this%localnscatterarr(5)
    
    !! With ISF we don`t care about the Fourier space and its dimensions
    !! We`ll put as in RS
    cube%fs_n_global(1) = cube%rs_n_global(1)
    cube%fs_n_global(2) = cube%rs_n_global(2)
    cube%fs_n_global(3) = cube%rs_n_global(3)
    cube%fs_n(1:2)      = cube%rs_n_global(1:2)
    cube%fs_n(3)        = this%localnscatterarr(1)
    cube%fs_istart(1:2) = 1
    cube%fs_istart(3)   = this%localnscatterarr(5)

    POP_SUB(libisf_get_dims)
#endif
  end subroutine poisson_libisf_get_dims
  
end module poisson_libisf_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
