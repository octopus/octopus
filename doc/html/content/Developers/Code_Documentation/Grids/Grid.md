---
Title: Grids
Weight: 1
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


The Grid in Octopus
====================

Introduction
------------


The fundamental element of Octopus is the real space grid, on which
all quantities, such as wafe functions, densities and potentials are
defined.

In general, a grid consists of

* the mesh
* the simulation box

The mesh is based on an regular mesh in _D_ dimensions, which can be distorted 
using so-called curvi-linear coordinates.

The underlying internal structure is a mesh of integer coordinates and the important
mapping arrays which 
* map the index to the _D_ integer coordinates
* map the integer coordinates to the index.

The grid is constructed in several steps:

1) Generation of a regular mesh filling a rectangular volume enclosing the simulation box.
2) Selecting inner points, enlargment points and ghost points. 
3) re-ordering the points and regeneration of the mapping arrays.




The data structures:
--------------------


The top level data structure, describing a grid is defined in
`grid/grid.F90`:

```Fortran
#include_type_def grid_t
```


The mesh descriptor is defined in 
`grid/mesh.F90`:

```Fortran
#include_type_def mesh_t
```
which uses the structure `index_t`, containing the range and the mapping arrays:
```Fortran
#include_type_def index_t
```

About lxyz and lxyz_inv maps:
-----------------------------

The direct map:

> lxyz(_d_, _i_) = R<sup>_(d)_</sup><sub>_i_</sub>  

where _d_ denotes the dimension (i.e. x, y or z) and _i_ is the index of the point.

The inverse map:

> lxyz_inv( R<sup>x</sup><sub>_i_</sub> , R<sup>y</sup><sub>_i_</sub> , R<sup>z</sup><sub>_i_</sub> ) 
= _i_

The points defined by lxyz define a rectangular box of _d_ dimensions.  
The real mesh vectors are related to R<sub>_i_</sub> by multiplication with spacing
and possibly distortion for curvilinear coordinates.

Note that this index array is the same in all parallel domains, and relates the integer coordinates to the global index of a given point.




```Fortran
#include_type_def mesh_cube_map_t
```


QUESTIONS:
----------


- What exactly is mesh_t::resolution(:,:,:)?

- How is the surface element(d) defined?

- What does the cube_map do?




Accessing functions
===================

Function, e.g. the density are given by an array rho(_i_,_spin_), where _i_ ranges over the mesh points associated with the given domain (see domain decomposition), representing

rho( **r**<sub>_i_</sub>,  _spin_ )


What is `mesh_x_global(mesh, ip)` ?

`mesh_x_global(mesh, ip)` returns a `FLOAT` vector, corresponding to the coordinates of point `ip` of the global mesh.


Integration
-----------

Integration is implemented as a straightforward summation:

```Fortran
 do ip = 1, mesh%np
      dd = dd + ff(ip)
 end do 
```

Differentiation
---------------

Taking derivatives is done by finite differences. The points involved
in a derivative are defined by the stencil (see below).

Derivatives are discussed in a separate document [Derivatives.md](Derivatives.md).


Note on packed states:
----------------------



Note on curvilinear meshes:
---------------------------

The `mesh::x(:,:)` array always contains a regular mesh, which gets 'distorded' to a curvilinear mesh by additional function calls.


Domain decomposition
--------------------

See [Domain_Decomposition.md](Domain_Decomposition.md).



--------------------


Setting up a grid
==================



The setup procedure is called in system/system.F90:system_init() and
consists of 3 steps which depend on other setup procedures taking place
in between those steps. The sequence is:

```Fortran

  SAFE_ALLOCATE(sys%gr)
  SAFE_ALLOCATE(sys%st)

  space_init(sys%space)                             ! define spatial dimensions
    
  geometry_init(sys%geo, sys%space)                 ! set up atomic positions

  grid_init_stage_0(sys%gr, sys%geo, sys%space)     ! set up simulation box from the geometry
  grid_init_stage_1(sys%gr, sys%geo)                ! 

  parallel_init()

  geometry_partition(sys%geo, sys%mc)
  kpoints_distribute(sys%st%d, sys%mc)
  states_distribute_nodes(sys%st, sys%mc)

  grid_init_stage_2(sys%gr, sys%mc, sys%geo)
```


In particular, we have:

* the "zero-th" stage of grid initialization. It initializes the simulation box.

```Fortran
  subroutine grid_init_stage_0(gr, geo, space)
    type(grid_t),          intent(inout) :: gr
    type(geometry_t),      intent(inout) :: geo
    type(space_t),         intent(in)    :: space

    call simul_box_init(gr%sb, geo, space)
      
  end subroutine grid_init_stage_0
```


* the "first stage":

  This sets up derivatives, double-grid and calls stages 1 and 2 of mesh_init.
At this stage we are still unaware or parallelism and domain decomposition.

```Fortran
  subroutine grid_init_stage_1(gr, geo)
    type(grid_t),      intent(inout) :: gr
    type(geometry_t),  intent(in)    :: geo

    ... [parse grid related input variables]

    ! initialize curvilinear coordinates
    call curvilinear_init(gr%cv, gr%sb, geo, grid_spacing)

    !> initialize derivatives (see Octopus_Derivatives.txt)
    !! this also defines the stencil and sets up the Laplacian and gradient
    
    call derivatives_init(gr%der, gr%sb, gr%cv%method /= CURV_METHOD_UNIFORM)


    call double_grid_init(gr%dgrid, gr%sb)

    enlarge = 0
    enlarge(1:gr%sb%dim) = 2
    enlarge = max(enlarge, double_grid_enlarge(gr%dgrid))
    enlarge = max(enlarge, gr%der%n_ghost)

    ! now we generate the mesh and the derivatives

    call mesh_init_stage_1(gr%mesh, gr%sb, gr%cv, grid_spacing, enlarge)

    ! the stencil used to generate the grid is a union of a cube (for
    ! multigrid) and the Laplacian.
    
    call stencil_cube_get_lapl(cube, gr%sb%dim, order = 2)
    call stencil_union(gr%sb%dim, cube, gr%der%lapl%stencil, gr%stencil)
    call stencil_end(cube)

    call mesh_init_stage_2(gr%mesh, gr%sb, geo, gr%cv, gr%stencil)

  end subroutine grid_init_stage_1
```

* the second stage:

```Fortran
  subroutine grid_init_stage_2(gr, mc, geo)
    type(grid_t), target, intent(inout) :: gr
    type(multicomm_t),    intent(in)    :: mc
    type(geometry_t),     intent(in)    :: geo

    PUSH_SUB(grid_init_stage_2)

    call mesh_init_stage_3(gr%mesh, gr%stencil, mc)

    call nl_operator_global_init()
    if(gr%have_fine_mesh) then
      message(1) = "Info: coarse mesh"
      call messages_info(1)
    end if
    call derivatives_build(gr%der, gr%mesh)

    ! initialize a finer mesh to hold the density, for this we use the
    ! multigrid routines
    
    if(gr%have_fine_mesh) then

      if(gr%mesh%parallel_in_domains) then
        message(1) = 'UseFineMesh does not work with domain parallelization.'
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(gr%fine%mesh)
      SAFE_ALLOCATE(gr%fine%der)
      
      call multigrid_mesh_double(geo, gr%cv, gr%mesh, gr%fine%mesh, gr%stencil)
      
      call derivatives_init(gr%fine%der, gr%mesh%sb, gr%cv%method /= CURV_METHOD_UNIFORM)
      
      call mesh_init_stage_3(gr%fine%mesh, gr%stencil, mc)
      
      call multigrid_get_transfer_tables(gr%fine%tt, gr%fine%mesh, gr%mesh)
      
      message(1) = "Info: fine mesh"
      call messages_info(1)
      call derivatives_build(gr%fine%der, gr%fine%mesh)

      gr%fine%der%coarser => gr%der
      gr%der%finer =>  gr%fine%der
      gr%fine%der%to_coarser => gr%fine%tt
      gr%der%to_finer => gr%fine%tt

    else
      gr%fine%mesh => gr%mesh
      gr%fine%der => gr%der
    end if

    call mesh_check_symmetries(gr%mesh, gr%mesh%sb)

    ! multigrids are not initialized by default
    nullify(gr%mgrid)

    ! print info concerning the grid
    call grid_write_info(gr, geo, stdout)

    POP_SUB(grid_init_stage_2)
  end subroutine grid_init_stage_2
```


Mesh initialization:
--------------------

* First stage: set up the mesh_t::idx index structure, defining the lower and upper bounds.

```Fortran
  subroutine mesh_init_stage_1(mesh, sb, cv, spacing, enlarge)
    type(mesh_t),                intent(inout) :: mesh
    type(simul_box_t),   target, intent(in)    :: sb
    type(curvilinear_t), target, intent(in)    :: cv
    FLOAT,                       intent(in)    :: spacing(1:MAX_DIM)
    integer,                     intent(in)    :: enlarge(MAX_DIM)

    integer :: idir, jj, delta
    FLOAT   :: x(MAX_DIM), chi(MAX_DIM), spacing_new(-1:1)
    logical :: out

    PUSH_SUB(mesh_init_stage_1)
    call profiling_in(mesh_init_prof, "MESH_INIT")

    mesh%sb => sb          ! keep an internal pointer
    mesh%spacing = spacing ! this number can change in the following
    mesh%use_curvilinear = cv%method /= CURV_METHOD_UNIFORM
    mesh%cv => cv

    ! multiresolution requires the curvilinear coordinates machinery
    mesh%use_curvilinear = mesh%use_curvilinear .or. sb%mr_flag

    mesh%idx%dim = sb%dim
    mesh%idx%is_hypercube = sb%box_shape == HYPERCUBE
    mesh%idx%enlarge = enlarge

    if(sb%mr_flag) mesh%idx%enlarge = mesh%idx%enlarge*(2**sb%hr_area%num_radii)

    ! adjust nr: determine the number of points inside a parallelepiped, containing the simulation box.

    mesh%idx%nr = 0

    do idir = 1, sb%dim            ! loop over dimensions
      chi = M_ZERO
      ! the upper border
      jj = 0
      out = .false.
      do while(.not.out)
        jj = jj + 1
	chi(idir) = real(jj, REAL_PRECISION)*mesh%spacing(idir)
	if ( mesh%use_curvilinear ) then
	  call curvilinear_chi2x(sb, cv, chi(1:sb%dim), x(1:sb%dim))
	  out = (x(idir) > nearest(sb%lsize(idir), M_ONE))
	else
	  out = (chi(idir) > nearest(sb%lsize(idir), M_ONE))
	end if
      end do
      mesh%idx%nr(2, idir) = jj - 1
    end do

    ! define lower bounds by mirroring:

    ! we have a symmetric mesh (for now)
    mesh%idx%nr(1,:) = -mesh%idx%nr(2,:)

    ! we have to adjust a couple of things for the periodic directions
    do idir = 1, sb%periodic_dim

      if(mesh%idx%nr(2, idir) == 0) then
        ! this happens if Spacing > box size
	mesh%idx%nr(2, idir) =  1
	mesh%idx%nr(1, idir) = -1
      end if

      ! We have to adjust the spacing to be commensurate with the box,
      ! for this we scan the possible values of the grid size around the
      ! one we selected. We choose the size that has the spacing closest
      ! to the requested one.
      do delta = -1, 1
        spacing_new(delta) = CNST(2.0)*sb%lsize(idir)/real(2*mesh%idx%nr(2, idir) + 1 - delta) 
	spacing_new(delta) = abs(spacing_new(delta) - spacing(idir))
      end do

      delta = minloc(spacing_new, dim = 1) - 2

      ASSERT(delta >= -1) 
      ASSERT(delta <=  1) 

      mesh%spacing(idir) = CNST(2.0)*sb%lsize(idir)/real(2*mesh%idx%nr(2, idir) + 1 - delta)

      ! we need to adjust the grid by adding or removing one point
      if(delta == -1) then
        mesh%idx%nr(1, idir) = mesh%idx%nr(1, idir) - 1
      else if (delta == 1) then
        mesh%idx%nr(2, idir) = mesh%idx%nr(2, idir) - 1
      end if

    end do ! idir

    if( any(abs(mesh%spacing(1:sb%periodic_dim) - spacing(1:sb%periodic_dim)) > CNST(1e-6)) ) then
      call messages_write('The spacing has been modified to make it commensurate with the periodicity of the system.')
      call messages_warning()
    end if

    do idir = sb%periodic_dim + 1, sb%dim
      if(mesh%idx%nr(2, idir) == 0) then
        write(message(1),'(a,i2)') 'Spacing > box size in direction ', idir
        call messages_fatal(1)
      end if
    end do

    mesh%idx%ll(1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) - mesh%idx%nr(1, 1:MAX_DIM) + 1

    call profiling_out(mesh_init_prof)
    POP_SUB(mesh_init_stage_1)

  end subroutine mesh_init_stage_1
```


* Second stage: 

```Fortran
! ---------------------------------------------------------
!> This subroutine checks if every grid point belongs to the internal
!! mesh, based on the global lxyz_inv matrix. Afterwards, it counts
!! how many points has the mesh and the enlargement.

subroutine mesh_init_stage_2(mesh, sb, geo, cv, stencil)
  type(mesh_t),        intent(inout) :: mesh
  type(simul_box_t),   intent(in)    :: sb
  type(geometry_t),    intent(in)    :: geo
  type(curvilinear_t), intent(in)    :: cv
  type(stencil_t),     intent(in)    :: stencil

  ! enlarge mesh for boundary points
  mesh%idx%nr(1, 1:MAX_DIM) = mesh%idx%nr(1, 1:MAX_DIM) - mesh%idx%enlarge(1:MAX_DIM)
  mesh%idx%nr(2, 1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) + mesh%idx%enlarge(1:MAX_DIM)
  
  if(mesh%idx%is_hypercube) then
    call hypercube_init(mesh%idx%hypercube, sb%dim, mesh%idx%nr, mesh%idx%enlarge(1))
    mesh%np_part_global = hypercube_number_total_points(mesh%idx%hypercube)
    mesh%np_global      = hypercube_number_inner_points(mesh%idx%hypercube)
    call profiling_out(mesh_init_prof)
    POP_SUB(mesh_init_stage_2)
    return
  end if

  nr = mesh%idx%nr

  ! allocate the xyz arrays
  SAFE_ALLOCATE(mesh%idx%lxyz_inv(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)))

  mesh%idx%lxyz_inv(:,:,:) = 0

  res = 1
  SAFE_ALLOCATE(xx(1:MAX_DIM, mesh%idx%nr(1,1):mesh%idx%nr(2,1)))
  SAFE_ALLOCATE(in_box(mesh%idx%nr(1,1):mesh%idx%nr(2,1)))
  chi = M_ZERO

  start_z = mesh%idx%nr(1, 3)
  end_z = mesh%idx%nr(2, 3)

  ! We label the points inside the mesh

  do iz = start_z, end_z
    chi(3) = real(iz, REAL_PRECISION) * mesh%spacing(3)
    do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      chi(2) = real(iy, REAL_PRECISION) * mesh%spacing(2)
      do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
        chi(1) = real(ix, REAL_PRECISION) * mesh%spacing(1)
        call curvilinear_chi2x(sb, cv, chi(:), xx(:, ix))   ! apply curvilinear distortion
      end do

      ! current mesh points in real cartesian coordinates are xx(1:dim, nr(1,1):nr(2,1))
      ! check whether xx is inside the simulation box (processing by lines for efficiency)

      call simul_box_in_box_vec(sb, geo, mesh%idx%nr(2,1) - mesh%idx%nr(1,1) + 1, xx, in_box)

      do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
        ! With multiresolution, only inner (not enlargement) points are marked now
        if(sb%mr_flag) then
	  ...
	else ! the usual way: mark both inner and enlargement points
	
          if (in_box(ix)) then

	    ! the point xx(ix) itself is inside the box: mark as INNER_POINT

            mesh%idx%lxyz_inv(ix, iy, iz) = ibset(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)

	    ! check whether the 'stencil surrounding' of xx(ix) is inside enlarged the 'lxyz box'
	    ! if so, mark that those surrounding points as ENLARGEMENT_POINT.

            do is = 1, stencil%size
              if(stencil%center == is) cycle
              ii = ix + stencil%points(1, is)
              jj = iy + stencil%points(2, is)
              kk = iz + stencil%points(3, is)
              if(any((/ii, jj, kk/) < mesh%idx%nr(1, 1:3)) .or. any((/ii, jj, kk/) >  mesh%idx%nr(2, 1:3))) cycle
              mesh%idx%lxyz_inv(ii, jj, kk) = ibset(mesh%idx%lxyz_inv(ii, jj, kk), ENLARGEMENT_POINT)
            end do
	    
          end if ! in_box(ix)
	  
        end if ! multi_resolution
	
      end do ! ix
    end do ! iy
  end do ! iz
  
  call profiling_out(prof)

  SAFE_DEALLOCATE_A(xx)
  SAFE_DEALLOCATE_A(in_box)


  ! count the points
  il = 0
  ik = 0
  
  do iz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
    do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
        if(btest(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)) ik = ik + 1
        if(mesh%idx%lxyz_inv(ix, iy, iz) /= 0) il = il + 1
      end do
    end do
  end do
  
  mesh%np_part_global = il ! inner and enlargement points
  mesh%np_global      = ik ! inner points only

  ASSERT(mesh%np_global > 0)
  ASSERT(mesh%np_part_global > 0)

  ... [multiresolution code removed for simplicity]

end subroutine mesh_init_stage_2
```



* Third stage:

  We can finally generate the list of mesh points mesh_t::x.
At this stage, domain decomposition will be considered and x(:,:) contains
only the domain local points.

```Fortran
! ---------------------------------------------------------
!> When running parallel in domains, stencil and np_stencil
!! are needed to compute the ghost points.
!! mpi_grp is the communicator group that will be used for
!! this mesh.
! ---------------------------------------------------------

subroutine mesh_init_stage_3(mesh, stencil, mc, parent)
  type(mesh_t),              intent(inout) :: mesh      !< current mesh to be completed
  type(stencil_t),           intent(in)    :: stencil   !< combined stencil for mesh generation
  type(multicomm_t),         intent(in)    :: mc        !< information about parallelism
  type(mesh_t),    optional, intent(in)    :: parent    !< ???


  call mpi_grp_init(mpi_grp, mc%group_comm(P_STRATEGY_DOMAINS))
  
  ! check if we are running in parallel in domains
  mesh%parallel_in_domains = (mpi_grp%size > 1)

  if(.not. mesh%parallel_in_domains) then
    ! When running parallel, x is computed later.
    SAFE_ALLOCATE(mesh%x(1:mesh%np_part_global, 1:MAX_DIM))  ! note, this only contains inner and enlargement points.
  end if
  
  if(.not. mesh%idx%is_hypercube) then
    call create_x_lxyz()                                             ! general case, e.g. for minimal simulation box.
    	 							     ! NOTE: This re-defines mesh_t:idx !!!
  else if(.not. mesh%parallel_in_domains) then
    do ip = 1, mesh%np_part_global
      mesh%x(ip, 1:MAX_DIM) = mesh_x_global(mesh, ip, force=.true.)
    end do
  end if

  mesh%mpi_grp = mpi_grp 
  
  if(mesh%parallel_in_domains) then

    call do_partition() ! this function is defined inside mesh_init_stage_3 and can access all variables.

  else

    ! When running serially those two are the same.
    mesh%np      = mesh%np_global
    mesh%np_part = mesh%np_part_global

    ! These must be initialized for vec_gather, vec_scatter to work
    ! as copy operations when running without domain parallelization.
    mesh%vp%np_global = mesh%np_global
    mesh%vp%npart = 1
    mesh%vp%xlocal = 1
  end if

  call mesh_cube_map_init(mesh%cube_map, mesh%idx, mesh%np_global)

  call mesh_get_vol_pp(mesh%sb)

  call profiling_out(mesh_init_prof)
  POP_SUB(mesh_init_stage_3)

end subroutine mesh_init_stage_3
```


This uses (for non-hypercubic meshes):

```Fortran
  subroutine create_x_lxyz()
    integer :: il, iin, ien, ix, iy, iz, point(1:3)
    integer(8) :: ihilbert
    integer :: ixb, iyb, izb, bsize(1:3)
    type(block_t) :: blk
    integer :: idir, nn, order, size, bits
    FLOAT :: chi(1:MAX_DIM), xx(1:MAX_DIM)


    ...

  end subroutine create_x_lxyz
```

which 'serializes' the mesh by a choice of space-filling curve methods.


