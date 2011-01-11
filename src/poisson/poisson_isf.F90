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
!! $Id: poisson.F90 2544 2006-11-03 17:41:04Z xavier $

#include "global.h"

module poisson_isf_m
  use cube_function_m
  use datasets_m
  use global_m
  use messages_m
  use mesh_m
  use mpi_m
  use par_vec_m
  use parser_m
  use profiling_m
  
  implicit none
  
  private
  
  public ::               &
    poisson_isf_t,        &
    poisson_isf_init,     &
    poisson_isf_solve,    & 
    poisson_isf_end

  ! Datatype to store kernel values to solve Poisson equation
  ! on different communicators (configurations).
  type isf_cnf_t
    real(8), pointer  :: kernel(:, :, :)
    integer           :: nfft1, nfft2, nfft3
    type(mpi_grp_t)   :: mpi_grp
    logical           :: all_nodes
  end type isf_cnf_t

  type poisson_isf_t
    integer                      :: all_nodes_comm
    type(dcf_t)                  :: rho_cf
    type(isf_cnf_t)              :: cnf(1:4)
  end type poisson_isf_t

  ! Indices for the cnf array
  integer, parameter :: serial = 1
  integer, parameter :: world = 2
  integer, parameter :: domain = 3
  integer, parameter :: n_cnf = 3

  integer, parameter           :: order_scaling_function = 8 

  ! Interface to the Poisson solver calls
  interface

    subroutine build_kernel(n01,n02,n03,nfft1,nfft2,nfft3,hgrid,itype_scf,karrayout)
      integer, intent(in)       :: n01,n02,n03,nfft1,nfft2,nfft3,itype_scf
      real(kind=8), intent(in)  :: hgrid
      real(kind=8), dimension(nfft1/2+1,nfft2/2+1,nfft3/2+1), intent(out) :: karrayout
    end subroutine build_kernel

    subroutine calculate_dimensions(n01,n02,n03,nfft1,nfft2,nfft3)
      integer, intent(in)  :: n01,n02,n03
      integer, intent(out) :: nfft1,nfft2,nfft3
    end subroutine calculate_dimensions

    subroutine psolver_kernel(n01,n02,n03,nfft1,nfft2,nfft3, hgrid,karray,rhopot)
      integer, intent(in)    :: n01,n02,n03,nfft1,nfft2,nfft3
      real(8), intent(in)    :: hgrid
      real(8), intent(in), dimension(nfft1/2+1,nfft2/2+1,nfft3/2+1) :: karray
      real(8), intent(inout), dimension(n01,n02,n03) :: rhopot
    end subroutine psolver_kernel

  end interface

contains

  ! ---------------------------------------------------------
  subroutine poisson_isf_init(this, mesh, all_nodes_comm, init_world)
    type(poisson_isf_t), intent(out)   :: this
    type(mesh_t),        intent(in)    :: mesh
    integer,             intent(in)    :: all_nodes_comm
    logical, optional,   intent(in)    :: init_world 

    integer :: n1, n2, n3
#if defined(HAVE_MPI)
    integer :: m1, m2, m3, md1, md2, md3
    integer :: n(3)
    integer :: i_cnf
    logical :: init_world_
    integer :: default_nodes
    integer :: ierr, world_grp, poisson_grp, ii
    integer, allocatable :: ranks(:)
    !data ranks /0, 1/
    integer :: world_size
#endif
    integer :: nodes

    PUSH_SUB(poisson_isf_init)

#ifdef HAVE_MPI
    init_world_ = .true.
    if(present(init_world)) init_world_ = init_world
#endif

    call dcf_new(mesh%idx%ll, this%rho_cf)

    if(.not. mesh%parallel_in_domains) then
      ! The serial version is always needed (as used, e.g., in the casida runmode)
      call calculate_dimensions(this%rho_cf%n(1), this%rho_cf%n(2), this%rho_cf%n(3), &
        this%cnf(serial)%nfft1, this%cnf(serial)%nfft2, this%cnf(serial)%nfft3)

      n1 = this%cnf(serial)%nfft1/2 + 1
      n2 = this%cnf(serial)%nfft2/2 + 1
      n3 = this%cnf(serial)%nfft3/2 + 1

      SAFE_ALLOCATE(this%cnf(serial)%kernel(1:n1, 1:n2, 1:n3))

      call build_kernel(this%rho_cf%n(1), this%rho_cf%n(2), this%rho_cf%n(3),   &
        this%cnf(serial)%nfft1, this%cnf(serial)%nfft2, this%cnf(serial)%nfft3, &
        real(mesh%spacing(1), 8), order_scaling_function, this%cnf(serial)%kernel)
    end if
#if defined(HAVE_MPI)

    ! Allocate to configurations. The initialisation, especially the kernel,
    ! depends on the number of nodes used for the calculations. To avoid
    ! recalculating the kernel on each call of poisson_isf_solve depending on
    ! the all_nodes argument, both kernels are calculated.
    this%cnf(domain)%mpi_grp = mesh%mpi_grp

    ! For the world configuration we build a new communicator

    default_nodes = 0 !All nodes

    !%Variable PoissonSolverNodes
    !%Type integer
    !%Section Hamiltonian::Poisson
    !%Default 0
    !%Description
    !% How many nodes to use to solve the Poisson equation. A value of
    !% 0, the default, implies that all available nodes are used.
    !%End
    call parse_integer(datasets_check('PoissonSolverNodes'), default_nodes, nodes)

    this%all_nodes_comm = all_nodes_comm

    call MPI_Comm_size(all_nodes_comm, world_size, mpi_err)

    if(nodes <= 0 .or. nodes > world_size) nodes = mpi_world%size
    this%cnf(world)%all_nodes = (nodes == mpi_world%size)

    SAFE_ALLOCATE(ranks(1:nodes))

    do ii = 1, nodes
      ranks(ii) = ii - 1
    end do

    !create a new communicator
    !Extract the original group handle and create new comm.
    call MPI_Comm_group(all_nodes_comm, world_grp, ierr)
    call MPI_Group_incl(world_grp, nodes, ranks(1), poisson_grp, ierr)
    call MPI_Comm_create(mpi_world%comm, poisson_grp, this%cnf(world)%mpi_grp%comm, ierr)

    SAFE_DEALLOCATE_A(ranks)

    !Fill the new data structure, for all nodes
    if (this%cnf(world)%mpi_grp%comm /= MPI_COMM_NULL) then
      call MPI_Comm_rank(this%cnf(world)%mpi_grp%comm, this%cnf(world)%mpi_grp%rank, ierr)
      call MPI_Comm_size(this%cnf(world)%mpi_grp%comm, this%cnf(world)%mpi_grp%size, ierr)
    else
      this%cnf(world)%mpi_grp%rank = -1
      this%cnf(world)%mpi_grp%size = -1
    end if

    ! Build the kernel for all configurations. At the moment, this is
    ! solving the poisson equation with all nodes (i_cnf == world) and
    ! with the domain nodes only (i_cnf == domain).
    do i_cnf = 2, n_cnf
      if( (i_cnf == world .and. .not. init_world_) &                  ! world is disabled
        .or. (i_cnf == domain .and. .not. mesh%parallel_in_domains) & ! not parallel in domains
        ) then
        nullify(this%cnf(i_cnf)%kernel)
        cycle
      end if
      if (this%cnf(i_cnf)%mpi_grp%rank /= -1 .or. i_cnf /= world ) then
        call par_calculate_dimensions(this%rho_cf%n(1), this%rho_cf%n(2), this%rho_cf%n(3), &
          m1, m2, m3, n1, n2, n3, md1, md2, md3, this%cnf(i_cnf)%nfft1, this%cnf(i_cnf)%nfft2, &
          this%cnf(i_cnf)%nfft3, this%cnf(i_cnf)%mpi_grp%size)

        ! Shortcuts to avoid to "line too long" errors.
        n(1) = this%cnf(i_cnf)%nfft1
        n(2) = this%cnf(i_cnf)%nfft2
        n(3) = this%cnf(i_cnf)%nfft3

        SAFE_ALLOCATE(this%cnf(i_cnf)%kernel(1:n(1), 1:n(2), 1:n(3)/this%cnf(i_cnf)%mpi_grp%size))

        call par_build_kernel(this%rho_cf%n(1), this%rho_cf%n(2), this%rho_cf%n(3), n1, n2, n3,     &
          this%cnf(i_cnf)%nfft1, this%cnf(i_cnf)%nfft2, this%cnf(i_cnf)%nfft3,                      &
          mesh%spacing(1), order_scaling_function,                                            &
          this%cnf(i_cnf)%mpi_grp%rank, this%cnf(i_cnf)%mpi_grp%size, this%cnf(i_cnf)%mpi_grp%comm, &
          this%cnf(i_cnf)%kernel)
      else
      	nullify(this%cnf(i_cnf)%kernel)
        cycle
      end if
    end do
#endif

    POP_SUB(poisson_isf_init)
  end subroutine poisson_isf_init

  ! ---------------------------------------------------------
  subroutine poisson_isf_solve(this, mesh, pot, rho, all_nodes)
    type(poisson_isf_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(out)   :: pot(:)
    FLOAT,               intent(in)    :: rho(:)
    logical,             intent(in)    :: all_nodes

    integer :: i_cnf
#if defined(HAVE_MPI)
    FLOAT, allocatable :: rho_global(:)
    FLOAT, allocatable :: pot_global(:)
#endif
    real(8), pointer :: rhop(:,:,:)
    
    !variables to meassure the global comunications time
    real :: t0,t1,tel
    !integer :: count1,count2,count_rate,count_max
    !integer :: sec1,sec2,usec1,usec2,sec1_t,sec2_t,usec1_t,usec2_t

    PUSH_SUB(poisson_isf_solve)
    !call loct_gettimeofday(sec1_t,usec1_t)

    call dcf_alloc_RS(this%rho_cf)

    if(mesh%parallel_in_domains) then
    
#if defined(HAVE_MPI)
      SAFE_ALLOCATE(rho_global(1:mesh%np_global))
      SAFE_ALLOCATE(pot_global(1:mesh%np_global))

      ! At this point, dvec_allgather is required because the ISF solver
      ! uses another data distribution algorithm than for the mesh functions
      ! for which every node requires the full data.

      !call loct_gettimeofday(sec1,usec1)
      call dvec_allgather(mesh%vp, rho_global, rho)
      !call loct_gettimeofday(sec2,usec2)
      !call time_diff(sec1,usec1,sec2,usec2)
     
      !write(78,*) 'PSolver: 1st DVEC_ALLGATHER TIME',sec2,'us',usec2
     
      call dmesh_to_cube(mesh, rho_global, this%rho_cf)
#endif
    else
      call dmesh_to_cube(mesh, rho, this%rho_cf)
    end if

    ! Choose configuration.
    i_cnf = serial

#if defined(HAVE_MPI)
    if(all_nodes) then
      i_cnf = world
    else if(mesh%parallel_in_domains) then
      i_cnf = domain
    end if
#endif

#if !defined(HAVE_MPI)
    ASSERT(i_cnf == serial)
#endif

    if(i_cnf == serial) then

#ifdef SINGLE_PRECISION
      SAFE_ALLOCATE(rhop(1:this%rho_cf%n(1), 1:this%rho_cf%n(2), 1:this%rho_cf%n(3)))
      rhop = this%rho_cf%RS
#else
      rhop => this%rho_cf%RS
#endif

      call psolver_kernel(this%rho_cf%n(1), this%rho_cf%n(2), this%rho_cf%n(3),    &
        this%cnf(serial)%nfft1, this%cnf(serial)%nfft2, this%cnf(serial)%nfft3, &
        real(mesh%spacing(1), 8), this%cnf(serial)%kernel, rhop)

#ifdef SINGLE_PRECISION
      this%rho_cf%RS = rhop
      SAFE_DEALLOCATE_P(rhop)
#endif

#if defined(HAVE_MPI)
    else
      if (this%cnf(i_cnf)%mpi_grp%size /= -1 .or. i_cnf /= world) then
        call par_psolver_kernel(this%rho_cf%n(1), this%rho_cf%n(2), this%rho_cf%n(3), &
          this%cnf(i_cnf)%nfft1, this%cnf(i_cnf)%nfft2, this%cnf(i_cnf)%nfft3,  &
          real(mesh%spacing(1), 8), this%cnf(i_cnf)%kernel, this%rho_cf%RS,                      &
          this%cnf(i_cnf)%mpi_grp%rank, this%cnf(i_cnf)%mpi_grp%size, this%cnf(i_cnf)%mpi_grp%comm)
      end if
      ! we need to be sure that the root of every domain-partition has a copy of the potential
      ! for the moment we broadcast to all nodes, but this is more than what we really need 
      if(i_cnf == world .and. .not. this%cnf(world)%all_nodes) then
        call MPI_Bcast(this%rho_cf%rs(1, 1, 1), this%rho_cf%n(1)*this%rho_cf%n(2)*this%rho_cf%n(3), &
          MPI_FLOAT, 0, this%all_nodes_comm, mpi_err)
      end if
#endif
    end if

    if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
       call dcube_to_mesh(mesh, this%rho_cf, pot_global)

       !call loct_gettimeofday(sec1,usec1)
       call dvec_scatter(mesh%vp, mesh%vp%root, pot_global, pot)
       !call loct_gettimeofday(sec2,usec2)
       !call time_diff(sec1,usec1,sec2,usec2)
       
       !write(78,*) 'PSolver: 1st DVEC_SCATTER TIME',sec2,'us',usec2

       SAFE_DEALLOCATE_A(rho_global)
       SAFE_DEALLOCATE_A(pot_global)
#endif
    else
       call dcube_to_mesh(mesh, this%rho_cf, pot)
    end if
    
    call dcf_free_RS(this%rho_cf)
    
    !call loct_gettimeofday(sec2_t,usec2_t)
    !call time_diff(sec1_t,usec1_t,sec2_t,usec2_t)
    !write(78,*) 'PoissonTime: iteration POISSON TIME',sec2_t,'us',usec2_t
    POP_SUB(poisson_isf_solve)
  end subroutine poisson_isf_solve

  ! ---------------------------------------------------------
  subroutine poisson_isf_end(this)
    type(poisson_isf_t), intent(inout) :: this

#if defined(HAVE_MPI)
    integer :: i_cnf
#endif

    PUSH_SUB(poisson_isf_end)

#if defined(HAVE_MPI)
    do i_cnf = 1, n_cnf
      SAFE_DEALLOCATE_P(this%cnf(i_cnf)%kernel)
    end do
    if(this%cnf(world)%mpi_grp%comm /= MPI_COMM_NULL) then
      call MPI_Comm_free(this%cnf(world)%mpi_grp%comm, mpi_err)
    end if
#else
    SAFE_DEALLOCATE_P(this%cnf(serial)%kernel)
#endif

    call dcf_free(this%rho_cf)

    POP_SUB(poisson_isf_end)
  end subroutine poisson_isf_end

end module poisson_isf_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
