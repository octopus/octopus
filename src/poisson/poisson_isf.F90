!! Copyright (C) 2002-2011 M. Marques, A. Castro, X. Andrade, J. Alberdi, M. Oliveira
!! Copyright (C) Luigi Genovese, Thierry Deutsch, CEA Grenoble, 2006
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
  use cube_m
  use datasets_m
  use global_m
  use io_m
  use messages_m
  use mesh_m
  use mpi_m
  use par_vec_m
  use parser_m
  use profiling_m
  use scaling_function_m
  use sgfft_m

  implicit none
  
  private
  
  public ::               &
    poisson_isf_t,        &
    isf_cnf_t,    &
    poisson_isf_init,     &
    poisson_isf_solve,    & 
    poisson_isf_end

  ! Indices for the cnf array
  integer, parameter :: SERIAL = 1
  integer, parameter :: WORLD = 2
  integer, parameter :: DOMAIN = 3
  integer, parameter :: N_CNF = 3

  ! Datatype to store kernel values to solve Poisson equation
  ! on different communicators (configurations).
  type isf_cnf_t
    real(8), pointer  :: kernel(:, :, :)
    integer           :: nfft1, nfft2, nfft3
    type(mpi_grp_t)   :: mpi_grp
    logical           :: all_nodes
  end type isf_cnf_t

  type poisson_isf_t
    integer         :: all_nodes_comm
    type(isf_cnf_t) :: cnf(1:N_CNF)
  end type poisson_isf_t

  integer, parameter :: order_scaling_function = 8 

contains

  ! ---------------------------------------------------------
  subroutine poisson_isf_init(this, mesh, cube, all_nodes_comm, init_world)
    type(poisson_isf_t), intent(out)   :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube
    integer,             intent(in)    :: all_nodes_comm
    logical, optional,   intent(in)    :: init_world 

    integer :: n1, n2, n3
    integer :: i_cnf
#if defined(HAVE_MPI)
    integer :: m1, m2, m3, md1, md2, md3
    integer :: n(3)
    logical :: init_world_
    integer :: default_nodes
    integer :: ierr, world_grp, poisson_grp, ii
    integer, allocatable :: ranks(:)
    !data ranks /0, 1/
    integer :: world_size
    integer :: nodes
#endif

    PUSH_SUB(poisson_isf_init)

#ifdef HAVE_MPI
    init_world_ = .true.
    if(present(init_world)) init_world_ = init_world
#endif

    ! we need to nullify the pointer so they can be deallocated safely
    ! afterwards
    do i_cnf = 1, N_CNF
      nullify(this%cnf(i_cnf)%kernel)
    end do

    if(.not. mesh%parallel_in_domains) then
      ! The serial version is always needed (as used, e.g., in the casida runmode)
      call calculate_dimensions(cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), &
        this%cnf(SERIAL)%nfft1, this%cnf(SERIAL)%nfft2, this%cnf(SERIAL)%nfft3)

      n1 = this%cnf(SERIAL)%nfft1/2 + 1
      n2 = this%cnf(SERIAL)%nfft2/2 + 1
      n3 = this%cnf(SERIAL)%nfft3/2 + 1

      SAFE_ALLOCATE(this%cnf(SERIAL)%kernel(1:n1, 1:n2, 1:n3))

      call build_kernel(cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3),   &
        this%cnf(SERIAL)%nfft1, this%cnf(SERIAL)%nfft2, this%cnf(SERIAL)%nfft3, &
        real(mesh%spacing(1), 8), order_scaling_function, this%cnf(SERIAL)%kernel)
    end if
#if defined(HAVE_MPI)

    ! Allocate to configurations. The initialisation, especially the kernel,
    ! depends on the number of nodes used for the calculations. To avoid
    ! recalculating the kernel on each call of poisson_isf_solve depending on
    ! the all_nodes argument, both kernels are calculated.
    this%cnf(DOMAIN)%mpi_grp = mesh%mpi_grp

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

    if(nodes <= 0 .or. nodes > world_size) nodes = world_size
    this%cnf(WORLD)%all_nodes = (nodes == mpi_world%size)

    SAFE_ALLOCATE(ranks(1:nodes))

    do ii = 1, nodes
      ranks(ii) = ii - 1
    end do

    !create a new communicator
    !Extract the original group handle and create new comm.
    call MPI_Comm_group(all_nodes_comm, world_grp, ierr)
    call MPI_Group_incl(world_grp, nodes, ranks(1), poisson_grp, ierr)
    call MPI_Comm_create(mpi_world%comm, poisson_grp, this%cnf(WORLD)%mpi_grp%comm, ierr)

    SAFE_DEALLOCATE_A(ranks)

    !Fill the new data structure, for all nodes
    if (this%cnf(WORLD)%mpi_grp%comm /= MPI_COMM_NULL) then
      call MPI_Comm_rank(this%cnf(WORLD)%mpi_grp%comm, this%cnf(WORLD)%mpi_grp%rank, ierr)
      call MPI_Comm_size(this%cnf(WORLD)%mpi_grp%comm, this%cnf(WORLD)%mpi_grp%size, ierr)
    else
      this%cnf(WORLD)%mpi_grp%rank = -1
      this%cnf(WORLD)%mpi_grp%size = -1
    end if

    ! Build the kernel for all configurations. At the moment, this is
    ! solving the poisson equation with all nodes (i_cnf == WORLD) and
    ! with the domain nodes only (i_cnf == DOMAIN).
    do i_cnf = 2, N_CNF
      if( (i_cnf == WORLD .and. .not. init_world_) &                  ! world is disabled
        .or. (i_cnf == DOMAIN .and. .not. mesh%parallel_in_domains) & ! not parallel in domains
        ) then
        nullify(this%cnf(i_cnf)%kernel)
        cycle
      end if
      if (this%cnf(i_cnf)%mpi_grp%rank /= -1 .or. i_cnf /= WORLD ) then
        call par_calculate_dimensions(cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), &
          m1, m2, m3, n1, n2, n3, md1, md2, md3, this%cnf(i_cnf)%nfft1, this%cnf(i_cnf)%nfft2, &
          this%cnf(i_cnf)%nfft3, this%cnf(i_cnf)%mpi_grp%size)

        ! Shortcuts to avoid to "line too long" errors.
        n(1) = this%cnf(i_cnf)%nfft1
        n(2) = this%cnf(i_cnf)%nfft2
        n(3) = this%cnf(i_cnf)%nfft3

        SAFE_ALLOCATE(this%cnf(i_cnf)%kernel(1:n(1), 1:n(2), 1:n(3)/this%cnf(i_cnf)%mpi_grp%size))

        call par_build_kernel(cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), n1, n2, n3,     &
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
  subroutine poisson_isf_solve(this, mesh, cube, pot, rho, all_nodes)
    type(poisson_isf_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(in)    :: cube
    FLOAT,               intent(out)   :: pot(:)
    FLOAT,               intent(in)    :: rho(:)
    logical,             intent(in)    :: all_nodes

    integer :: i_cnf
    real(8), pointer :: rhop(:,:,:)
    type(cube_function_t) :: rho_cf
    
    PUSH_SUB(poisson_isf_solve)

    call cube_function_null(rho_cf)
    call dcube_function_alloc_RS(cube, rho_cf)

    if(mesh%parallel_in_domains) then
      call dmesh_to_cube(mesh, rho, cube, rho_cf, local=.true.)
    else
      call dmesh_to_cube(mesh, rho, cube, rho_cf)
    end if

    ! Choose configuration.
    i_cnf = SERIAL

#if defined(HAVE_MPI)
    if(all_nodes) then
      i_cnf = WORLD
    else if(mesh%parallel_in_domains) then
      i_cnf = DOMAIN
    end if
#endif

#if !defined(HAVE_MPI)
    ASSERT(i_cnf == SERIAL)
#endif

    if(i_cnf == SERIAL) then

#ifdef SINGLE_PRECISION
      SAFE_ALLOCATE(rhop(1:cube%rs_n_global(1), 1:cube%rs_n_global(2), 1:cube%rs_n_global(3)))
      rhop = rho_cf%dRS
#else
      rhop => rho_cf%dRS
#endif

      call psolver_kernel(cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3),    &
        this%cnf(SERIAL)%nfft1, this%cnf(SERIAL)%nfft2, this%cnf(SERIAL)%nfft3, &
        real(mesh%spacing(1), 8), this%cnf(SERIAL)%kernel, rhop)

#ifdef SINGLE_PRECISION
      rho_cf%dRS = rhop
      SAFE_DEALLOCATE_P(rhop)
#endif

#if defined(HAVE_MPI)
    else
      if (this%cnf(i_cnf)%mpi_grp%size /= -1 .or. i_cnf /= WORLD) then
        call par_psolver_kernel(cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), &
          this%cnf(i_cnf)%nfft1, this%cnf(i_cnf)%nfft2, this%cnf(i_cnf)%nfft3,  &
          real(mesh%spacing(1), 8), this%cnf(i_cnf)%kernel, rho_cf%dRS,                      &
          this%cnf(i_cnf)%mpi_grp%rank, this%cnf(i_cnf)%mpi_grp%size, this%cnf(i_cnf)%mpi_grp%comm)
      end if
      ! we need to be sure that the root of every domain-partition has a copy of the potential
      ! for the moment we broadcast to all nodes, but this is more than what we really need 
      if(i_cnf == WORLD .and. .not. this%cnf(WORLD)%all_nodes) then
        call MPI_Bcast(rho_cf%drs(1, 1, 1), cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3), &
          MPI_FLOAT, 0, this%all_nodes_comm, mpi_err)
      end if
#endif
    end if

    if(mesh%parallel_in_domains) then
      call dcube_to_mesh(cube, rho_cf, mesh, pot, local=.true.)
    else
       call dcube_to_mesh(cube, rho_cf, mesh, pot)
    end if
    
    call dcube_function_free_RS(cube, rho_cf)

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
    do i_cnf = 1, N_CNF
      SAFE_DEALLOCATE_P(this%cnf(i_cnf)%kernel)
    end do
    if(this%cnf(WORLD)%mpi_grp%comm /= MPI_COMM_NULL) then
      call MPI_Comm_free(this%cnf(WORLD)%mpi_grp%comm, mpi_err)
    end if
#else
    SAFE_DEALLOCATE_P(this%cnf(SERIAL)%kernel)
#endif

    POP_SUB(poisson_isf_end)
  end subroutine poisson_isf_end

  ! --------------------------------------------------------------

  !!****h* BigDFT/psolver_kernel
  !! NAME
  !!   psolver_kernel
  !!
  !! FUNCTION
  !!    Solver of Poisson equation applying a kernel
  !!
  !! SYNOPSIS
  !!    Poisson solver applying a kernel and 
  !!    using Fourier transform for the convolution.
  !!    rhopot : input  -> the density
  !!             output -> the Hartree potential + pot_ion
  !!    The potential pot_ion is ADDED in the array rhopot.
  !!    Calculate also the Hartree potential
  !!
  !!    Replaces the charge density contained in rhopot 
  !!    by the Hartree stored as well in rhopot.
  !!    If xc_on is true, it also adds the XC potential and 
  !!    ionic potential pot_ion
  !!
  !!    We double the size of the mesh except in one dimension
  !!    in order to use the property of the density to be real.
  !! WARNING
  !!    For the use of FFT routine
  !!        inzee=1: first part of Z is data (output) array, 
  !!                 second part work array
  !!        inzee=2: first part of Z is work array, second part data array
  !!                 real(F(i1,i2,i3))=Z(1,i1,i2,i3,inzee)
  !!                 imag(F(i1,i2,i3))=Z(2,i1,i2,i3,inzee)
  !!        inzee on output is in general different from inzee on input
  !!
  !! AUTHOR
  !!    Thierry Deutsch, Luigi Genovese
  !! COPYRIGHT
  !!    Copyright (C) 2005 CEA
  !! CREATION DATE
  !!    13/07/2005
  !!
  !! MODIFICATION HISTORY
  !!    12/2005 Kernel stored into memory
  !!    12/2005 Real Kernel FFT and use less memory
  !!
  !! SOURCE
  !!
  subroutine psolver_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, karray, rhopot)
    integer, intent(in)    :: n01
    integer, intent(in)    :: n02
    integer, intent(in)    :: n03
    integer, intent(in)    :: nfft1
    integer, intent(in)    :: nfft2
    integer, intent(in)    :: nfft3
    real(8), intent(in)    :: hgrid
    real(8), intent(in)    :: karray(nfft1/2 + 1,nfft2/2 + 1, nfft3/2 + 1)
    real(8), intent(inout) :: rhopot(n01, n02, n03)

    real(8), allocatable :: zarray(:,:,:)
    real(8) :: factor
    integer :: n1, n2, n3, nd1, nd2, nd3, n1h, nd1h
    integer :: inzee, i_sign

    !Dimension of the FFT
    call calculate_dimensions(n01, n02, n03, n1, n2, n3)

    !Half size of nd1
    n1h=n1/2
    nd1 = n1 + modulo(n1+1,2)
    nd2 = n2 + modulo(n2+1,2)
    nd3 = n3 + modulo(n3+1,2)
    nd1h=(nd1+1)/2

    SAFE_ALLOCATE(zarray(1:2, 1:nd1h*nd2*nd3, 1:2))

    !Set zarray
    call zarray_in(n01,n02,n03,nd1h,nd2,nd3,rhopot,zarray)

    !FFT
    !print *,"Do a 3D HalFFT for the density"
    i_sign=1
    inzee=1
    call fft(n1h,n2,n3,nd1h,nd2,nd3,zarray,i_sign,inzee)

    !print *, "Apply the kernel"
    call kernel_application(n1,n2,n3,nd1h,nd2,nd3,nfft1,nfft2,nfft3,zarray,karray,inzee)

    !Inverse FFT
    i_sign=-1
    !print *,"Do a 3D inverse HalFFT"
    call fft(n1h,n2,n3,nd1h,nd2,nd3,zarray,i_sign,inzee)

    !Recollect the result
    !We have to multiply by a factor
    factor = hgrid**3/(n1*n2*n3)

    ! Calling this routine gives only the Hartree potential
    call zarray_out(n01, n02, n03, nd1h, nd2, nd3, rhopot, zarray(1, 1, inzee), factor)

    SAFE_DEALLOCATE_A(zarray)
  end subroutine psolver_kernel
  !!***

  !!****h* BigDFT/kernel_application
  !! NAME
  !!   kernel_application
  !!
  !! FUNCTION
  !!    Multiply the FFT of the density by the FFT of the kernel
  !!
  !! SYNOPSIS
  !!    zarray(:,:,:,:,inzee) : IN -> FFT of the density with the x dimension divided by two
  !!                            (HalFFT), OUT -> FFT of the potential
  !!    karray                : kernel FFT (real, 1/8 of the total grid)
  !!    n1h,n2,n3             : dimension of the FFT grid for zarray
  !!    nd1h,nd2,nd3          : dimensions of zarray
  !!    nfft1,nfft2,nfft3     : original FFT grid dimensions, to be used for karray dimensions
  !!
  !! WARNING
  !!    We use all the zarray vector, storing the auxiliary part using ouzee=3-inzee
  !!    All the loop are unrolled such to avoid different conditions
  !!    the "min" functions are substituted by kink computed with absolute values
  !! AUTHOR
  !!    Luigi Genovese
  !! CREATION DATE
  !!    March 2006
  !!
  !! SOURCE
  !!
  subroutine kernel_application(n1,n2,n3,nd1h,nd2,nd3,nfft1,nfft2,nfft3,zarray,karray,inzee)
    integer, intent(in)    :: n1,n2,n3,nd1h,nd2,nd3,nfft1,nfft2,nfft3,inzee
    real(8), intent(in)    :: karray(1:nfft1/2 + 1, 1:nfft2/2 + 1, 1:nfft3/2 + 1)
    real(8), intent(inout) :: zarray(1:2, 1:nd1h, 1:nd2, 1:nd3, 1:2)

    real(8), dimension(:), allocatable :: cos_array,sin_array
    real(8) :: a,b,c,d,pi2,g1,cp,sp
    real(8) :: rfe,ife,rfo,ifo,rk,ik,rk2,ik2,re,ro,ie,io,rhk,ihk
    integer :: i1,i2,i3,j1,j2,j3, ouzee,n1h,n2h,n3h
    integer :: si1,si2,si3

    n1h = n1/2
    n2h = n2/2
    n3h = n3/2

    SAFE_ALLOCATE(cos_array(1:n1h + 1))
    SAFE_ALLOCATE(sin_array(1:n1h + 1))

    pi2=8.d0*datan(1.d0)
    pi2=pi2/real(n1,8)
    do i1=1,n1h+1
      cos_array(i1)=dcos(pi2*(i1-1))
      sin_array(i1)=-dsin(pi2*(i1-1))
    end do

    ouzee=3-inzee

    !--------------------------------------------!
    !--- Starting reconstruction half -> full ---!
    !--------------------------------------------!   

    !-------------Case i3 = 1
    i3=1
    j3=1
    si3=1

    !-------------Case i2 = 1, i3 = 1
    i2=1
    j2=1
    si2=1

    !Case i1 == 1
    i1=1
    si1=1
    a=zarray(1,i1,i2,i3,inzee)
    b=zarray(2,i1,i2,i3,inzee)
    c=zarray(1,si1,si2,si3,inzee)
    d=zarray(2,si1,si2,si3,inzee)
    rfe=.5d0*(a+c)
    ife=.5d0*(b-d)
    rfo=.5d0*(a-c)
    ifo=.5d0*(b+d) 
    cp=cos_array(i1)
    sp=sin_array(i1)
    rk=rfe+cp*ifo-sp*rfo
    ik=ife-cp*rfo-sp*ifo
    g1=karray(i1,j2,j3)
    rk2=rk*g1
    ik2=ik*g1

    zarray(1,1,i2,i3,ouzee) = rk2
    zarray(2,1,i2,i3,ouzee) = ik2

    !Case i1=2,n1h
    do i1=2,n1h
      si1=n1h+2-i1

      a=zarray(1,i1,i2,i3,inzee)
      b=zarray(2,i1,i2,i3,inzee)
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,i1,i2,i3,ouzee) = rk2
      zarray(2,i1,i2,i3,ouzee) = ik2
    end do

    !Case i1=n1h+1
    i1=n1h+1
    si1=n1h+2-i1

    a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
    b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
    c=zarray(1,si1,si2,si3,inzee)
    d=zarray(2,si1,si2,si3,inzee)
    rfe=.5d0*(a+c)
    ife=.5d0*(b-d)
    rfo=.5d0*(a-c)
    ifo=.5d0*(b+d) 
    cp=cos_array(i1)
    sp=sin_array(i1)
    rk=rfe+cp*ifo-sp*rfo
    ik=ife-cp*rfo-sp*ifo
    g1=karray(i1,j2,j3)
    rk2=rk*g1
    ik2=ik*g1

    zarray(1,i1,i2,i3,ouzee) = rk2
    zarray(2,i1,i2,i3,ouzee) = ik2
    !-------------END case i2 = 1 , i3=1

    !case i2 >=2
    do i2=2,n2
      j2=n2h+1-abs(n2h+1-i2)
      si2=n2+2-i2 !if i2 /=1, otherwise si2=1

      !Case i1 == 1
      i1=1
      si1=1
      a=zarray(1,i1,i2,i3,inzee)
      b=zarray(2,i1,i2,i3,inzee)
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,1,i2,i3,ouzee) = rk2
      zarray(2,1,i2,i3,ouzee) = ik2

      !Case i1=2,n1h
      do i1=2,n1h
        si1=n1h+2-i1

        a=zarray(1,i1,i2,i3,inzee)
        b=zarray(2,i1,i2,i3,inzee)
        c=zarray(1,si1,si2,si3,inzee)
        d=zarray(2,si1,si2,si3,inzee)
        rfe=.5d0*(a+c)
        ife=.5d0*(b-d)
        rfo=.5d0*(a-c)
        ifo=.5d0*(b+d) 
        cp=cos_array(i1)
        sp=sin_array(i1)
        rk=rfe+cp*ifo-sp*rfo
        ik=ife-cp*rfo-sp*ifo
        g1=karray(i1,j2,j3)
        rk2=rk*g1
        ik2=ik*g1

        zarray(1,i1,i2,i3,ouzee) = rk2
        zarray(2,i1,i2,i3,ouzee) = ik2
      end do

      !Case i1=n1h+1
      i1=n1h+1
      si1=n1h+2-i1

      a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
      b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,i1,i2,i3,ouzee) = rk2
      zarray(2,i1,i2,i3,ouzee) = ik2
    end do
    !-------------END Case i3 = 1

    !case i3 >=2
    do i3=2,n3
      j3=n3h+1-abs(n3h+1-i3)
      si3=n3+2-i3 !if i3 /=1, otherwise si3=1

      !-------------Case i2 = 1
      i2=1
      j2=1
      si2=1

      !Case i1 == 1
      i1=1
      si1=1
      a=zarray(1,i1,i2,i3,inzee)
      b=zarray(2,i1,i2,i3,inzee)
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,1,i2,i3,ouzee) = rk2
      zarray(2,1,i2,i3,ouzee) = ik2

      !Case i1=2,n1h
      do i1=2,n1h
        si1=n1h+2-i1

        a=zarray(1,i1,i2,i3,inzee)
        b=zarray(2,i1,i2,i3,inzee)
        c=zarray(1,si1,si2,si3,inzee)
        d=zarray(2,si1,si2,si3,inzee)
        rfe=.5d0*(a+c)
        ife=.5d0*(b-d)
        rfo=.5d0*(a-c)
        ifo=.5d0*(b+d) 
        cp=cos_array(i1)
        sp=sin_array(i1)
        rk=rfe+cp*ifo-sp*rfo
        ik=ife-cp*rfo-sp*ifo
        g1=karray(i1,j2,j3)
        rk2=rk*g1
        ik2=ik*g1

        zarray(1,i1,i2,i3,ouzee) = rk2
        zarray(2,i1,i2,i3,ouzee) = ik2
      end do

      !Case i1=n1h+1
      i1=n1h+1
      si1=n1h+2-i1

      a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
      b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,i1,i2,i3,ouzee) = rk2
      zarray(2,i1,i2,i3,ouzee) = ik2
      !-------------END case i2 = 1      

      !case i2 >=2
      do i2=2,n2
        j2=n2h+1-abs(n2h+1-i2)
        si2=n2+2-i2 !if i2 /=1, otherwise si2=1

        !Case i1 == 1
        i1=1
        si1=1
        a=zarray(1,i1,i2,i3,inzee)
        b=zarray(2,i1,i2,i3,inzee)
        c=zarray(1,si1,si2,si3,inzee)
        d=zarray(2,si1,si2,si3,inzee)
        rfe=.5d0*(a+c)
        ife=.5d0*(b-d)
        rfo=.5d0*(a-c)
        ifo=.5d0*(b+d) 
        cp=cos_array(i1)
        sp=sin_array(i1)
        rk=rfe+cp*ifo-sp*rfo
        ik=ife-cp*rfo-sp*ifo
        g1=karray(i1,j2,j3)
        rk2=rk*g1
        ik2=ik*g1

        zarray(1,1,i2,i3,ouzee) = rk2
        zarray(2,1,i2,i3,ouzee) = ik2

        !Case i1=2,n1h
        do i1=2,n1h
          si1=n1h+2-i1

          a=zarray(1,i1,i2,i3,inzee)
          b=zarray(2,i1,i2,i3,inzee)
          c=zarray(1,si1,si2,si3,inzee)
          d=zarray(2,si1,si2,si3,inzee)
          rfe=.5d0*(a+c)
          ife=.5d0*(b-d)
          rfo=.5d0*(a-c)
          ifo=.5d0*(b+d) 
          cp=cos_array(i1)
          sp=sin_array(i1)
          rk=rfe+cp*ifo-sp*rfo
          ik=ife-cp*rfo-sp*ifo
          g1=karray(i1,j2,j3)
          rk2=rk*g1
          ik2=ik*g1

          zarray(1,i1,i2,i3,ouzee) = rk2
          zarray(2,i1,i2,i3,ouzee) = ik2
        end do

        !Case i1=n1h+1
        i1=n1h+1
        si1=n1h+2-i1

        a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
        b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
        c=zarray(1,si1,si2,si3,inzee)
        d=zarray(2,si1,si2,si3,inzee)
        rfe=.5d0*(a+c)
        ife=.5d0*(b-d)
        rfo=.5d0*(a-c)
        ifo=.5d0*(b+d) 
        cp=cos_array(i1)
        sp=sin_array(i1)
        rk=rfe+cp*ifo-sp*rfo
        ik=ife-cp*rfo-sp*ifo
        g1=karray(i1,j2,j3)
        rk2=rk*g1
        ik2=ik*g1

        zarray(1,i1,i2,i3,ouzee) = rk2
        zarray(2,i1,i2,i3,ouzee) = ik2
      end do

    end do


    !--------------------------------------------!
    !--- Starting reconstruction full -> half ---!
    !--------------------------------------------!   

    !case i3=1
    i3=1
    j3=1
    !case i2=1
    i2=1
    j2=1
    do i1=1,n1h
      j1=n1h+2-i1

      a=zarray(1,i1,i2,i3,ouzee)
      b=zarray(2,i1,i2,i3,ouzee)
      c=zarray(1,j1,j2,j3,ouzee)
      d=-zarray(2,j1,j2,j3,ouzee)
      cp=cos_array(i1)
      sp=sin_array(i1)
      re=(a+c)
      ie=(b+d)
      ro=(a-c)*cp-(b-d)*sp
      io=(a-c)*sp+(b-d)*cp
      rhk=re-io 
      ihk=ie+ro

      zarray(1,i1,i2,i3,inzee)=rhk
      zarray(2,i1,i2,i3,inzee)=ihk
    end do
    !case i2 >= 2
    do i2=2,n2
      j2=nd2+1-i2
      do i1=1,n1h
        j1=n1h+2-i1

        a=zarray(1,i1,i2,i3,ouzee)
        b=zarray(2,i1,i2,i3,ouzee)
        c=zarray(1,j1,j2,j3,ouzee)
        d=-zarray(2,j1,j2,j3,ouzee)
        cp=cos_array(i1)
        sp=sin_array(i1)
        re=(a+c)
        ie=(b+d)
        ro=(a-c)*cp-(b-d)*sp
        io=(a-c)*sp+(b-d)*cp
        rhk=re-io 
        ihk=ie+ro

        zarray(1,i1,i2,i3,inzee)=rhk
        zarray(2,i1,i2,i3,inzee)=ihk
      end do
    end do


    !case i3 >=2
    do i3=2,n3
      j3=nd3+1-i3
      !case i2=1
      i2=1
      j2=1
      do i1=1,n1h
        j1=n1h+2-i1

        a=zarray(1,i1,i2,i3,ouzee)
        b=zarray(2,i1,i2,i3,ouzee)
        c=zarray(1,j1,j2,j3,ouzee)
        d=-zarray(2,j1,j2,j3,ouzee)
        cp=cos_array(i1)
        sp=sin_array(i1)
        re=(a+c)
        ie=(b+d)
        ro=(a-c)*cp-(b-d)*sp
        io=(a-c)*sp+(b-d)*cp
        rhk=re-io 
        ihk=ie+ro

        zarray(1,i1,i2,i3,inzee)=rhk
        zarray(2,i1,i2,i3,inzee)=ihk
      end do
      !case i2 >= 2
      do i2=2,n2
        j2=nd2+1-i2
        do i1=1,n1h
          j1=n1h+2-i1

          a=zarray(1,i1,i2,i3,ouzee)
          b=zarray(2,i1,i2,i3,ouzee)
          c=zarray(1,j1,j2,j3,ouzee)
          d=-zarray(2,j1,j2,j3,ouzee)
          cp=cos_array(i1)
          sp=sin_array(i1)
          re=(a+c)
          ie=(b+d)
          ro=(a-c)*cp-(b-d)*sp
          io=(a-c)*sp+(b-d)*cp
          rhk=re-io 
          ihk=ie+ro

          zarray(1,i1,i2,i3,inzee)=rhk
          zarray(2,i1,i2,i3,inzee)=ihk
        end do
      end do

    end do

    !De-allocations
    SAFE_DEALLOCATE_A(cos_array)
    SAFE_DEALLOCATE_A(sin_array)

  end subroutine kernel_application

  !!****h* BigDFT/norm_ind
  !! NAME
  !!   norm_ind
  !!
  !! FUNCTION
  !!   Index in zarray
  !!
  !! SOURCE
  !!
  subroutine norm_ind(nd1,nd2,nd3,i1,i2,i3,ind)
    integer :: nd1,nd2,nd3,i1,i2,i3
    integer :: ind

    !Local variables
    integer :: a1,a2,a3
    if ( i1 == nd1 ) then
      a1=1
    else
      a1=i1
    end if
    if ( i2 == nd2 ) then
      a2=1
    else
      a2=i2
    end if
    if ( i3 == nd3 ) then
      a3=1
    else
      a3=i3
    end if
    ind=a1+nd1*(a2-1)+nd1*nd2*(a3-1)
  end subroutine norm_ind
  !!***


  !!****h* BigDFT/symm_ind
  !! NAME
  !!   symm_ind
  !!
  !! FUNCTION
  !!   Index in zarray for -g vector
  !!
  !! SOURCE
  !!
  subroutine symm_ind(nd1,nd2,nd3,i1,i2,i3,ind)
    integer :: nd1,nd2,nd3,i1,i2,i3
    integer :: ind

    integer ::  a1,a2,a3
    if (i1 /= 1) then 
      a1=nd1+1-i1
    else
      a1=i1
    end if
    if (i2 /= 1) then 
      a2=nd2+1-i2
    else
      a2=i2
    end if
    if (i3 /= 1) then 
      a3=nd3+1-i3
    else
      a3=i3
    end if
    ind=a1+nd1*(a2-1)+nd1*nd2*(a3-1)
  end subroutine symm_ind
  !!***
  
  !!****h* BigDFT/zarray_in
  !! NAME
  !!   zarray_in
  !!
  !! FUNCTION
  !!   Put the density into zarray
  !!
  !! SOURCE
  !!
  subroutine zarray_in(n01,n02,n03,nd1,nd2,nd3,density,zarray)
    integer :: n01,n02,n03,nd1,nd2,nd3
    real(8), dimension(n01,n02,n03) :: density
    real(8), dimension(2,nd1,nd2,nd3) :: zarray

    integer :: i1,i2,i3,n01h,nd1hm,nd3hm,nd2hm

    !Half the size of n01
    n01h=n01/2
    nd1hm=(nd1-1)/2
    nd2hm=(nd2-1)/2
    nd3hm=(nd3-1)/2
    !Set to zero
    do i3=1,nd3
      do i2=1,nd2
        do i1=1,nd1
          zarray(1,i1,i2,i3) = 0.0_8
          zarray(2,i1,i2,i3) = 0.0_8
        end do
      end do
    end do
    !Set zarray
    do i3=1,n03
      do i2=1,n02
        do i1=1,n01h
          zarray(1,i1+nd1hm,i2+nd2hm,i3+nd3hm) = density(2*i1-1,i2,i3)
          zarray(2,i1+nd1hm,i2+nd2hm,i3+nd3hm) = density(2*i1,i2,i3)
        end do
      end do
    end do
    if(modulo(n01,2) == 1) then
      do i3=1,n03
        do i2=1,n02
          zarray(1,n01h+1+nd1hm,i2+nd2hm,i3+nd3hm) = density(n01,i2,i3)
        end do
      end do
    end if
  end subroutine zarray_in
  !!***


  !!****h* BigDFT/zarray_out
  !! NAME
  !!   zarray_out
  !!
  !! FUNCTION
  !!   Set the potential (rhopot) from zarray
  !!   Calculate the Hartree energy.
  !!
  !! SOURCE
  !!
  subroutine zarray_out(n01, n02, n03, nd1, nd2, nd3, rhopot, zarray, factor)
    integer, intent(in)    :: n01,n02,n03,nd1,nd2,nd3
    real(8), intent(out)   :: rhopot(n01,n02,n03)
    real(8), intent(in)    :: zarray(2*nd1,nd2,nd3)  ! Convert zarray(2,nd1,nd2,nd3) -> zarray(2*nd1,nd2,nd3) to use i1=1,n01
                                                     ! instead of i1=1,n1h + special case for modulo(n01,2)
    real(8), intent(in)    :: factor
    
    integer :: i1,i2,i3
    !
    do i3=1,n03
      do i2=1,n02
        do i1=1,n01
          rhopot(i1, i2, i3) = factor*zarray(i1,i2,i3)
        end do
      end do
    end do

  end subroutine zarray_out
  !!***

  !!****h* BigDFT/build_kernel
  !! NAME
  !!   build_kernel
  !!
  !! FUNCTION
  !!    Build the kernel of a gaussian function 
  !!    for interpolating scaling functions.
  !!
  !! SYNOPSIS
  !!    Build the kernel (karrayout) of a gaussian function 
  !!    for interpolating scaling functions
  !!    $$ K(j) = \int \int \phi(x) g(x`-x) \delta(x`- j) dx dx` $$
  !!
  !!    n01,n02,n03     Mesh dimensions of the density
  !!    n1k,n2k,n3k     Dimensions of the kernel
  !!    hgrid           Mesh step
  !!    itype_scf       Order of the scaling function (8,14,16)
  !!
  !! AUTHORS
  !!    T. Deutsch, L. Genovese
  !! COPYRIGHT
  !!    Copyright (C) 2005 CEA
  !! CREATION DATE
  !!    13/07/2005
  !!
  !! MODIFICATION HISTORY
  !!    13/09/2005 Use hgrid instead of acell
  !!    09/12/2005 Real kernel, stocked only half components
  !!    13/12/2005 Add routines to simplify the port into Stefan`s code
  !!
  !! SOURCE
  !!
  subroutine build_kernel(n01,n02,n03,nfft1,nfft2,nfft3,hgrid,itype_scf,karrayout)
    integer, intent(in) :: n01,n02,n03,nfft1,nfft2,nfft3,itype_scf
    real(8), intent(in) :: hgrid
    real(8), dimension(nfft1/2+1,nfft2/2+1,nfft3/2+1), intent(out) :: karrayout

    !Do not touch !!!!
    integer, parameter :: N_GAUSS = 89
    !Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
    integer, parameter :: n_points = 2**6

    !Better p_gauss for calculation
    !(the support of the exponential should be inside [-n_range/2,n_range/2])
    real(8), parameter :: p0_ref = 1.d0
    real(8), dimension(N_GAUSS) :: p_gauss,w_gauss

    real(8), allocatable :: kernel_scf(:), kern_1_scf(:)
    real(8), allocatable :: x_scf(:), y_scf(:)
    real(8), allocatable :: karrayhalf(:, :, :)

    real(8) :: ur_gauss,dr_gauss,acc_gauss,pgauss,kern,a_range
    real(8) :: factor,factor2,dx,absci,p0gauss,p0_cell
    real(8) :: a1,a2,a3
    integer :: nd1,nd2,nd3,n1k,n2k,n3k,n_scf
    integer :: i_gauss,n_range,n_cell
    integer :: i,n_iter,i1,i2,i3,i_kern
    integer :: i01,i02,i03,inkee,n1h,n2h,n3h,nd1h

    !Number of integration points : 2*itype_scf*n_points
    n_scf=2*itype_scf*n_points
    !Dimensions of Kernel
    n1k=nfft1/2+1 ; n2k=nfft2/2+1 ; n3k=nfft3/2+1
    n1h=nfft1/2  ; n2h=nfft2/2 ; n3h=nfft3/2
    nd1 = nfft1 + modulo(nfft1+1,2)
    nd2 = nfft2 + modulo(nfft2+1,2)
    nd3 = nfft3 + modulo(nfft3+1,2)

    !Half size for the half FFT
    nd1h=(nd1+1)/2

    !Allocations
    SAFE_ALLOCATE(x_scf(0:n_scf))
    SAFE_ALLOCATE(y_scf(0:n_scf))

    !Build the scaling function
    call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
    !Step grid for the integration
    dx = real(n_range,8)/real(n_scf,8)
    !Extend the range (no more calculations because fill in by 0.0_8)
    n_cell = max(n01,n02,n03)
    n_range = max(n_cell,n_range)

    !Allocations
    SAFE_ALLOCATE(kernel_scf(-n_range:n_range))
    SAFE_ALLOCATE(kern_1_scf(-n_range:n_range))

    !Lengthes of the box (use FFT dimension)
    a1 = hgrid * real(n01,8)
    a2 = hgrid * real(n02,8)
    a3 = hgrid * real(n03,8)

    x_scf(:) = hgrid * x_scf(:)
    y_scf(:) = 1.d0/hgrid * y_scf(:)
    dx = hgrid * dx
    !To have a correct integration
    p0_cell = p0_ref/(hgrid*hgrid)

    !Initialisation of the gaussian (Beylkin)
    call gequad(N_GAUSS,p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
    !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3) 
    !(biggest length in the cube)
    !We divide the p_gauss by a_range**2 and a_gauss by a_range
    a_range = sqrt(a1*a1+a2*a2+a3*a3)
    factor = 1.d0/a_range
    !factor2 = factor*factor
    factor2 = 1.d0/(a1*a1+a2*a2+a3*a3)
    do i_gauss=1,N_GAUSS
      p_gauss(i_gauss) = factor2*p_gauss(i_gauss)
    end do
    do i_gauss=1,N_GAUSS
      w_gauss(i_gauss) = factor*w_gauss(i_gauss)
    end do

    karrayout(:,:,:) = 0.0_8

    !Use in this order (better for accuracy).
    loop_gauss: do i_gauss=N_GAUSS,1,-1
      !Gaussian
      pgauss = p_gauss(i_gauss)

      !We calculate the number of iterations to go from pgauss to p0_ref
      n_iter = nint((log(pgauss) - log(p0_cell))/log(4.d0))
      if (n_iter <= 0)then
        n_iter = 0
        p0gauss = pgauss
      else
        p0gauss = pgauss/4.d0**n_iter
      end if

      !Stupid integration
      !Do the integration with the exponential centered in i_kern
      kernel_scf(:) = 0.0_8
      do i_kern=0,n_range
        kern = 0.0_8
        do i=0,n_scf
          absci = x_scf(i) - real(i_kern,8)*hgrid
          absci = absci*absci
          kern = kern + y_scf(i)*exp(-p0gauss*absci)*dx
        end do
        kernel_scf(i_kern) = kern
        kernel_scf(-i_kern) = kern
        if (abs(kern) < 1.d-18) then
          !Too small not useful to calculate
          exit
        end if
      end do

      !Start the iteration to go from p0gauss to pgauss
      call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,kern_1_scf)

      !Add to the kernel.
      do i3=1,n03
        i03 = i3-1
        do i2=1,n02
          i02 = i2-1
          do i1=1,n01
            i01 = i1-1
            karrayout(i1,i2,i3) = karrayout(i1,i2,i3) + w_gauss(i_gauss)* & 
              kernel_scf(i01)*kernel_scf(i02)*kernel_scf(i03)
          end do
        end do
      end do

    end do loop_gauss

    SAFE_DEALLOCATE_A(kernel_scf)
    SAFE_DEALLOCATE_A(kern_1_scf)
    SAFE_DEALLOCATE_A(x_scf)
    SAFE_DEALLOCATE_A(y_scf)

    !Set karray
    SAFE_ALLOCATE(karrayhalf(1:2, 1:nd1h*nd2*nd3, 1:2))

    !Set karray : use mirror symmetries
    inkee=1
    call karrayhalf_in(n01,n02,n03,n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,&
      karrayout,karrayhalf)
    call fft(n1h,nfft2,nfft3,nd1h,nd2,nd3,karrayhalf,1,inkee)
    !Reconstruct the real kernel
    call kernel_recon(n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,&
      karrayhalf(1,1,inkee),karrayout)

    SAFE_DEALLOCATE_A(karrayhalf)
  end subroutine Build_Kernel
  !!***

  
  !!****h* BigDFT/calculate_dimensions
  !! NAME
  !!   calculate_dimensions
  !!
  !! FUNCTION
  !!   Give the dimensions of the FFT
  !!
  !! SOURCE
  !!
  subroutine calculate_dimensions(n01,n02,n03,nfft1,nfft2,nfft3)
    integer, intent(in) :: n01,n02,n03
    integer, intent(out) :: nfft1,nfft2,nfft3

    integer :: i1,i2,i3,l1
    !Test 2*n01, 2*n02, 2*n03
    !write(*,*) 'in dimensions_fft',n01,n02,n03
    i1=2*n01
    i2=2*n02
    i3=2*n03
    do
      call fourier_dim(i1,nfft1)
      call fourier_dim(nfft1/2,l1)
      if (modulo(nfft1,2) == 0 .and. modulo(l1,2) == 0 .and. 2*l1 == nfft1) then
        exit
      end if
      i1=i1+1
    end do
    do
      call fourier_dim(i2,nfft2)
      if (modulo(nfft2,2) == 0) then
        exit
      end if
      i2=i2+1
    end do
    do
      call fourier_dim(i3,nfft3)
      if (modulo(nfft3,2) == 0) then
        exit
      end if
      i3=i3+1
    end do
    !nd1 = nfft1 + modulo(nfft1+1,2)
    !nd2 = nfft2 + modulo(nfft2+1,2)
    !nd3 = nfft3 + modulo(nfft3+1,2)
    !write(*,*) 'out dimensions_fft',nfft1,nfft2,nfft3
  end subroutine calculate_dimensions
  !!***


  !!****h* BigDFT/karrayhalf_in
  !! NAME
  !!   karrayhalf_in
  !!
  !! FUNCTION
  !!    Put in the array for4446666444 FFT
  !!
  !! SOURCE
  !!
  subroutine karrayhalf_in(n01,n02,n03,n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3, kernel,karrayhalf)
    integer, intent(in) :: n01,n02,n03,n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3
    real(8), dimension(n1k,n2k,n3k), intent(in) :: kernel
    real(8), dimension(2,(nd1+1)/2,nd2,nd3), intent(out) :: karrayhalf

    real(8), dimension(:), allocatable :: karray
    integer :: i1,i2,i3,nd1h,n1h,n2h,n3h
    !Body
    n1h=nfft1/2
    n2h=nfft2/2
    n3h=nfft3/2

    SAFE_ALLOCATE(karray(1:nfft1))

    nd1h=(nd1+1)/2
    karrayhalf(:,:,:,:) = 0.0_8
    do i3=1,n03
      do i2=1,n02
        karray(:) = 0.0_8
        do i1=1,n01
          karray(i1+n1h) = kernel(i1,i2,i3)
        end do
        do i1=2,n01
          karray(n1h-i1+1+nd1-nfft1) = kernel(i1,i2,i3)
        end do
        do i1=1,n1h
          karrayhalf(1,i1,i2+n2h,i3+n3h) = karray(2*i1-1)
          karrayhalf(2,i1,i2+n2h,i3+n3h) = karray(2*i1)
        end do
      end do
      do i2=2,n02
        do i1=1,nd1h
          karrayhalf(:,i1,n2h-i2+1+nd2-nfft2,i3+n3h) = &
            karrayhalf(:,i1,i2+n2h,i3+n3h)
        end do
      end do
    end do
    do i3=2,n03
      do i2=1,nd2
        do i1=1,nd1h
          karrayhalf(:,i1,i2,n3h-i3+1+nd3-nfft3) = karrayhalf(:,i1,i2,i3+n3h)
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(karray)
  end subroutine karrayhalf_in
  !!***

  
  !!****h* BigDFT/kernel_recon
  !! NAME
  !!   kernel_recon
  !!
  !! FUNCTION
  !!    Reconstruction of the kernel from the FFT array zarray
  !!    We keep only the half kernel in each direction (x,y,z).
  !!
  !! SOURCE
  !!
  subroutine kernel_recon(n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,zarray,karray)
    integer, intent(in) :: n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3
    real(8), dimension(2,(nd1+1)/2*nd2*nd3), intent(in) :: zarray
    real(8), dimension(n1k,n2k,n3k), intent(out) :: karray

    real(8), dimension(:), allocatable :: cos_array,sin_array
    integer :: i1,i2,i3,ind1,ind2,nd1h,n1h,n2h,n3h
    real(8) :: rfe,ife,rfo,ifo,cp,sp,rk,ik,a,b,c,d,pi2
    !Body
    n1h=nfft1/2
    n2h=nfft2/2
    n3h=nfft3/2
    nd1h=(nd1+1)/2
    pi2=8.d0*datan(1.d0)
    pi2=pi2/real(nfft1,8)

    SAFE_ALLOCATE(cos_array(1:nd1h))
    SAFE_ALLOCATE(sin_array(1:nd1h))

    do i1=1,nd1h
      cos_array(i1)= dcos(pi2*(i1-1))
      sin_array(i1)=-dsin(pi2*(i1-1))
    end do
    do i3=1,n3h+1
      do i2=1,n2h+1
        do i1=1,nd1h
          call norm_ind(nd1h,nd2,nd3,i1,i2,i3,ind1)
          call symm_ind(nd1h,nd2,nd3,i1,i2,i3,ind2)
          a=zarray(1,ind1)
          b=zarray(2,ind1)
          c=zarray(1,ind2)
          d=zarray(2,ind2)
          rfe=0.5d0*(a+c)
          ife=0.5d0*(b-d)
          rfo=0.5d0*(a-c)
          ifo=0.5d0*(b+d) 
          cp=cos_array(i1)
          sp=sin_array(i1)
          rk=rfe+cp*ifo-sp*rfo
          ik=ife-cp*rfo-sp*ifo
          !For big dimension 1.d-9 otherwise 1.d-10
          !Remove the test
          !if(abs(ik) >= 1.d-10) then
          !   print *,"non real kernel FFT",i1,i2,i3,ik  
          !   stop
          !end if
          !Build the intermediate FFT convolution (full)
          !call norm_ind(nd1,nd2,nd3,i1,i2,i3,indA)
          karray(i1,i2,i3)=rk 
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(cos_array)
    SAFE_DEALLOCATE_A(sin_array)
  end subroutine kernel_recon
  !!***

  !!****h* BigDFT/par_calculate_dimensions
  !! NAME
  !!   par_calculate_dimensions
  !!
  !! FUNCTION
  !!    Calculate four sets of dimension needed for the calculation of the
  !!    zero-padded convolution
  !!
  !! SYNOPSIS
  !!    n01,n02,n03 original real dimensions (input)
  !!
  !!    m1,m2,m3 original real dimension with the dimension 2 and 3 exchanged
  !!
  !!    n1,n2 the first FFT even dimensions greater that 2*m1, 2*m2
  !!    n3    the double of the first FFT even dimension greater than m3
  !!          (improved for the HalFFT procedure)
  !!
  !!    md1,md2,md3 half of n1,n2,n3 dimension. They contain the real unpadded space,
  !!                properly enlarged to be compatible with the FFT dimensions n_i.
  !!                md2 is further enlarged to be a multiple of nproc
  !!
  !!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
  !!                formally 1/8 of the fourier grid. Here the dimension nd3 is
  !!                enlarged to be a multiple of nproc
  !!                
  !! WARNING
  !!    The dimension m2 and m3 correspond to n03 and n02 respectively
  !!    this is needed since the convolution routine manage arrays of dimension
  !!    (md1,md3,md2/nproc)
  !!
  !! AUTHOR
  !!    Luigi Genovese
  !! CREATION DATE
  !!    February 2006
  !!
  !! SOURCE
  !!
  subroutine par_calculate_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3, md1,md2,md3,nd1,nd2,nd3,nproc)
    integer, intent(in) :: n01,n02,n03,nproc
    integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3

    integer :: l1,l2,l3

    !dimensions of the density in the real space, inverted for convenience

    m1=n01
    m2=n03
    m3=n02

    ! real space grid dimension (suitable for number of processors)

    !     n1=2*m1; n2=2*m2; n3=2*m3

    l1=2*m1
    l2=2*m2
    l3=m3 !beware of the half dimension
    do
      !this is for the FFT of the kernel
      !one can erase it when the kernel is parallelized
      call fourier_dim(l1,n1)
      !call fourier_dim(n1/2,l1A)
      if (modulo(n1,2) == 0&! .and. 2*l1A == n1
        ) then
        exit
      end if
      l1=l1+1
    end do
    do
      call fourier_dim(l2,n2)
      if (modulo(n2,2) == 0) then
        exit
      end if
      l2=l2+1
    end do
    do
      call fourier_dim(l3,n3)
      !call fourier_dim(n3/2,l3A)
      if (modulo(n3,2) == 0 &!.and. 2*l3A == n3 .and. modulo(l3A,2) == 0
        ) then
        exit
      end if
      l3=l3+1
    end do
    n3=2*n3

    !dimensions that contain the unpadded real space,
    ! compatible with the number of processes
    md1=n1/2
    md2=n2/2
    md3=n3/2
151 if (nproc*(md2/nproc).lt.n2/2) then
      md2=md2+1
      goto 151
    endif


    !dimensions of the kernel, 1/8 of the total volume,
    !compatible with nproc

    nd1=n1/2+1;  nd2=n2/2+1
    nd3=n3/2+1
250 if (modulo(nd3,nproc) .ne. 0) then
      nd3=nd3+1
      goto 250
    endif

  end subroutine par_calculate_dimensions
  !!***


  !!****h* BigDFT/par_psolver_kernel
  !! NAME
  !!   par_psolver_kernel
  !!
  !! FUNCTION
  !!    Solver of Poisson equation applying a kernel, parallel computation
  !!
  !! SYNOPSIS
  !!    Poisson solver applying a kernel and
  !!    using Fourier transform for the convolution.
  !!    rhopot : input  -> the density
  !!             output -> the Hartree potential + pot_ion
  !!    All the processes manage the same global rhopot array
  !!    The potential pot_ion is ADDED in the array rhopot.
  !!    Calculate also the Hartree potential
  !!
  !!    Replaces the charge density contained in rhopot
  !!    by the Hartree stored as well in rhopot.
  !!    If xc_on is true, it also adds the XC potential and
  !!    ionic potential pot_ion
  !!
  !!    kernelLOC: the kernel in fourier space, calculated from ParBuil_Kernel routine
  !!               it is a local vector (each process have its own part)
  !!
  !!    comm: MPI communicator to use
  !!
  !!    We double the size of the mesh except in one dimension
  !!    in order to use the property of the density to be real.
  !! WARNING
  !!
  !! AUTHOR
  !!    Luigi Genovese
  !! CREATION DATE
  !!    February 2006
  !!
  !! SOURCE
  !!
  subroutine par_psolver_kernel(n01, n02, n03, nd1, nd2, nd3, hgrid, kernelLOC, rhopot, iproc, nproc, comm)
    integer, intent(in)  :: n01,n02,n03,iproc,nproc
    integer, intent(inout) :: nd1,nd2,nd3
    real(8), intent(in) :: hgrid
    real(8), intent(in), dimension(nd1,nd2,nd3/nproc) :: kernelLOC
    real(8), intent(inout), dimension(n01,n02,n03) :: rhopot
    integer, intent(in) :: comm
    
    integer :: m1,m2,m3,n1,n2,n3,md1,md2,md3
    
    call par_calculate_dimensions(n01, n02, n03, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
    call pconvxc_off(m1, m2, m3, n1, n2, n3, nd1, nd2, nd3, md1, md2, md3, iproc, nproc, rhopot, kernelLOC, hgrid, comm)
    
  end subroutine par_psolver_kernel
  !!***


  !!****h* BigDFT/pconvxc_off
  !! NAME
  !!   pconvxc_off
  !!
  !! FUNCTION
  !!    Calculate the parallel convolution with the kernel
  !!    without the exchange-correlation part
  !!
  !! SYNOPSIS
  !!    Poisson solver applying a kernel and
  !!    using Fourier transform for the convolution.
  !!    rhopot : input  -> the density
  !!             output -> the Hartree potential + pot_ion
  !!    All the processes manage the same global rhopot array
  !!    The potential pot_ion is ADDED in the array rhopot.
  !!    Calculate also the Hartree potential
  !!
  !!    Replaces the charge density contained in rhopot
  !!    by the Hartree stored as well in rhopot.
  !!
  !!    kernelLOC: the kernel in fourier space, calculated from ParBuild_Kernel routine
  !!               it is a local vector (each process have its own part)
  !!
  !!    comm: MPI communicator to use
  !!
  !!    We double the size of the mesh except in one dimension
  !!    in order to use the property of the density to be real.
  !! WARNING
  !!
  !! AUTHOR
  !!    Luigi Genovese
  !! CREATION DATE
  !!    February 2006
  !!
  !! SOURCE
  !!
  subroutine pconvxc_off(m1, m2, m3, n1, n2, n3, nd1, nd2, nd3, md1, md2, md3, iproc, nproc, rhopot, kernelloc, hgrid, comm)
    integer, intent(in) :: m1,m2,m3,n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,iproc,nproc
    real(8), dimension(nd1,nd2,nd3/nproc), intent(in) :: kernelloc
    real(8), dimension(m1,m3,m2), intent(inout) :: rhopot
    real(8), intent(in) :: hgrid
    integer, intent(in) :: comm

#if defined(HAVE_MPI)
    !Local variables
    integer :: istart,iend,jend,jproc
    real(8) :: scal
    real(8), dimension(:,:,:), allocatable :: zf, lrhopot(:, :, :)
    integer, dimension(:,:), allocatable :: gather_arr
    type(profile_t), save :: prof

    !factor to be used to keep unitarity
    scal=hgrid**3/real(n1*n2*n3,8)

    SAFE_ALLOCATE(zf(1:md1, 1:md3, 1:md2/nproc))
    SAFE_ALLOCATE(gather_arr(0:nproc-1, 1:2))

    !Here we insert the process-related values of the density, starting from the total density
    call enterdensity(rhopot(1,1,1), m1, m2, m3, md1, md2, md3, iproc, nproc, zf(1,1,1))

    !this routine builds the values for each process of the potential (zf), multiplying by the factor
    call convolxc_off(n1, n2, n3, nd1, nd2, nd3, md1, md2, md3, nproc, iproc, kernelloc, zf, scal, hgrid, comm)

    !building the array of the data to be sent from each process
    !and the array of the displacement
    do jproc = 0, nproc - 1
      istart = min(jproc*(md2/nproc),m2-1)
      jend = max(min(md2/nproc, m2 - md2/nproc*jproc), 0)
      gather_arr(jproc, 1) = m1*m3*jend
      gather_arr(jproc, 2) = istart*m1*m3
    end do

    !assign the distributed density to the rhopot array
    istart=min(iproc*(md2/nproc),m2-1)
    jend=max(min(md2/nproc,m2-md2/nproc*iproc),0)
    iend=istart+jend

    if(jend == 0) jend = 1

    SAFE_ALLOCATE(lrhopot(1:m1, 1:m3, 1:jend))

    lrhopot(1:m1, 1:m3, 1:jend) = zf(1:m1, 1:m3, 1:jend)

    call profiling_in(prof, "ISF_GATHER")
    call MPI_Allgatherv(lrhopot(1, 1, 1), gather_arr(iproc, 1), MPI_DOUBLE_PRECISION, rhopot(1, 1, 1), gather_arr(0, 1),&
      gather_arr(0, 2), MPI_DOUBLE_PRECISION, comm, mpi_err)
    call profiling_out(prof)

    SAFE_DEALLOCATE_A(zf)
    SAFE_DEALLOCATE_A(lrhopot)
    SAFE_DEALLOCATE_A(gather_arr)

#endif
  end subroutine pconvxc_off
!!***


  !!****h* BigDFT/enterdensity
  !! NAME
  !!   enterdensity
  !!
  !! FUNCTION
  !!
  !!   Define a real space process-dependent vector zf with the global dimensions that are half of the FFT grid
  !!   in order to perform convolution. The dimension md2 is a multiple of nproc
  !!   Can be used also to define the local part of pot_ion
  !!
  !! AUTHOR
  !!    L. Genovese
  !! CREATION DATE
  !!    February 2006
  !!
  !! SOURCE
  !!
  subroutine enterdensity(rhopot,m1,m2,m3,md1,md2,md3,iproc,nproc,zf)
    integer, intent(in) :: m1,m2,m3,md1,md2,md3,iproc,nproc
    real(8), dimension(0:md1-1,0:md3-1,0:md2/nproc-1), intent(out) :: zf
    real(8), dimension(0:m1-1,0:m3-1,0:m2-1), intent(in) :: rhopot

    integer :: j1,j2,j3,jp2

    !Body
    do jp2=0,md2/nproc-1
      j2=iproc*(md2/nproc)+jp2
      if (j2.le.m2-1) then
        do j3=0,m3-1
          do j1=0,m1-1
            zf(j1,j3,jp2)=rhopot(j1,j3,j2)
          end do
          do j1=m1,md1-1
            zf(j1,j3,jp2)=0.0_8
          end do
        end do
        do j3=m3,md3-1
          do j1=0,md1-1
            zf(j1,j3,jp2)=0.0_8
          end do
        end do
      else
        do j3=0,md3-1
          do j1=0,md1-1
            zf(j1,j3,jp2)=0.0_8
          end do
        end do
      end if
    end do

  end subroutine enterdensity
  !!***

  
  !!****h* BigDFT/par_build_kernel
  !! NAME
  !!   par_build_kernel
  !!
  !! FUNCTION
  !!    Build the kernel of a gaussian function
  !!    for interpolating scaling functions.
  !!    Do the parallel HalFFT of the symmetrized function and stores into
  !!    memory only 1/8 of the grid divided by the number of processes nproc
  !!
  !! SYNOPSIS
  !!    Build the kernel (karray) of a gaussian function
  !!    for interpolating scaling functions
  !!    $$ K(j) = \sum_k \omega_k \int \int \phi(x) g_k(y-x) \delta(y-j) dx dy $$
  !!
  !!    n01,n02,n03        Mesh dimensions of the density
  !!    nfft1,nfft2,nfft3  Dimensions of the FFT grid (HalFFT in the third direction)
  !!    n1k,n2k,n3k        Dimensions of the kernel FFT
  !!    hgrid              Mesh step
  !!    itype_scf          Order of the scaling function (8,14,16)
  !!    comm               MPI communicator to use
  !!
  !! AUTHORS
  !!    T. Deutsch, L. Genovese
  !! CREATION DATE
  !!    February 2006
  !!
  !! SOURCE
  !!
  subroutine par_build_kernel(n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,hgrid,itype_scf, &
    iproc,nproc,comm,karrayoutLOC)
    integer, intent(in)    :: n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,itype_scf,iproc,nproc
    real(8), intent(in)    :: hgrid
    real(8), intent(out)   :: karrayoutLOC(1:n1k, 1:n2k, 1:n3k/nproc)
    integer, intent(in)    :: comm
 
    !Do not touch !!!!
    integer, parameter :: N_GAUSS = 89
    !Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
    integer, parameter :: N_POINTS = 2**6

    !Better p_gauss for calculation
    !(the support of the exponential should be inside [-n_range/2,n_range/2])
    real(8), parameter :: p0_ref = 1.d0
    real(8) :: p_gauss(N_GAUSS), w_gauss(N_GAUSS)

    real(8), dimension(:), allocatable :: kernel_scf,kern_1_scf
    real(8), dimension(:), allocatable :: x_scf ,y_scf
    real(8), dimension(:,:,:,:), allocatable :: karrayfour
    real(8), dimension(:,:,:), allocatable :: karray

    real(8) :: ur_gauss,dr_gauss,acc_gauss,pgauss,kern,a_range
    real(8) :: factor,factor2,dx,absci,p0gauss,p0_cell
    real(8) :: a1,a2,a3
    integer :: n_scf,nker1,nker2,nker3
    integer :: i_gauss,n_range,n_cell,istart,iend,istart1,istart2,iend1,iend2
    integer :: i,n_iter,i1,i2,i3,i_kern
    integer :: i01,i02,i03,n1h,n2h,n3h

    !Number of integration points : 2*itype_scf*n_points
    n_scf=2*itype_scf*n_points
    !Set karray
    !dimension test

    !here we must set the dimensions for the fft part, starting from the nfft
    !remember that actually nfft2 is associated to n03 and viceversa

    !dimensions that define the center of symmetry
    n1h=nfft1/2
    n2h=nfft2/2
    n3h=nfft3/2

    !Auxiliary dimensions only for building the FFT part
    nker1=nfft1
    nker2=nfft2
    nker3=nfft3/2+1

    !adjusting the last two dimensions to be multiples of nproc
    do
      if(modulo(nker2,nproc) == 0) exit
      nker2=nker2+1
    end do
    do
      if(modulo(nker3,nproc) == 0) exit
      nker3=nker3+1
    end do

    !this will be the array of the kernel in the real space
    allocate(karray(nker1,nfft3,nker2/nproc))

    !defining proper extremes for the calculation of the
    !local part of the kernel

    istart=iproc*nker2/nproc+1
    iend=min((iproc+1)*nker2/nproc,n2h+n03)

    istart1=istart
    if(iproc .eq. 0) istart1=n2h-n03+2

    iend2=iend

    iend1=n2h
    istart2=n2h+1
    if(istart .gt. n2h) then
      iend1=istart1-1
      istart2=istart
    end if
    if(iend .le. n2h) then
      istart2=iend2+1
      iend1=iend
    end if

!!!!!START KERNEL CONSTRUCTION
    !  if(iproc .eq. 0) then
    !     write(unit=*,fmt="(1x,a,i0,a)") &
    !          "Build the kernel in parallel using a sum of ",N_GAUSS," gaussians"
    !     write(unit=*,fmt="(1x,a,i0,a)") &
    !          "Use interpolating scaling functions of ",itype_scf," order"
    !  end if

    SAFE_ALLOCATE(x_scf(0:n_scf))
    SAFE_ALLOCATE(y_scf(0:n_scf))

    !Build the scaling function
    call scaling_function(itype_scf, n_scf, n_range, x_scf, y_scf)
    !Step grid for the integration
    dx = real(n_range,8)/real(n_scf,8)
    !Extend the range (no more calculations because fill in by 0.0_8)
    n_cell = max(n01,n02,n03)
    n_range = max(n_cell,n_range)

    !Allocations
    SAFE_ALLOCATE(kernel_scf(-n_range:n_range))
    SAFE_ALLOCATE(kern_1_scf(-n_range:n_range))

    !Lengthes of the box (use FFT dimension)
    a1 = hgrid * real(n01,8)
    a2 = hgrid * real(n02,8)
    a3 = hgrid * real(n03,8)

    x_scf(:) = hgrid * x_scf(:)
    y_scf(:) = 1.d0/hgrid * y_scf(:)
    dx = hgrid * dx
    !To have a correct integration
    p0_cell = p0_ref/(hgrid*hgrid)

    !Initialization of the gaussian (Beylkin)
    call gequad(N_GAUSS,p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
    !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3)
    !(biggest length in the cube)
    !We divide the p_gauss by a_range**2 and a_gauss by a_range
    a_range = sqrt(a1*a1+a2*a2+a3*a3)
    factor = 1.d0/a_range
    !factor2 = factor*factor
    factor2 = 1.d0/(a1*a1+a2*a2+a3*a3)
    do i_gauss=1,N_GAUSS
      p_gauss(i_gauss) = factor2*p_gauss(i_gauss)
    end do
    do i_gauss=1,N_GAUSS
      w_gauss(i_gauss) = factor*w_gauss(i_gauss)
    end do

    karray(:,:,:) = 0.0_8
    !Use in this order (better for accuracy).
    loop_gauss: do i_gauss = N_GAUSS, 1, -1
      !Gaussian
      pgauss = p_gauss(i_gauss)

      !We calculate the number of iterations to go from pgauss to p0_ref
      n_iter = nint((log(pgauss) - log(p0_cell))/log(4.d0))
      if (n_iter <= 0)then
        n_iter = 0
        p0gauss = pgauss
      else
        p0gauss = pgauss/4.d0**n_iter
      end if

      !Stupid integration
      !Do the integration with the exponential centered in i_kern
      kernel_scf(:) = 0.0_8
      do i_kern=0,n_range
        kern = 0.0_8
        do i=0,n_scf
          absci = x_scf(i) - real(i_kern,8)*hgrid
          absci = absci*absci
          kern = kern + y_scf(i)*exp(-p0gauss*absci)*dx
        end do
        kernel_scf(i_kern) = kern
        kernel_scf(-i_kern) = kern
        if (abs(kern) < 1.d-18) then
          !Too small not useful to calculate
          exit
        end if
      end do

      !Start the iteration to go from p0gauss to pgauss
      call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,kern_1_scf)

      !Add to the kernel (only the local part)

      do i3=istart1,iend1
        i03 =  n2h - i3 + 1
        do i2=1,n02
          i02 = i2-1
          do i1=1,n01
            i01 = i1-1
            karray(i1+n1h,i2+n3h,i3-istart+1) = karray(i1+n1h,i2+n3h,i3-istart+1) + w_gauss(i_gauss)* &
              kernel_scf(i01)*kernel_scf(i02)*kernel_scf(i03)
          end do
        end do
      end do
      do i3=istart2,iend2
        i03 = i3 - n2h -1
        do i2=1,n02
          i02 = i2-1
          do i1=1,n01
            i01 = i1-1
            karray(i1+n1h,i2+n3h,i3-istart+1) = karray(i1+n1h,i2+n3h,i3-istart+1) + w_gauss(i_gauss)* &
              kernel_scf(i01)*kernel_scf(i02)*kernel_scf(i03)
          end do
        end do
      end do


    end do loop_gauss

    !Build the kernel in the real space as an even function, thus having a real FFT

    do i3=istart1,iend2
      do i2=1,n02
        do i1=2,n01
          karray(n1h+2-i1,i2+n3h,i3-istart+1) = karray(i1+n1h,i2+n3h,i3-istart+1)
        end do
      end do
      do i2=2,n02
        do i1=1,nker1
          karray(i1,n3h+2-i2,i3-istart+1) = karray(i1,i2+n3h,i3-istart+1)
        end do
      end do
    end do


    !De-allocations
    SAFE_DEALLOCATE_A(kernel_scf)
    SAFE_DEALLOCATE_A(kern_1_scf)
    SAFE_DEALLOCATE_A(x_scf)
    SAFE_DEALLOCATE_A(y_scf)

!!!!END KERNEL CONSTRUCTION

    SAFE_ALLOCATE(karrayfour(1:2, 1:nker1, 1:nker2, 1:nker3/nproc))

    ! if(iproc .eq. 0) print *,"Do a 3D PHalFFT for the kernel"

    call kernelfft(nfft1,nfft2,nfft3,nker1,nker2,nker3,nproc,iproc,karray,karrayfour,comm)

    !Reconstruct the real kernel FFT
    do i3=1,n3k/nproc
      do i2=1,n2k
        do i1=1,n1k
          karrayoutLOC(i1,i2,i3)=karrayfour(1,i1,i2,i3)
        end do
      end do
    end do

    !De-allocations
    SAFE_DEALLOCATE_A(karray)
    SAFE_DEALLOCATE_A(karrayfour)
    
  end subroutine par_build_kernel
  !!***

  ! -------------------------------------------------------------------------

  subroutine gequad(n_gauss, p_gauss, w_gauss, ur_gauss, dr_gauss, acc_gauss)
    integer, intent(in)    :: n_gauss
    real(8), intent(out)   :: p_gauss(:)
    real(8), intent(out)   :: w_gauss(:)
    real(8), intent(out)   :: ur_gauss
    real(8), intent(out)   :: dr_gauss
    real(8), intent(out)   :: acc_gauss

    integer :: iunit, i, idx

    ur_gauss = 1.0_8
    dr_gauss = 1.0e-08_8
    acc_gauss = 1.0e-08_8
    
    iunit = io_open(trim(conf%share)//'/gequad.data', action = 'read', status = 'old', die = .true.)

    do i = 1, n_gauss
      read(iunit, *) idx, p_gauss(i), w_gauss(i)
    end do
    
    call io_close(iunit)

  end subroutine gequad

end module poisson_isf_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
