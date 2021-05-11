#include "global.h"

! Interface to libvdwxc.
! Naming:
!  * Functions that start with libvdwxc_ are public, to be called from other parts of Octopus.
!  * Interfaces that start with vdwxc_ are actual functions of libvdwxc.


module libvdwxc_oct_m
  use box_parallelepiped_oct_m
  use cube_oct_m
  use cube_function_oct_m
  use derivatives_oct_m
  use fft_oct_m
  use global_oct_m
  use grid_oct_m
  use io_function_oct_m
  use mesh_oct_m
  use mesh_cube_parallel_map_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pfft_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use unit_system_oct_m

  implicit none

  integer, parameter ::       &
    LIBVDWXC_MODE_AUTO = 1,   &
    LIBVDWXC_MODE_SERIAL = 2, &
    LIBVDWXC_MODE_MPI = 3!, &
    !LIBVDWXC_MODE_PFFT = 3

  private

  public ::                &
    libvdwxc_t,            &
    libvdwxc_init,         &
    libvdwxc_print,        &
    libvdwxc_write_info,   &
    libvdwxc_set_geometry, &
    libvdwxc_calculate,    &
    libvdwxc_end

  type libvdwxc_t
    private
    integer, pointer               :: libvdwxc_ptr
    type(mesh_t)                   :: mesh
    type(cube_t)                   :: cube
    type(mesh_cube_parallel_map_t) :: mesh_cube_map
    integer                        :: functional
    logical                        :: debug
    FLOAT, public                  :: energy
    FLOAT                          :: vdw_factor
  end type libvdwxc_t

#ifdef HAVE_LIBVDWXC
  include "vdwxcfort.f90"
#endif

contains

  subroutine libvdwxc_init(libvdwxc, namespace, functional)
    type(libvdwxc_t),          intent(out) :: libvdwxc
    type(namespace_t),         intent(in)  :: namespace
    integer,                   intent(in)  :: functional

    PUSH_SUB(libvdwxc_init)
#ifdef HAVE_LIBVDWXC
    call vdwxc_new(functional, libvdwxc%libvdwxc_ptr)
#else
    message(1) = "Octopus not compiled with libvdwxc"
    call messages_fatal(1, namespace=namespace)
#endif
    ASSERT(associated(libvdwxc%libvdwxc_ptr))
    libvdwxc%functional = functional

    !%Variable libvdwxcDebug
    !%Type logical
    !%Section Hamiltonian::XC
    !%Description
    !% Dump libvdwxc inputs and outputs to files.
    !%End
    call parse_variable(namespace, 'libvdwxcDebug', .false., libvdwxc%debug)
    POP_SUB(libvdwxc_init)

    !%Variable libvdwxcVDWFactor
    !%Type float
    !%Section Hamiltonian::XC
    !%Description
    !% Prefactor of non-local van der Waals functional.
    !% Setting a prefactor other than one is wrong, but useful
    !% for debugging.
    !%End
    call parse_variable(namespace, 'libvdwxcVDWFactor', M_ONE, libvdwxc%vdw_factor)
  end subroutine libvdwxc_init

  subroutine libvdwxc_print(this)
    type(libvdwxc_t), intent(in) :: this
    PUSH_SUB(libvdwxc_print)

#ifdef HAVE_LIBVDWXC
    call vdwxc_print(this%libvdwxc_ptr)
#endif

    POP_SUB(libvdwxc_print)
  end subroutine libvdwxc_print

  subroutine libvdwxc_write_info(this, iunit)
    type(libvdwxc_t), intent(in) :: this
    integer,          intent(in) :: iunit

    PUSH_SUB(libvdwxc_write_info)

    write(message(1), '(2x,a)') 'Correlation'
    if(this%functional == 1) then
      write(message(2), '(4x,a)') 'vdW-DF from libvdwxc'
    else if(this%functional == 2) then
      write(message(2), '(4x,a)') 'vdW-DF2 from libvdwxc'
    else if(this%functional == 3) then
      write(message(2), '(4x,a)') 'vdW-DF-cx from libvdwxc'
    else
      write(message(2), '(4x,a)') 'unknown libvdwxc functional'
    end if
    call messages_info(2, iunit)

    POP_SUB(libvdwxc_write_info)
  end subroutine libvdwxc_write_info

  subroutine libvdwxc_set_geometry(this, namespace, space, mesh)
    type(libvdwxc_t),  intent(inout) :: this
    type(namespace_t), intent(in)    :: namespace
    type(space_t),     intent(in)    :: space
    type(mesh_t),      intent(inout) :: mesh

    integer :: blocksize
    integer :: libvdwxc_mode

    PUSH_SUB(libvdwxc_set_geometry)
    this%mesh = mesh

    ! libvdwxc can use either FFTW-MPI or PFFT and Octopus should not
    ! care too much, as long as this particular cube is properly
    ! configured.  Unfortunately most of the args to cube_init are
    ! very entangled with FFT libraries.  All we want and need is to
    ! specify our own decomposition and pass that to libvdwxc, but
    ! we can only say "we want to use PFFT" and then it does the
    ! decomposition.
    libvdwxc_mode = LIBVDWXC_MODE_AUTO

    !%Variable libvdwxcMode
    !%Type integer
    !%Section Hamiltonian::XC
    !%Description
    !% Whether libvdwxc should run with serial fftw3, fftw3-mpi, or pfft.
    !% to specify fftw3-mpi in serial for debugging.
    !% pfft is not implemented at the moment.
    !%Option libvdwxc_mode_auto 1
    !% Use serial fftw3 if actually running in serial, else fftw3-mpi.
    !%Option libvdwxc_mode_serial 2
    !% Run with serial fftw3.  Works only when not parallelizing over domains.
    !%Option libvdwxc_mode_mpi 3
    !% Run with fftw3-mpi.  Works only if Octopus is compiled with MPI.
    !%End
    call parse_variable(namespace, 'libvdwxcMode', LIBVDWXC_MODE_AUTO, libvdwxc_mode)

    if(libvdwxc_mode == LIBVDWXC_MODE_AUTO) then
      if(mesh%mpi_grp%size == 1) then
        libvdwxc_mode = LIBVDWXC_MODE_SERIAL
      else
        libvdwxc_mode = LIBVDWXC_MODE_MPI
      end if
    end if

    ! TODO implement.  Should dump quantities to files.

    blocksize = mesh%idx%ll(3) / mesh%mpi_grp%size
    if(mod(mesh%idx%ll(3), mesh%mpi_grp%size) /= 0) then
      blocksize = blocksize + 1
    end if

    if(libvdwxc_mode == LIBVDWXC_MODE_SERIAL) then
      call cube_init(this%cube, mesh%idx%ll, namespace, space)
    else
#ifdef HAVE_MPI
      call cube_init(this%cube, mesh%idx%ll, namespace, space, mpi_grp = mesh%mpi_grp, &
        need_partition = .true., blocksize = blocksize)
      call mesh_cube_parallel_map_init(this%mesh_cube_map, mesh, this%cube)
#endif
    end if

    ! There is some low-level implementation issue where with PFFT,
    ! the leading dimension is one smaller than normal for some reason.
    ! Therefore we cannot use the PFFT stuff without a frightful mess.

#ifdef HAVE_LIBVDWXC
    select type (box => mesh%sb%box)
    type is (box_parallelepiped_t)
      call vdwxc_set_unit_cell(this%libvdwxc_ptr, &
        this%cube%rs_n_global(3), this%cube%rs_n_global(2), this%cube%rs_n_global(1), &
        mesh%sb%latt%rlattice(3, 3), mesh%sb%latt%rlattice(2, 3), mesh%sb%latt%rlattice(1, 3), &
        mesh%sb%latt%rlattice(3, 2), mesh%sb%latt%rlattice(2, 2), mesh%sb%latt%rlattice(1, 2), &
        mesh%sb%latt%rlattice(3, 1), mesh%sb%latt%rlattice(2, 1), mesh%sb%latt%rlattice(1, 1))
    class default
      call vdwxc_set_unit_cell(this%libvdwxc_ptr, &
        this%cube%rs_n_global(3), this%cube%rs_n_global(2), this%cube%rs_n_global(1), &
        mesh%spacing(3) * this%cube%rs_n_global(3), 0.0_8, 0.0_8, &
        0.0_8, mesh%spacing(2) * this%cube%rs_n_global(2), 0.0_8, &
        0.0_8, 0.0_8, mesh%spacing(1) * this%cube%rs_n_global(1))
    end select

    if(libvdwxc_mode == LIBVDWXC_MODE_SERIAL) then
      call vdwxc_init_serial(this%libvdwxc_ptr)
    else
#ifdef HAVE_LIBVDWXC_MPI
      call vdwxc_init_mpi(this%libvdwxc_ptr, mesh%mpi_grp%comm)
#else
      message(1) = "libvdwxc was not compiled with MPI"
      message(2) = "Recompile libvdwxc with MPI for vdW with domain decomposition"
      call messages_fatal(2, namespace=namespace)
#endif
    end if
    call vdwxc_print(this%libvdwxc_ptr)
#endif

    POP_SUB(libvdwxc_set_geometry)
  end subroutine libvdwxc_set_geometry

  subroutine libvdwxc_calculate(this, namespace, space, rho, gradrho, dedd, dedgd)
    type(libvdwxc_t),         intent(inout) :: this
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    FLOAT, dimension(:,:),    intent(inout) :: rho !!! data type
    FLOAT, dimension(:,:,:),  intent(in)    :: gradrho
    FLOAT, dimension(:,:),    intent(inout) :: dedd
    FLOAT, dimension(:,:,:),  intent(inout) :: dedgd

    type(cube_function_t) :: cf!rhocf, sigmacf, dedrhocf, dedsigmacf

    real(8), allocatable :: workbuffer(:)
    real(8), allocatable :: cube_rho(:,:,:), cube_sigma(:,:,:), cube_dedrho(:,:,:), cube_dedsigma(:,:,:)
    real(8), dimension(3) :: energy_and_integrals_buffer
    integer :: ii
#ifdef HAVE_MPI
    integer :: ierr
#endif

    PUSH_SUB(libvdwxc_calculate)

    ASSERT(size(rho, 2) == 1)
    ASSERT(size(gradrho, 3) == 1)
    ASSERT(size(dedd, 2) == 1)
    ASSERT(size(dedgd, 3) == 1)

    ! Well.  I thought we would be using four different cube functions
    ! for rho, sigma, dedrho and dedsigma.
    !
    ! But since the cube has PFFT associated, any attempt to use it
    ! (update: We do not actually use PFFT anymore, so maybe this will
    !  not be that broken now)
    ! results in some kind of behind-the-scenes use of a possibly
    ! global FFT-related buffer.
    ! (Note: actually we now use cubes without FFT lib)
    !
    ! cube functions take a force_alloc variable which appears to make
    ! them allocate their own buffer (like we want), but then Octopus
    ! segfaults on cube function free, which I suspect is a bug in the
    ! code, or at the very least due to something so undocumented that
    ! I cannot reasonably figure it out.
    !
    ! We could just disable PFFT for the cube (since we will never
    ! actually call any FFT from Octopus) and do our own
    ! redistribution, but the redistribution code is loaded with
    ! references to FFT library, so this would probably be a bit
    ! optimistic.
    !
    ! So we create here our own arrays over which we have reasonable control.

    SAFE_ALLOCATE(workbuffer(1:this%mesh%np))
    SAFE_ALLOCATE(cube_rho(1:this%cube%rs_n(1), 1:this%cube%rs_n(2), 1:this%cube%rs_n(3)))
    SAFE_ALLOCATE(cube_sigma(1:this%cube%rs_n(1), 1:this%cube%rs_n(2), 1:this%cube%rs_n(3)))
    SAFE_ALLOCATE(cube_dedrho(1:this%cube%rs_n(1), 1:this%cube%rs_n(2), 1:this%cube%rs_n(3)))
    SAFE_ALLOCATE(cube_dedsigma(1:this%cube%rs_n(1), 1:this%cube%rs_n(2), 1:this%cube%rs_n(3)))
    ! This is sigma, the absolute-squared density gradient:
    workbuffer(:) = sum(gradrho(:, :, 1)**2, 2)

    if(this%debug) then
      call libvdwxc_write_array(rho(:, 1), 'rho')
      call libvdwxc_write_array(workbuffer, 'gradrho')
    end if

    cube_rho = M_ZERO
    cube_sigma = M_ZERO
    cube_dedrho = M_ZERO
    cube_dedsigma = M_ZERO

    call dcube_function_alloc_RS(this%cube, cf, in_device = .false.)

    call tocube(rho(:, 1), cube_rho)
    call tocube(workbuffer, cube_sigma)

    this%energy = M_ZERO
#ifdef HAVE_LIBVDWXC
    call vdwxc_calculate(this%libvdwxc_ptr, cube_rho, cube_sigma, cube_dedrho, cube_dedsigma, this%energy)
#endif
    this%energy = this%energy * this%vdw_factor
    cube_dedrho = cube_dedrho * this%vdw_factor
    cube_dedsigma = cube_dedsigma * this%vdw_factor

    call fromcube(cube_dedrho, workbuffer)
    ! dedd is 1:mesh%np_part for some reason
    if(this%debug) then
      call libvdwxc_write_array(workbuffer, 'dedrho')
    end if
    dedd(1:this%mesh%np, 1) = dedd(1:this%mesh%np, 1) + workbuffer
    call fromcube(cube_dedsigma, workbuffer)
    if(this%debug) then
      call libvdwxc_write_array(workbuffer, 'dedsigma')
    end if
    do ii = 1, this%mesh%np
      dedgd(ii, :, 1) = dedgd(ii, :, 1) + M_TWO * workbuffer(ii) * gradrho(ii, :, 1)
    end do

    energy_and_integrals_buffer(1) = this%energy
    energy_and_integrals_buffer(2) = sum(rho(1:this%mesh%np,:) * dedd(1:this%mesh%np,:)) * this%mesh%volume_element
    energy_and_integrals_buffer(3) = sum(gradrho(1:this%mesh%np,:,:) * dedgd(1:this%mesh%np,:,:)) * this%mesh%volume_element

#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, energy_and_integrals_buffer, 3, MPI_DOUBLE_PRECISION, MPI_SUM, this%mesh%mpi_grp%comm, ierr)
    this%energy = energy_and_integrals_buffer(1)
#endif
    write(message(1), '(a,f18.10,a)') 'libvdwxc non-local correlation energy: ', energy_and_integrals_buffer(1), ' Ha'
    write(message(2), '(a,f18.10)')   '                      n-dedn integral: ', energy_and_integrals_buffer(2)
    write(message(3), '(a,f18.10)')   '              gradn-dedgradn integral: ', energy_and_integrals_buffer(3)
    call messages_info(3)

    SAFE_DEALLOCATE_A(workbuffer)
    SAFE_DEALLOCATE_A(cube_rho)
    SAFE_DEALLOCATE_A(cube_sigma)
    SAFE_DEALLOCATE_A(cube_dedrho)
    SAFE_DEALLOCATE_A(cube_dedsigma)
    call dcube_function_free_RS(this%cube, cf)

    POP_SUB(libvdwxc_calculate)

    contains

      subroutine tocube(array, cubearray)
        FLOAT,                 intent(in)    :: array(:)
        FLOAT,                 intent(out)   :: cubearray(:,:,:)

        PUSH_SUB(libvdwxc_calculate.tocube)

        if (this%cube%parallel_in_domains) then
          call dmesh_to_cube_parallel(this%mesh, array, this%cube, cf, this%mesh_cube_map)
        else
          if(this%mesh%parallel_in_domains) then
            call dmesh_to_cube(this%mesh, array, this%cube, cf, local = .true.)
          else
            call dmesh_to_cube(this%mesh, array, this%cube, cf)
          end if
        end if
        cubearray(:,:,:) = cf%dRS
        POP_SUB(libvdwxc_calculate.tocube)
      end subroutine tocube

      subroutine fromcube(cubearray, array)
        FLOAT,                 intent(in)  :: cubearray(:,:,:)
        FLOAT,                 intent(out) :: array(:)

        PUSH_SUB(libvdwxc_calculate.fromcube)
        cf%dRS = cubearray
        if (this%cube%parallel_in_domains) then
          call dcube_to_mesh_parallel(this%cube, cf, this%mesh, array, this%mesh_cube_map)
        else
          if(this%mesh%parallel_in_domains) then
            call dcube_to_mesh(this%cube, cf, this%mesh, array, local=.true.)
          else
            call dcube_to_mesh(this%cube, cf, this%mesh, array)
          end if
        end if
        POP_SUB(libvdwxc_calculate.fromcube)
      end subroutine fromcube

      subroutine libvdwxc_write_array(arr, fname)
        FLOAT,            intent(in) :: arr(:)
        character(len=*), intent(in) :: fname
        integer :: ierr

        call dio_function_output(OPTION__OUTPUTFORMAT__DX,  'libvdwxc-debug', &
          fname, namespace, space, this%mesh, arr, unit_one, ierr)
      end subroutine libvdwxc_write_array

    end subroutine libvdwxc_calculate

  subroutine libvdwxc_end(this)
    type(libvdwxc_t), intent(inout) :: this
    PUSH_SUB(libvdwxc_end)

#ifdef HAVE_LIBVDWXC
    call vdwxc_finalize(this%libvdwxc_ptr)
#endif

    POP_SUB(libvdwxc_end)
  end subroutine libvdwxc_end

end module libvdwxc_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
