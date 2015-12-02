#include "global.h"

! Interface to libvdwxc.
! Naming:
!  * Functions that start with libvdwxc_ are public, to be called from other parts of Octopus.
!  * Interfaces that start with vdwxc_ are actual functions of libvdwxc.

module libvdwxc_m
  use cube_m
  use cube_function_m
  use derivatives_m
  use fft_m
  use global_m
  use grid_m
  use mesh_m
  use mesh_cube_parallel_map_m
  use messages_m
  use mpi_m
  use parser_m
  use pfft_m
  use profiling_m

  implicit none

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
    integer, pointer               :: libvdwxc_ptr
    type(mesh_t)                   :: mesh
    type(cube_t)                   :: cube
    type(mesh_cube_parallel_map_t) :: mesh_cube_map
    integer                        :: functional
    FLOAT                          :: energy
  end type libvdwxc_t

#ifdef HAVE_LIBVDWXC
  ! These interfaces correspond to functions in libvdwxc
  interface
    subroutine vdwxc_new(functional, vdw)
      integer,          intent(in)  :: functional
      integer, pointer, intent(out) :: vdw
    end subroutine vdwxc_new
  end interface

  interface
    subroutine vdwxc_print(vdw)
      integer, pointer, intent(in) :: vdw
    end subroutine vdwxc_print
  end interface

  interface
    subroutine vdwxc_calculate(vdw, rho, sigma, dedn, dedsigma, energy)
      integer, pointer, intent(inout) :: vdw
      real(8),          intent(in)    :: rho(:,:,:)
      real(8),          intent(in)    :: sigma(:,:,:)
      real(8),          intent(inout) :: dedn(:,:,:)
      real(8),          intent(inout) :: dedsigma(:,:,:)
      real(8),          intent(inout) :: energy
      !real(8), intent(inout) :: rho, sigma, dedn, dedsigma, energy
    end subroutine vdwxc_calculate
  end interface

  interface
    subroutine vdwxc_set_unit_cell(vdw, nx, ny, nz, C00, C01, C02, &
      C10, C11, C12, C20, C21, C22)
      integer, pointer, intent(inout) :: vdw
      integer,          intent(in)    :: nx, ny, nz
      real(8),          intent(in)    :: C00, C01, C02, C10, C11, C12, C20, C21, C22
    end subroutine vdwxc_set_unit_cell
  end interface

  interface
    subroutine vdwxc_init_serial(vdw)
      integer, pointer, intent(inout) :: vdw
    end subroutine vdwxc_init_serial
  end interface

#ifdef HAVE_MPI
  interface
    subroutine vdwxc_init_mpi(vdw, comm)
      integer, pointer, intent(inout) :: vdw
      integer,          intent(in)    :: comm
    end subroutine vdwxc_init_mpi
  end interface

  interface
    subroutine vdwxc_init_pfft(vdw, comm, ncpu1, ncpu2)
      integer, pointer, intent(inout) :: vdw
      integer,          intent(in)    :: comm
      integer,          intent(in)    :: ncpu1
      integer,          intent(in)    :: ncpu2
    end subroutine vdwxc_init_pfft
  end interface
#endif

  interface
    subroutine vdwxc_finalize(vdw)
      integer, pointer, intent(inout) :: vdw
    end subroutine vdwxc_finalize
  end interface
#endif

contains

  subroutine libvdwxc_init(libvdwxc, functional)
    type(libvdwxc_t), intent(out) :: libvdwxc
    integer,          intent(in)  :: functional

    PUSH_SUB(libvdwxc_init)
#ifdef HAVE_LIBVDWXC
    call vdwxc_new(functional, libvdwxc%libvdwxc_ptr)
#else
    message(1) = "Octopus not compiled with libvdwxc"
    call messages_fatal(1)
#endif
    ASSERT(associated(libvdwxc%libvdwxc_ptr))
    libvdwxc%functional = functional
    POP_SUB(libvdwxc_init)
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
      write(message(2), '(4x,a)') 'vdW-DF-CX from libvdwxc'
    else
      write(message(2), '(4x,a)') 'unknown libvdwxc functional'
    end if
    call messages_info(2, iunit)
    POP_SUB(libvdwxc_write_info)
  end subroutine libvdwxc_write_info

  subroutine libvdwxc_set_geometry(this, mesh)
    type(libvdwxc_t), intent(inout) :: this
    type(mesh_t),     intent(inout)    :: mesh

    integer :: pfftx, pffty
    integer :: adjusted_leading_dim
    integer :: fft_lib
    logical :: use_pfft

    PUSH_SUB(libvdwxc_set_geometry)
    this%mesh = mesh

    ! libvdwxc can use either FFTW-MPI or PFFT and Octopus should not
    ! care too much, as long as this particular cube is properly
    ! configured.  Unfortunately most of the args to cube_init are
    ! very entangled with FFT libraries.  All we want and need is to
    ! specify our own decomposition and pass that to libvdwxc, but
    ! we can only say "we want to use PFFT" and then it does the
    ! decomposition.

    call parse_variable('FFTLibrary', FFTLIB_NONE, fft_lib)

    if(fft_lib == FFTLIB_PFFT) then
      use_pfft = .true.
    else if(fft_lib == FFTLIB_NONE) then
      use_pfft = mesh%parallel_in_domains
    else
      write(message(1), '(a)') 'libvdwxc/fftlib conflict'
      call messages_fatal(1)
    end if

    if(use_pfft) then
      ! For parallel libvdwxc one needs PFFT as of now.
      ! Passing FFTLIB_PFFT causes reasonable error if library is not present
      call cube_init(this%cube, mesh%idx%ll, mesh%sb, mpi_grp = mesh%mpi_grp, &
        fft_type = FFT_REAL, fft_library = FFTLIB_PFFT, need_partition = .true.)
      call mesh_cube_parallel_map_init(this%mesh_cube_map, mesh, this%cube)
      ASSERT(this%cube%parallel_in_domains)
    else
      call cube_init(this%cube, mesh%idx%ll, mesh%sb)
    end if

    ! There is some low-level implementation issue where with PFFT,
    ! the leading dimension is one smaller than normal for some reason.
    ! We therefore use rs_n as the global size.  (We never parallelize
    ! over the leading dimension anyway)
    adjusted_leading_dim = this%cube%rs_n(1)

#ifdef HAVE_LIBVDWXC
    ! XXX need to be sure that we get the correct physical cell
    call vdwxc_set_unit_cell(this%libvdwxc_ptr, &
      this%cube%rs_n_global(3), this%cube%rs_n_global(2), adjusted_leading_dim, &
      mesh%spacing(3) * this%cube%rs_n_global(3), 0.0_8, 0.0_8, &
      0.0_8, mesh%spacing(2) * this%cube%rs_n_global(2), 0.0_8, &
      0.0_8, 0.0_8, mesh%spacing(1) * this%cube%rs_n_global(1))

    if(use_pfft) then
#ifdef HAVE_MPI
#ifdef HAVE_PFFT
      ASSERT(adjusted_leading_dim == this%cube%rs_n_global(1) - 1)
      call pfft_decompose(mesh%mpi_grp%size, pfftx, pffty)
      call vdwxc_init_pfft(this%libvdwxc_ptr, mesh%mpi_grp%comm, pfftx, pffty)
#endif
#endif
    else
      ASSERT(adjusted_leading_dim == this%cube%rs_n_global(1))
      call vdwxc_init_serial(this%libvdwxc_ptr)
    end if
    call vdwxc_print(this%libvdwxc_ptr)
#endif
    POP_SUB(libvdwxc_set_geometry)
  end subroutine libvdwxc_set_geometry

  subroutine libvdwxc_calculate(this, rho, gradrho, dedd, dedgd)
    type(libvdwxc_t),         intent(inout) :: this
    FLOAT, dimension(:,:),    intent(inout) :: rho !!! data type
    FLOAT, dimension(:,:,:),  intent(in)    :: gradrho
    FLOAT, dimension(:,:),    intent(inout) :: dedd
    FLOAT, dimension(:,:,:),  intent(inout) :: dedgd

    type(cube_function_t) :: cf!rhocf, sigmacf, dedrhocf, dedsigmacf

    real(8), allocatable :: workbuffer(:)
    real(8), allocatable :: cube_rho(:,:,:), cube_sigma(:,:,:), cube_dedrho(:,:,:), cube_dedsigma(:,:,:)
    real(8), dimension(1) :: tmp_energy
    integer :: ii, ierr, magic

    PUSH_SUB(libvdwxc_calculate)

    ASSERT(size(rho, 2) == 1)
    ASSERT(size(gradrho, 3) == 1)
    ASSERT(size(dedd, 2) == 1)
    ASSERT(size(dedgd, 3) == 1)

    ! Well.  I thought we would be using four different cube functions
    ! for rho, sigma, dedrho and dedsigma.
    !
    ! But since the cube has PFFT associated, any attempt to use it
    ! results in some kind of behind-the-scenes use of a possibly
    ! global FFT-related buffer.
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
    magic = this%cube%rs_n(1)! + 1
    SAFE_ALLOCATE(cube_rho(1:magic, 1:this%cube%rs_n(2), 1:this%cube%rs_n(3)))
    SAFE_ALLOCATE(cube_sigma(1:magic, 1:this%cube%rs_n(2), 1:this%cube%rs_n(3)))
    SAFE_ALLOCATE(cube_dedrho(1:magic, 1:this%cube%rs_n(2), 1:this%cube%rs_n(3)))
    SAFE_ALLOCATE(cube_dedsigma(1:magic, 1:this%cube%rs_n(2), 1:this%cube%rs_n(3)))
    ! This is sigma, the absolute-squared density gradient:
    workbuffer(:) = sum(gradrho(:, :, 1)**2, 2)

    cube_rho = M_ZERO
    cube_sigma = M_ZERO
    cube_dedrho = M_ZERO
    cube_dedsigma = M_ZERO

    call cube_function_null(cf)
    call dcube_function_alloc_RS(this%cube, cf, in_device = .false.)

    call tocube(rho(:, 1), cube_rho)
    call tocube(workbuffer, cube_sigma)

    this%energy = M_ZERO
#ifdef HAVE_LIBVDWXC
    !print*, 'inp sum', sum(cube_rho), sum(cube_sigma)
    call vdwxc_calculate(this%libvdwxc_ptr, cube_rho, cube_sigma, cube_dedrho, cube_dedsigma, this%energy)
    !print*, 'pot sum', sum(cube_dedrho), sum(cube_dedsigma)
    !print*, 'energy', this%energy
#endif
#ifdef HAVE_MPI
    !if(this%cube%parallel_in_domains) then
      !tmp_energy(1) = this%energy
      !call MPI_Allreduce(MPI_IN_PLACE, tmp_energy, 1, MPI_FLOAT, MPI_SUM, this%mesh%mpi_grp%comm, ierr) ! XXX check ierr?
      !this%energy = tmp_energy(1)
    !end if
#endif
    ! XXXXXXX energy may require MPI sum
    call fromcube(cube_dedrho, workbuffer)
    ! dedd is 1:mesh%np_part for some reason
    dedd(1:this%mesh%np, 1) = dedd(1:this%mesh%np, 1) + workbuffer
    call fromcube(cube_dedsigma, workbuffer)
    do ii=1, this%mesh%np
      dedgd(ii, :, 1) = dedgd(ii, :, 1) + M_TWO * workbuffer(ii) * gradrho(ii, :, 1)
    end do

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
        !type(cube_function_t), intent(inout) :: cf

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
        cubearray(1:magic,:,:) = cf%dRS
        POP_SUB(libvdwxc_calculate.tocube)
      end subroutine tocube

      subroutine fromcube(cubearray, array)
        FLOAT,                 intent(in)  :: cubearray(:,:,:)
        !type(cube_function_t), intent(in)  :: cf
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

  end subroutine libvdwxc_calculate

  subroutine libvdwxc_end(this)
    type(libvdwxc_t), intent(inout) :: this
    PUSH_SUB(libvdwxc_end)
#ifdef HAVE_LIBVDWXC
    call vdwxc_finalize(this%libvdwxc_ptr)
#endif
    POP_SUB(libvdwxc_end)
  end subroutine libvdwxc_end

end module libvdwxc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
