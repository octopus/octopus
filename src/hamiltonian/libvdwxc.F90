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
      ! XXX libvdwxc actually adjusts annoying (<1e-20) values of density so inout for now.
      ! These parameters are really arrays
      real(8), intent(inout) :: rho(:,:,:), sigma(:,:,:), dedn(:,:,:), dedsigma(:,:,:), energy
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

!  interface
!    subroutine vdwxc_init_mpi(vdw, comm)
!      integer, pointer, intent(inout) :: vdw
!      integer,          intent(in)    :: comm
!    end subroutine vdwxc_init_mpi
!  end interface

!  interface
!    subroutine vdwxc_init_pfft(vdw, comm, ncpu1, ncpu2)
!      integer, pointer, intent(inout) :: vdw
!      integer,          intent(in)    :: comm
!      integer,          intent(in)    :: ncpu1
!      integer,          intent(in)    :: ncpu2
!    end subroutine vdwxc_init_pfft
!  end interface

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
    !ASSERT(.not.associated(libvdwxc%libvdwxc_ptr))
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

    integer :: fft_library

    PUSH_SUB(libvdwxc_set_geometry)
    this%mesh = mesh
    ! XXX will not work with all FFTLlibraries.  We will do FFTW(/MPI) for now.
    !call parse_variable('FFTLibrary', FFTLIB_FFTW, fft_library)
    call cube_init(this%cube, mesh%idx%ll, mesh%sb, fft_type = FFT_REAL, verbose = .true., &
      need_partition=.not.mesh%parallel_in_domains)
    fft_library = cube_getFFTLibrary(this%cube)

    if (this%cube%parallel_in_domains) then
      call mesh_cube_parallel_map_init(this%mesh_cube_map, mesh, this%cube)
    end if

    ASSERT(this%cube%parallel_in_domains .eqv. (fft_library == FFTLIB_PFFT))

#ifdef HAVE_LIBVDWXC
    ! XXX need to be sure that we get the correct physical cell
    call vdwxc_set_unit_cell(this%libvdwxc_ptr, &
      this%cube%rs_n_global(3), this%cube%rs_n_global(2), this%cube%rs_n_global(1), &
      mesh%spacing(3) * this%cube%rs_n_global(3), 0.0_8, 0.0_8, &
      0.0_8, mesh%spacing(2) * this%cube%rs_n_global(2), 0.0_8, &
      0.0_8, 0.0_8, mesh%spacing(1) * this%cube%rs_n_global(1))

    !if(mesh%parallel_in_domains) then
      !call vdwxc_init_mpi(this%libvdwxc_ptr, mesh%mpi_grp%comm)
    !else
    call vdwxc_init_serial(this%libvdwxc_ptr)
    !end if
    !call vdwxc_print(this%libvdwxc_ptr)
#else
    ASSERT(.false.)
#endif
    POP_SUB(libvdwxc_set_geometry)
  end subroutine libvdwxc_set_geometry

  subroutine libvdwxc_calculate(this, rho, gradrho, dedd, dedgd)
    type(libvdwxc_t),         intent(inout) :: this
    FLOAT, dimension(:,:),    intent(inout) :: rho !!! data type
    FLOAT, dimension(:,:,:),  intent(in)    :: gradrho
    FLOAT, dimension(:,:),    intent(inout) :: dedd
    FLOAT, dimension(:,:,:),  intent(inout) :: dedgd

    type(cube_function_t) :: rhocf, sigmacf, dedrhocf, dedsigmacf

    real(8), allocatable :: workbuffer(:)
    integer :: ii

    PUSH_SUB(libvdwxc_calculate)

    ASSERT(size(rho, 2) == 1)
    ASSERT(size(gradrho, 3) == 1)
    ASSERT(size(dedd, 2) == 1)
    ASSERT(size(dedgd, 3) == 1)

    SAFE_ALLOCATE(workbuffer(1:this%mesh%np))
    ! This is sigma, the absolute-squared density gradient:
    workbuffer(:) = sum(gradrho(:, :, 1)**2, 2)

    !call libvdwxc_print(this)

    call cube_function_null(rhocf)
    call cube_function_null(sigmacf)
    call cube_function_null(dedrhocf)
    call cube_function_null(dedsigmacf)

    call dcube_function_alloc_RS(this%cube, rhocf, in_device = .false.)
    call dcube_function_alloc_RS(this%cube, sigmacf, in_device = .false.)
    call dcube_function_alloc_RS(this%cube, dedrhocf, in_device = .false.)
    call dcube_function_alloc_RS(this%cube, dedsigmacf, in_device = .false.)

    call tocube(rho(:, 1), rhocf)
    call tocube(workbuffer, sigmacf)

    this%energy = M_ZERO
#ifdef HAVE_LIBVDWXC
    call vdwxc_calculate(this%libvdwxc_ptr, rhocf%dRS, sigmacf%dRS, dedrhocf%dRS, dedsigmacf%dRS, this%energy)
#endif
    ! XXXXXXX energy may require MPI sum

    call fromcube(dedrhocf, workbuffer)
    ! dedd is 1:mesh%np_part for some reason
    dedd(1:this%mesh%np, 1) = dedd(1:this%mesh%np, 1) + workbuffer
    call fromcube(dedsigmacf, workbuffer)
    do ii=1, this%mesh%np
      dedgd(ii, :, 1) = dedgd(ii, :, 1) + M_TWO * workbuffer(ii) * gradrho(ii, :, 1)
    end do

    SAFE_DEALLOCATE_A(workbuffer)
    call dcube_function_free_RS(this%cube, rhocf)
    call dcube_function_free_RS(this%cube, sigmacf)
    call dcube_function_free_RS(this%cube, dedrhocf)
    call dcube_function_free_RS(this%cube, dedsigmacf)

    POP_SUB(libvdwxc_calculate)

    contains

      subroutine tocube(array, cf)
        FLOAT,                 intent(in)    :: array(:)
        type(cube_function_t), intent(inout) :: cf

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
        POP_SUB(libvdwxc_calculate.tocube)
      end subroutine tocube

      subroutine fromcube(cf, array)
        type(cube_function_t), intent(in)  :: cf
        FLOAT,                 intent(out) :: array(:)

        PUSH_SUB(libvdwxc_calculate.fromcube)
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
