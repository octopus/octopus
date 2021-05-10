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

#include "global.h"

module poisson_fft_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use fft_oct_m
  use fourier_space_oct_m
  use global_oct_m
  use lattice_vectors_oct_m
  use loct_math_oct_m
  use math_oct_m
  use mesh_cube_parallel_map_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_cutoff_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use splines_oct_m
  use submesh_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none
  private
  public ::                  &
    poisson_fft_t,           &
    poisson_fft_init,        &
    poisson_fft_get_kernel,  &
    poisson_fft_end,         &
    dpoisson_fft_solve,      &
    zpoisson_fft_solve

  integer, public, parameter ::                &
       POISSON_FFT_KERNEL_NONE      = -1,      &
       POISSON_FFT_KERNEL_SPH       =  0,      &
       POISSON_FFT_KERNEL_CYL       =  1,      &
       POISSON_FFT_KERNEL_PLA       =  2,      &
       POISSON_FFT_KERNEL_NOCUT     =  3,      &
       POISSON_FFT_KERNEL_CORRECTED =  4,      &
       POISSON_FFT_KERNEL_HOCKNEY   =  5

  type poisson_fft_t
    ! Components are public by default
    type(fourier_space_op_t) :: coulb  !< object for Fourier space operations
    integer                  :: kernel !< choice of kernel, one of options above
    FLOAT                    :: soft_coulb_param !< Soft-Coulomb parameter
  end type poisson_fft_t
contains

  subroutine poisson_fft_init(this, namespace, space, mesh, cube, kernel, soft_coulb_param, fullcube)
    type(poisson_fft_t), intent(out)   :: this
    type(namespace_t),   intent(in)    :: namespace
    type(space_t),       intent(in)    :: space
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube
    integer,             intent(in)    :: kernel
    FLOAT, optional,     intent(in)    :: soft_coulb_param
    type(cube_t), optional, intent(in) :: fullcube !< needed for Hockney kernel

    PUSH_SUB(poisson_fft_init)

    this%kernel = kernel
    this%soft_coulb_param = optional_default(soft_coulb_param, M_ZERO)

    this%coulb%qq(1:space%periodic_dim) = M_ZERO
    this%coulb%singularity = M_ZERO
    this%coulb%mu = M_ZERO

    call poisson_fft_get_kernel(namespace, space, mesh, cube, this%coulb, kernel, soft_coulb_param, fullcube) 

    POP_SUB(poisson_fft_init)
  end subroutine poisson_fft_init

  subroutine poisson_fft_get_kernel(namespace, space, mesh, cube, coulb, kernel, soft_coulb_param, fullcube)
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb
    integer,                  intent(in)    :: kernel
    FLOAT,        optional,   intent(in)    :: soft_coulb_param
    type(cube_t), optional,   intent(in)    :: fullcube !< needed for Hockney kernel 

    PUSH_SUB(poisson_fft_get_kernel)

    if(coulb%mu > M_EPSILON) then
      if (space%dim /= 3 .or. kernel /= POISSON_FFT_KERNEL_NOCUT) then
        message(1) = "The screened Coulomb potential is only implemented in 3D for PoissonFFTKernel=fft_nocut."
        call messages_fatal(1)
      end if
    end if


    if(kernel == POISSON_FFT_KERNEL_HOCKNEY) then
      if (.not. present(fullcube)) then
        message(1) = "Hockney's FFT-kernel needs cube of full unit cell "
        call messages_fatal(1)
      else
        if (.not. allocated(fullcube%fft)) then
          message(1) = "Hockney's FFT-kernel needs PoissonSolver=fft"
          call messages_fatal(1)
        end if
      end if
    end if


    select case (space%dim)
    case(1)
      ASSERT(present(soft_coulb_param))
      select case(kernel)
      case(POISSON_FFT_KERNEL_SPH)
        call poisson_fft_build_1d_0d(namespace, mesh, cube, coulb, soft_coulb_param)
      case(POISSON_FFT_KERNEL_NOCUT)
        call poisson_fft_build_1d_1d(mesh, cube, coulb, soft_coulb_param)
      case default
        message(1) = "Invalid Poisson FFT kernel for 1D."
        call messages_fatal(1)
      end select

    case(2)
      select case(kernel)
      case(POISSON_FFT_KERNEL_SPH)
        call poisson_fft_build_2d_0d(namespace, mesh, cube, coulb)
      case(POISSON_FFT_KERNEL_CYL)
        call poisson_fft_build_2d_1d(namespace, mesh, cube, coulb)
      case(POISSON_FFT_KERNEL_NOCUT)
        call poisson_fft_build_2d_2d(mesh, cube, coulb)
      case default
        message(1) = "Invalid Poisson FFT kernel for 2D."
        call messages_fatal(1)
      end select

    case(3)
      select case(kernel)
      case(POISSON_FFT_KERNEL_SPH, POISSON_FFT_KERNEL_CORRECTED)
        call poisson_fft_build_3d_0d(namespace,  mesh, cube, kernel, coulb)

      case(POISSON_FFT_KERNEL_CYL)
        call poisson_fft_build_3d_1d(namespace, space, mesh, cube, coulb)

      case(POISSON_FFT_KERNEL_PLA)
        call poisson_fft_build_3d_2d(namespace, mesh, cube, coulb)

      case(POISSON_FFT_KERNEL_NOCUT)
        call poisson_fft_build_3d_3d(mesh, cube, coulb)

      case(POISSON_FFT_KERNEL_HOCKNEY)
        call poisson_fft_build_3d_3d_hockney(mesh, cube, coulb, fullcube)

      case default
        message(1) = "Invalid Poisson FFT kernel for 3D."
        call messages_fatal(1)
      end select
    end select

    POP_SUB(poisson_fft_get_kernel)
  end subroutine poisson_fft_get_kernel

  !-----------------------------------------------------------------

  subroutine get_cutoff(namespace, default_r_c, r_c)
    type(namespace_t),   intent(in)  :: namespace
    FLOAT,               intent(in)  :: default_r_c
    FLOAT,               intent(out) :: r_c

    PUSH_SUB(get_cutoff)

    call parse_variable(namespace, 'PoissonCutoffRadius', default_r_c, r_c, units_inp%length)

    call messages_write('Info: Poisson Cutoff Radius     =')
    call messages_write(r_c, units = units_out%length, fmt = '(f6.1)')
    call messages_info()

    if ( r_c > default_r_c + M_EPSILON) then
      call messages_write('Poisson cutoff radius is larger than cell size.', new_line = .true.)
      call messages_write('You can see electrons in neighboring cell(s).')
      call messages_warning()
    end if

    POP_SUB(get_cutoff)
  end subroutine get_cutoff

  !-----------------------------------------------------------------
  subroutine poisson_fft_gg_transform(gg_in, temp, periodic_dim, latt, qq, gg, modg2)
    integer,                 intent(in)    :: gg_in(:)
    FLOAT,                   intent(in)    :: temp(:)
    integer,                 intent(in)    :: periodic_dim
    type(lattice_vectors_t), intent(in)    :: latt
    FLOAT,                   intent(in)    :: qq(:)
    FLOAT,                   intent(inout) :: gg(:)
    FLOAT,                   intent(out)   :: modg2

!    integer :: idir

    ! no PUSH_SUB, called too frequently

    gg(1:3) = gg_in(1:3)
    gg(1:periodic_dim) = gg(1:periodic_dim) + qq(1:periodic_dim)
    gg(1:3) = gg(1:3) * temp(1:3)
    gg(1:3) = matmul(latt%klattice_primitive(1:3,1:3),gg(1:3))
! MJV 27 jan 2015 this should not be necessary
!    do idir = 1, 3
!      gg(idir) = gg(idir) / lalg_nrm2(3, latt%klattice_primitive(1:3, idir))
!    end do

    modg2 = sum(gg(1:3)**2)

  end subroutine poisson_fft_gg_transform

  !-----------------------------------------------------------------
  subroutine poisson_fft_build_3d_3d(mesh, cube, coulb)
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb

    integer :: ix, iy, iz, ixx(3), db(3), n1, n2, n3, lx, ly, lz
    FLOAT :: temp(3), modg2, inv_four_mu2
    FLOAT :: gg(3)
    FLOAT, allocatable :: fft_Coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_3d)

    db(1:3) = cube%rs_n_global(1:3)

    if(coulb%mu > M_EPSILON) then
      inv_four_mu2 = M_ONE/((M_TWO*coulb%mu)**2)
    end if

    n1 = max(1, cube%fs_n(1))
    n2 = max(1, cube%fs_n(2))
    n3 = max(1, cube%fs_n(3))
    ! store the Fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:n1, 1:n2, 1:n3))
    fft_Coulb_FS = M_ZERO

    temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))

    do lx = 1, n1
      ix = cube%fs_istart(1) + lx - 1
      ixx(1) = pad_feq(ix, db(1), .true.)
      do ly = 1, n2
        iy = cube%fs_istart(2) + ly - 1
        ixx(2) = pad_feq(iy, db(2), .true.)
        do lz = 1, n3
          iz = cube%fs_istart(3) + lz - 1
          ixx(3) = pad_feq(iz, db(3), .true.)

         call poisson_fft_gg_transform(ixx, temp, 3, mesh%sb%latt, coulb%qq, gg, modg2)

         !HH not very elegant
         if(cube%fft%library.eq.FFTLIB_NFFT) modg2=cube%Lfs(ix,1)**2+cube%Lfs(iy,2)**2+cube%Lfs(iz,3)**2

         if(abs(modg2) > CNST(1e-6)) then
           !Screened coulomb potential (erfc function)
           if(coulb%mu > M_EPSILON) then
             fft_Coulb_FS(lx, ly, lz) = M_FOUR*M_PI/modg2*(M_ONE-exp(-modg2*inv_four_mu2))
           else
             fft_Coulb_FS(lx, ly, lz) = M_FOUR*M_PI/modg2
           end if
         else
           !Screened coulomb potential (erfc function)
           if(coulb%mu > M_EPSILON) then
             !Analytical limit of 1/|q|^2*(1-exp(-|q|^2/4mu^2))
             fft_Coulb_FS(lx, ly, lz) =  M_FOUR*M_PI*inv_four_mu2
           else
             !We use the user-defined value of the singularity
             fft_Coulb_FS(lx, ly, lz) = coulb%singularity
           end if
         end if

        end do
      end do

    end do

    call dfourier_space_op_init(coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_3d_3d)
  end subroutine poisson_fft_build_3d_3d
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  !> Kernel for Hockneys algorithm that solves the poisson equation
  !! in a small box while respecting the periodicity of a larger box
  !! A. Damle, L. Lin, L. Ying, JCTC, 2015
  !! DOI: 10.1021/ct500985f, supplementary info  
  subroutine poisson_fft_build_3d_3d_hockney(mesh, cube, coulb, fullcube)
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb
    type(cube_t),             intent(in)    :: fullcube

    integer :: ix, iy, iz, ixx(3), db(3), nfs(3), nrs(3), nfs_s(3), nrs_s(3)
    FLOAT :: temp(3), modg2
    FLOAT :: gg(3)
    FLOAT, allocatable :: fft_Coulb_small_RS(:,:,:)
    FLOAT, allocatable :: fft_Coulb_RS(:,:,:)
    CMPLX, allocatable :: fft_Coulb_small_FS(:,:,:), fft_Coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_3d_hockney)

    ! dimensions of large boxes
    nfs(1:3) = fullcube%fs_n_global(1:3)
    nrs(1:3) = fullcube%rs_n_global(1:3)

    SAFE_ALLOCATE(fft_Coulb_FS(1:nfs(1),1:nfs(2),1:nfs(3)))
    SAFE_ALLOCATE(fft_Coulb_RS(1:nrs(1),1:nrs(2),1:nrs(3)))

    ! dimensions of small boxes x_s
    nfs_s(1:3) = cube%fs_n_global(1:3)
    nrs_s(1:3) = cube%rs_n_global(1:3)

    SAFE_ALLOCATE(fft_Coulb_small_FS(1:nfs_s(1),1:nfs_s(2),1:nfs_s(3)))
    SAFE_ALLOCATE(fft_Coulb_small_RS(1:nrs_s(1),1:nrs_s(2),1:nrs_s(3)))

    ! build full periodic Coulomb potenital in Fourier space
    fft_Coulb_FS = M_ZERO

    db(1:3) = fullcube%rs_n_global(1:3)
    temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))
    
    do ix = 1, nfs(1)
      ixx(1) = pad_feq(ix, db(1), .true.)
      do iy = 1, nfs(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, nfs(3)
          ixx(3) = pad_feq(iz, db(3), .true.)
          
          call poisson_fft_gg_transform(ixx, temp, 3, mesh%sb%latt, coulb%qq, gg, modg2)
          
          if(abs(modg2) > M_EPSILON) then
            fft_Coulb_FS(ix, iy, iz) = M_ONE/modg2
          else
            fft_Coulb_FS(ix, iy, iz) = M_ZERO
          end if
        end do
      end do
    end do
    
    do iz = 1, nfs(3)
      do iy = 1, nfs(2)
        do ix = 1, nfs(1)
          fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
        end do
      end do
    end do

    ! get periodic Coulomb potential in real space
    call dfft_backward(fullcube%fft, fft_Coulb_FS, fft_Coulb_RS)

    !  copy to small box by respecting this pattern
    !  full periodic coulomb: |abc--------------------------xyz|
    !                Hockney: |abcxyz|
    do ix = 1, nrs_s(1)
      if(ix.le.nrs_s(1)/2+1) then
        ixx(1)=ix
      else
        ixx(1) = nrs(1) - nrs_s(1)+ix
      end if
      do iy = 1, nrs_s(2)
        if(iy.le.nrs_s(2)/2+1) then
          ixx(2)=iy
        else
          ixx(2) = nrs(2) - nrs_s(2)+iy
        end if
        do iz = 1, nrs_s(3)
          if(iz.le.nrs_s(3)/2+1) then
            ixx(3)=iz
          else
            ixx(3) = nrs(3) - nrs_s(3)+iz
          end if
          fft_Coulb_small_RS(ix,iy,iz) = fft_Coulb_RS(ixx(1),ixx(2),ixx(3))
        end do
      end do
    end do
    ! make Hockney kernel in Fourier space
    call dfft_forward(cube%fft, fft_Coulb_small_RS, fft_Coulb_small_FS)
    !dummy copy for type conversion
    fft_Coulb_small_RS(1:nfs_s(1),1:nfs_s(2),1:nfs_s(3)) = &
                                                 TOFLOAT( fft_Coulb_small_FS(1:nfs_s(1),1:nfs_s(2),1:nfs_s(3)))


    ! Restrict array to local part to support pfft
    ! For FFTW this reduces simply to the full array
    call dfourier_space_op_init(coulb, cube, &
      fft_Coulb_small_RS(cube%fs_istart(1):cube%fs_istart(1)+cube%fs_n(1), &
                         cube%fs_istart(2):cube%fs_istart(2)+cube%fs_n(2), &
                         cube%fs_istart(3):cube%fs_istart(3)+cube%fs_n(3)))

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    SAFE_DEALLOCATE_A(fft_Coulb_RS)
    SAFE_DEALLOCATE_A(fft_Coulb_small_FS)
    SAFE_DEALLOCATE_A(fft_Coulb_small_RS)

    POP_SUB(poisson_fft_build_3d_3d_hockney)

  end subroutine poisson_fft_build_3d_3d_hockney

  !-----------------------------------------------------------------
  !> C. A. Rozzi et al., Phys. Rev. B 73, 205119 (2006), Table I
  subroutine poisson_fft_build_3d_2d(namespace, mesh, cube, coulb)
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb

    integer :: ix, iy, iz, ixx(3), db(3)
    integer :: lx, ly, lz, n1, n2, n3
    FLOAT :: temp(3), modg2
    FLOAT :: gpar, gz, r_c, gg(3), default_r_c
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_2d)

    db(1:3) = cube%rs_n_global(1:3)

    !%Variable PoissonCutoffRadius
    !%Type float
    !%Section Hamiltonian::Poisson
    !%Description
    !% When <tt>PoissonSolver = fft</tt> and <tt>PoissonFFTKernel</tt> is neither <tt>multipole_corrections</tt>
    !% nor <tt>fft_nocut</tt>,
    !% this variable controls the distance after which the electron-electron interaction goes to zero.
    !% A warning will be written if the value is too large and will cause spurious interactions between images.
    !% The default is half of the FFT box max dimension in a finite direction.
    !%End

    default_r_c = db(3)*mesh%spacing(3)/M_TWO
    call get_cutoff(namespace, default_r_c, r_c)

    n1 = max(1, cube%fs_n(1))
    n2 = max(1, cube%fs_n(2))
    n3 = max(1, cube%fs_n(3))
    ! store the Fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:n1, 1:n2, 1:n3))
    fft_Coulb_FS = M_ZERO

    temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))

    do lx = 1, n1
      ix = cube%fs_istart(1) + lx - 1
      ixx(1) = pad_feq(ix, db(1), .true.)
      do ly = 1, n2
        iy = cube%fs_istart(2) + ly - 1
        ixx(2) = pad_feq(iy, db(2), .true.)
        do lz = 1, n3
          iz = cube%fs_istart(3) + lz - 1
          ixx(3) = pad_feq(iz, db(3), .true.)

          call poisson_fft_gg_transform(ixx, temp, 2, mesh%sb%latt, coulb%qq, gg, modg2)

          if(abs(modg2) > M_EPSILON) then
            gz = abs(gg(3))
            gpar = hypot(gg(1), gg(2))
            ! note: if gpar = 0, then modg2 = gz**2
            fft_Coulb_FS(lx, ly, lz) = poisson_cutoff_3D_2D(gpar,gz,r_c)/modg2
          else
            fft_Coulb_FS(lx, ly, lz) = -M_HALF*r_c**2
          end if
          fft_Coulb_FS(lx, ly, lz) = M_FOUR*M_PI*fft_Coulb_FS(lx, ly, lz)
        end do
      end do

    end do

    call dfourier_space_op_init(coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_3d_2d)
  end subroutine poisson_fft_build_3d_2d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> C. A. Rozzi et al., Phys. Rev. B 73, 205119 (2006), Table I
  subroutine poisson_fft_build_3d_1d(namespace, space, mesh, cube, coulb)
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb

    type(spline_t)     :: cylinder_cutoff_f
    FLOAT, allocatable :: x(:), y(:)
    integer :: ix, iy, iz, ixx(3), db(3), k, ngp
    integer :: lx, ly, lz, n1, n2, n3, lxx(3)
    FLOAT :: temp(3), modg2, xmax
    FLOAT :: gperp, gx, gy, gz, r_c, gg(3), default_r_c
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_1d)

    db(1:3) = cube%rs_n_global(1:3)

    default_r_c = maxval(db(2:3)*mesh%spacing(2:3)/M_TWO)
    call get_cutoff(namespace, default_r_c, r_c)

    n1 = max(1, cube%fs_n(1))
    n2 = max(1, cube%fs_n(2))
    n3 = max(1, cube%fs_n(3))
    ! store the Fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:n1, 1:n2, 1:n3))
    fft_Coulb_FS = M_ZERO

    temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))

    if (.not. space%is_periodic()) then
      ngp = 8*db(2)
      SAFE_ALLOCATE(x(1:ngp))
      SAFE_ALLOCATE(y(1:ngp))
    end if


    do lx = 1, n1
      ix = cube%fs_istart(1) + lx - 1
      ixx(1) = pad_feq(ix, db(1), .true.)
      lxx(1) = ixx(1) - cube%fs_istart(1) + 1
      gx = temp(1)*ixx(1)

      if (.not. space%is_periodic()) then
        call spline_init(cylinder_cutoff_f)
        xmax = sqrt((temp(2)*db(2)/2)**2 + (temp(3)*db(3)/2)**2)
        do k = 1, ngp
          x(k) = (k-1)*(xmax/(ngp-1))
          y(k) = poisson_cutoff_3D_1D_finite(gx, x(k), M_TWO*mesh%sb%lsize(1), maxval(M_TWO*mesh%sb%lsize(2:3)))
        end do
        call spline_fit(ngp, x, y, cylinder_cutoff_f)
      end if

      do ly = 1, n2
        iy = cube%fs_istart(2) + ly - 1
        ixx(2) = pad_feq(iy, db(2), .true.)
        lxx(2) = ixx(2) - cube%fs_istart(2) + 1
        do lz = 1, n3
          iz = cube%fs_istart(3) + lz - 1
          ixx(3) = pad_feq(iz, db(3), .true.)
          lxx(3) = ixx(3) - cube%fs_istart(3) + 1

          call poisson_fft_gg_transform(ixx, temp, 1, mesh%sb%latt, coulb%qq, gg, modg2)

          if(abs(modg2) > M_EPSILON) then
            gperp = hypot(gg(2), gg(3))
            if (space%periodic_dim == 1) then
              fft_Coulb_FS(lx, ly, lz) = poisson_cutoff_3D_1D(abs(gx), gperp, r_c)/modg2
            else if (.not. space%is_periodic()) then
              gy = gg(2)
              gz = gg(3)
              if ((gz >= M_ZERO) .and. (gy >= M_ZERO)) then
                fft_Coulb_FS(lx, ly, lz) = spline_eval(cylinder_cutoff_f, gperp)
              end if
              if ((gz >= M_ZERO) .and. (gy < M_ZERO)) then
                fft_Coulb_FS(lx, ly, lz) = fft_Coulb_FS(lx, -lxx(2) + 1, lz)
              end if
              if ((gz < M_ZERO) .and. (gy >= M_ZERO)) then
                fft_Coulb_FS(lx, ly, lz) = fft_Coulb_FS(lx, ly, -lxx(3) + 1)
              end if
              if ((gz < M_ZERO) .and. (gy < M_ZERO) ) then
                fft_Coulb_FS(lx, ly, lz) = fft_Coulb_FS(lx, -lxx(2) + 1, -lxx(3) + 1)
              end if
            end if

          else
            if (space%periodic_dim == 1) then
              fft_Coulb_FS(lx, ly, lz) = -(M_HALF*log(r_c) - M_FOURTH)*r_c**2
            else if (.not. space%is_periodic()) then
              fft_Coulb_FS(lx, ly, lz) = poisson_cutoff_3D_1D_finite(M_ZERO, M_ZERO, &
                M_TWO*mesh%sb%lsize(1), maxval(M_TWO*mesh%sb%lsize(2:3)))
            end if

          end if
          fft_Coulb_FS(lx, ly, lz) = M_FOUR*M_PI*fft_Coulb_FS(lx, ly, lz)
        end do
      end do

      if (.not. space%is_periodic()) then
        call spline_end(cylinder_cutoff_f)
      end if
    end do

    call dfourier_space_op_init(coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)
    POP_SUB(poisson_fft_build_3d_1d)
  end subroutine poisson_fft_build_3d_1d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> C. A. Rozzi et al., Phys. Rev. B 73, 205119 (2006), Table I
  subroutine poisson_fft_build_3d_0d(namespace, mesh, cube, kernel, coulb)
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    integer,                  intent(in)    :: kernel
    type(fourier_space_op_t), intent(inout) :: coulb

    integer :: ix, iy, iz, ixx(3), db(3), lx, ly, lz, n1, n2, n3
    FLOAT :: temp(3), modg2
    FLOAT :: r_c, gg(3), default_r_c
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_0d)

    db(1:3) = cube%rs_n_global(1:3)

    if (kernel /= POISSON_FFT_KERNEL_CORRECTED) then
      default_r_c = M_ZERO
      do ix = 1, 3  
        temp(1:3) = M_ZERO
        temp(ix) = db(ix)*mesh%spacing(ix)/M_TWO
        temp(1:3) = matmul(mesh%sb%latt%klattice_primitive(1:3,1:3),temp(1:3)) 
        default_r_c = maxval(temp(1:3))
     end do
     call get_cutoff(namespace, default_r_c, r_c)
    end if

    n1 = max(1, cube%fs_n(1))
    n2 = max(1, cube%fs_n(2))
    n3 = max(1, cube%fs_n(3))

    ! store the fourier transform of the Coulomb interaction
    ! store only the relevant part if PFFT is used
    SAFE_ALLOCATE(fft_Coulb_FS(1:n1,1:n2,1:n3))
    fft_Coulb_FS = M_ZERO

    temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))
    do lx = 1, n1
      ix = cube%fs_istart(1) + lx - 1
      ixx(1) = pad_feq(ix, db(1), .true.)
      do ly = 1, n2
        iy = cube%fs_istart(2) + ly - 1
        ixx(2) = pad_feq(iy, db(2), .true.)
        do lz = 1, n3
          iz = cube%fs_istart(3) + lz - 1
          ixx(3) = pad_feq(iz, db(3), .true.)

          call poisson_fft_gg_transform(ixx, temp, 0, mesh%sb%latt, coulb%qq, gg, modg2)

          !HH
          if(cube%fft%library.eq.FFTLIB_NFFT) then
             modg2=cube%Lfs(ix,1)**2+cube%Lfs(iy,2)**2+cube%Lfs(iz,3)**2
             r_c = default_r_c*M_TWO
          end if

          if(abs(modg2) > M_EPSILON) then
            select case(kernel)
            case(POISSON_FFT_KERNEL_SPH)
              fft_Coulb_FS(lx, ly, lz) = poisson_cutoff_3D_0D(sqrt(modg2),r_c)/modg2
            case(POISSON_FFT_KERNEL_CORRECTED)
              fft_Coulb_FS(lx, ly, lz) = M_ONE/modg2
            end select
          else
            select case(kernel)
            case(POISSON_FFT_KERNEL_SPH)
              fft_Coulb_FS(lx, ly, lz) = r_c**2/M_TWO
            case (POISSON_FFT_KERNEL_CORRECTED)
              fft_Coulb_FS(lx, ly, lz) = M_ZERO
            end select
          end if
        end do
      end do
    end do

    do iz = 1, cube%fs_n(3)
      do iy = 1, cube%fs_n(2)
        do ix = 1, cube%fs_n(1)
          fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
        end do
      end do
    end do

    call dfourier_space_op_init(coulb, cube, fft_coulb_fs, in_device = (kernel /= POISSON_FFT_KERNEL_CORRECTED))

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_3d_0d)
  end subroutine poisson_fft_build_3d_0d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> A. Castro et al., Phys. Rev. B 80, 033102 (2009)
  subroutine poisson_fft_build_2d_0d(namespace, mesh, cube, coulb)
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb

    type(spline_t) :: besselintf
    integer :: i, ix, iy, ixx(2), db(2), npoints
    FLOAT :: temp(2), vec, r_c, maxf, dk, default_r_c
    FLOAT, allocatable :: x(:), y(:)
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_2d_0d)

    db(1:2) = cube%rs_n_global(1:2)

    default_r_c = maxval(db(1:2)*mesh%spacing(1:2)/M_TWO)
    call get_cutoff(namespace, default_r_c, r_c)

    call spline_init(besselintf)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_Coulb_FS = M_ZERO
    temp(1:2) = M_TWO*M_PI/(db(1:2)*mesh%spacing(1:2))

    maxf = r_c * sqrt((temp(1)*db(1)/2)**2 + (temp(2)*db(2)/2)**2)
    dk = CNST(0.25) ! This seems to be reasonable.
    npoints = nint(maxf/dk)
    SAFE_ALLOCATE(x(1:npoints))
    SAFE_ALLOCATE(y(1:npoints))
    x(1) = M_ZERO
    y(1) = M_ZERO
    do i = 2, npoints
      x(i) = (i-1) * maxf / (npoints-1)
      y(i) = y(i-1) + poisson_cutoff_2D_0D(x(i-1), x(i))
    end do
    call spline_fit(npoints, x, y, besselintf)

    do iy = 1, cube%fs_n_global(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      do ix = 1, cube%fs_n_global(1)
        ixx(1) = pad_feq(ix, db(1), .true.)
        vec = sqrt( (temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
        if (vec > M_ZERO) then
           fft_coulb_fs(ix, iy, 1) = (M_TWO * M_PI / vec) * spline_eval(besselintf, vec*r_c)
        else
           fft_coulb_fs(ix, iy, 1) = M_TWO * M_PI * r_c
        end if
      end do
    end do

    call dfourier_space_op_init(coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)
    call spline_end(besselintf)
    POP_SUB(poisson_fft_build_2d_0d)
  end subroutine poisson_fft_build_2d_0d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> A. Castro et al., Phys. Rev. B 80, 033102 (2009)
  subroutine poisson_fft_build_2d_1d(namespace, mesh, cube, coulb)
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb

    integer :: ix, iy, ixx(2), db(2)
    FLOAT :: temp(2), r_c, gx, gy, default_r_c
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_2d_1d)

    db(1:2) = cube%rs_n_global(1:2)

    default_r_c = db(2)*mesh%spacing(2)/M_TWO
    call get_cutoff(namespace, default_r_c, r_c)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_Coulb_FS = M_ZERO
    temp(1:2) = M_TWO*M_PI/(db(1:2)*mesh%spacing(1:2))

    ! First, the term ix = 0 => gx = 0.
    fft_coulb_fs(1, 1, 1) = -M_FOUR * r_c * (log(r_c)-M_ONE)
    do iy = 2, cube%fs_n_global(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      gy = temp(2)*ixx(2)
      fft_coulb_fs(1, iy, 1) = -M_FOUR * poisson_cutoff_intcoslog(r_c, gy, M_ONE )
    end do

    do ix = 2, cube%fs_n_global(1)
      ixx(1) = pad_feq(ix, db(1), .true.)
      gx = temp(1)*ixx(1)
      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        gy = temp(2)*ixx(2)
        fft_coulb_fs(ix, iy, 1) = poisson_cutoff_2d_1d(gy, gx, r_c)
      end do
    end do

    call dfourier_space_op_init(coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)

    POP_SUB(poisson_fft_build_2d_1d)
  end subroutine poisson_fft_build_2d_1d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> A. Castro et al., Phys. Rev. B 80, 033102 (2009)
  subroutine poisson_fft_build_2d_2d(mesh, cube, coulb)
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb

    integer :: ix, iy, ixx(2), db(2)
    FLOAT :: temp(2), vec
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_2d_2d)

    db(1:2) = cube%rs_n_global(1:2)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_Coulb_FS = M_ZERO
    temp(1:2) = M_TWO*M_PI/(db(1:2)*mesh%spacing(1:2))

    do iy = 1, cube%fs_n_global(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      do ix = 1, cube%fs_n_global(1)
        ixx(1) = pad_feq(ix, db(1), .true.)
        vec = sqrt( (temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
        if (vec > M_ZERO) fft_coulb_fs(ix, iy, 1) = M_TWO * M_PI / vec
      end do
    end do

    call dfourier_space_op_init(coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_2d_2d)
  end subroutine poisson_fft_build_2d_2d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_1d_1d(mesh, cube, coulb, poisson_soft_coulomb_param)
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb
    FLOAT,                    intent(in)    :: poisson_soft_coulomb_param

    integer            :: ix, ixx
    FLOAT              :: g
    FLOAT, allocatable :: fft_coulb_fs(:, :, :)

    PUSH_SUB(poisson_fft_build_1d_1d)

    SAFE_ALLOCATE(fft_coulb_fs(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_coulb_fs = M_ZERO

    ! Fourier transform of Soft Coulomb interaction.
    do ix = 1, cube%fs_n_global(1)
      ixx = pad_feq(ix, cube%rs_n_global(1), .true.)
      g = (ixx + coulb%qq(1))*M_TWO*M_PI/abs(mesh%sb%latt%rlattice(1,1))
      if(abs(g) > CNST(1e-6)) then
        fft_coulb_fs(ix, 1, 1) = M_TWO * loct_bessel_k0(poisson_soft_coulomb_param*abs(g))
      end if
    end do

    call dfourier_space_op_init(coulb, cube, fft_coulb_fs)
    SAFE_DEALLOCATE_A(fft_coulb_fs)

    POP_SUB(poisson_fft_build_1d_1d)
  end subroutine poisson_fft_build_1d_1d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_1d_0d(namespace, mesh, cube, coulb, poisson_soft_coulomb_param)
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    type(cube_t),             intent(in)    :: cube
    type(fourier_space_op_t), intent(inout) :: coulb
    FLOAT,                    intent(in)    :: poisson_soft_coulomb_param

    integer            :: box(1), ixx(1), ix
    FLOAT              :: temp(1), g, r_c, default_r_c
    FLOAT, allocatable :: fft_coulb_fs(:, :, :)

    PUSH_SUB(poisson_fft_build_1d_0d)

    box(1:1) = cube%rs_n_global(1:1)

    default_r_c = box(1)*mesh%spacing(1)/M_TWO
    call get_cutoff(namespace, default_r_c, r_c)

    SAFE_ALLOCATE(fft_coulb_fs(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_coulb_fs = M_ZERO
    temp(1:1) = M_TWO*M_PI/(box(1:1)*mesh%spacing(1:1))

    ! Fourier transform of Soft Coulomb interaction.
    do ix = 1, cube%fs_n_global(1)
      ixx(1) = pad_feq(ix, box(1), .true.)
      g = temp(1)*ixx(1)
      fft_coulb_fs(ix, 1, 1) = poisson_cutoff_1D_0D(g, poisson_soft_coulomb_param, r_c)
    end do

    call dfourier_space_op_init(coulb, cube, fft_coulb_fs)
    SAFE_DEALLOCATE_A(fft_coulb_fs)

    POP_SUB(poisson_fft_build_1d_0d)
  end subroutine poisson_fft_build_1d_0d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_end(this)
    type(poisson_fft_t), intent(inout) :: this

    PUSH_SUB(poisson_fft_end)

    call fourier_space_op_end(this%coulb)

    POP_SUB(poisson_fft_end)
  end subroutine poisson_fft_end

#include "undef.F90"
#include "real.F90"
#include "poisson_fft_inc.F90"
#include "undef.F90"
#include "complex.F90"
#include "poisson_fft_inc.F90"


end module poisson_fft_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
