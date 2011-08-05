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
!! $Id: poisson.F90 2660 2007-01-23 15:11:54Z lorenzen $

#include "global.h"

module poisson_fft_m
  use cube_function_m
  use datasets_m
  use fft_m
  use fourier_space_m
  use geometry_m
  use global_m
  use lalg_basic_m
  use loct_m
  use loct_math_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use par_vec_m
  use parser_m
  use poisson_cutoff_m
  use profiling_m
  use simul_box_m
  use splines_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::                  &
    poisson_fft_build_1d_0d, &
    poisson_fft_build_1d_1d, &
    poisson_fft_build_2d_0d, &
    poisson_fft_build_2d_1d, &
    poisson_fft_build_2d_2d, &
    poisson_fft_build_3d_0d, &
    poisson_fft_build_3d_1d, &
    poisson_fft_build_3d_2d, &
    poisson_fft_build_3d_3d, &
    poisson_fft_end,         &
    poisson_fft

  integer, public, parameter :: &
       POISSON_FFT_SPH       =  0,      &
       POISSON_FFT_CYL       =  1,      &
       POISSON_FFT_PLA       =  2,      &
       POISSON_FFT_NOCUT     =  3,      &
       POISSON_FFT_CORRECTED =  4

  type(cube_function_t), public      :: fft_cf
  type(fourier_space_op_t) :: coulb

contains

  !-----------------------------------------------------------------
  subroutine poisson_fft_build_3d_3d(mesh)
    type(mesh_t), intent(inout) :: mesh
    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), idim
    FLOAT :: temp(MAX_DIM), modg2
    FLOAT :: gg(MAX_DIM)
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)
    
    PUSH_SUB(poisson_fft_build_3d_3d)
    
    ! double the box (or not) to perform the fourier transforms
    call mesh_double_box(mesh%sb, mesh, db)                 ! get dimensions of the double box
    ! allocate cube function where we will perform the ffts
    call cube_function_init(fft_cf, db)
    ! initialize the fft
    call dcube_function_fft_init(fft_cf, mesh%sb)
    ! dimensions may have been optimized
    db(1:3) = fft_cf%n(1:3)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
    fft_Coulb_FS = M_ZERO
    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))

    do ix = 1, fft_cf%nx
      ixx(1) = pad_feq(ix, db(1), .true.)
      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, db(3)
          ixx(3) = pad_feq(iz, db(3), .true.)

          gg(:) = temp(:)*ixx(:)
          gg(:) = matmul(gg, mesh%sb%klattice_primitive)
          do idim = 1, mesh%sb%dim
            gg(idim) = gg(idim) / lalg_nrm2(mesh%sb%dim, mesh%sb%klattice_primitive(:, idim))
          end do

          modg2 = sum(gg(:)**2)

          if(abs(modg2) > M_EPSILON) then
            fft_Coulb_FS(ix, iy, iz) = M_ONE/modg2
          else
            fft_Coulb_FS(ix, iy, iz) = M_ZERO
          end if
        end do
      end do

    end do

    forall(iz=1:fft_cf%n(3), iy=1:fft_cf%n(2), ix=1:fft_cf%nx)
      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
    end forall

    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_3d_3d)
  end subroutine poisson_fft_build_3d_3d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_3d_2d(mesh)
    type(mesh_t), intent(inout) :: mesh

    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), idim
    FLOAT :: temp(MAX_DIM), modg2
    FLOAT :: gpar, gx, gz, r_c, gg(MAX_DIM)
    FLOAT :: DELTA_R = CNST(1.0e-12)
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_2d)
    
    ! double the box to perform the fourier transforms
    call mesh_double_box(mesh%sb, mesh, db)                 ! get dimensions of the double box

    ! allocate cube function where we will perform the ffts
    call cube_function_init(fft_cf, db)

    ! initialize the fft
    call dcube_function_fft_init(fft_cf, mesh%sb)

    ! dimensions may have been optimized
    db(1:3) = fft_cf%n(1:3)

    call parse_float(datasets_check('PoissonCutoffRadius'),&
      maxval(db(:)*mesh%spacing(:)/M_TWO), r_c, units_inp%length)

    write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
      trim(units_abbrev(units_out%length)), '] = ',       &
      units_from_atomic(units_out%length, r_c)
    call messages_info(1)
    if ( r_c > maxval(db(:)*mesh%spacing(:)/M_TWO) + DELTA_R) then
      message(1) = 'Poisson cutoff radius is larger than cell size.'
      message(2) = 'You can see electrons in next cell(s).'
      call messages_warning(2)
    end if

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
    fft_Coulb_FS = M_ZERO

    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))

    do ix = 1, fft_cf%nx
      ixx(1) = pad_feq(ix, db(1), .true.)
      gx = temp(1)*ixx(1)
      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, db(3)
          ixx(3) = pad_feq(iz, db(3), .true.)

          gg(:) = temp(:)*ixx(:)

          gg(:) = matmul(gg, mesh%sb%klattice_primitive)
          do idim = 1, mesh%sb%dim
            gg(idim) = gg(idim) / lalg_nrm2(mesh%sb%dim, mesh%sb%klattice_primitive(:, idim))
          end do

          modg2 = sum(gg(:)**2)

          if(abs(modg2) > M_EPSILON) then
            gz = abs(gg(3))
            gpar = hypot(gg(1), gg(2))
            fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_2D(gpar,gz,r_c)/modg2
          else
            fft_Coulb_FS(ix, iy, iz) = -M_HALF*r_c**2
          end if
        end do
      end do

    end do

    forall(iz=1:fft_cf%n(3), iy=1:fft_cf%n(2), ix=1:fft_cf%nx)
      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
    end forall

    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_3d_2d)
  end subroutine poisson_fft_build_3d_2d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_3d_1d(mesh)
    type(mesh_t), intent(inout) :: mesh

    type(spline_t)     :: cylinder_cutoff_f
    FLOAT, allocatable :: x(:), y(:)
    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), k, ngp, idim
    FLOAT :: temp(MAX_DIM), modg2, xmax
    FLOAT :: gperp, gx, gy, gz, r_c, gg(MAX_DIM)
    FLOAT :: DELTA_R = CNST(1.0e-12)
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_1d)
    
    call mesh_double_box(mesh%sb, mesh, db)

    ! allocate cube function where we will perform the ffts
    call cube_function_init(fft_cf, db)

    ! initialize the fft
    call dcube_function_fft_init(fft_cf, mesh%sb)

    ! dimensions may have been optimized
    db(1:3) = fft_cf%n(1:3)

    call parse_float(datasets_check('PoissonCutoffRadius'),&
      maxval(db(:)*mesh%spacing(:)/M_TWO), r_c, units_inp%length)

    write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
      trim(units_abbrev(units_out%length)), '] = ',       &
      units_from_atomic(units_out%length, r_c)
    call messages_info(1)
    if ( r_c > maxval(db(:)*mesh%spacing(:)/M_TWO) + DELTA_R) then
      message(1) = 'Poisson cutoff radius is larger than cell size.'
      message(2) = 'You can see electrons in next cell(s).'
      call messages_warning(2)
    end if

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
    fft_Coulb_FS = M_ZERO

    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))

    if( mesh%sb%periodic_dim == 0 ) then
      ngp = 8*db(2)
      SAFE_ALLOCATE(x(1:ngp))
      SAFE_ALLOCATE(y(1:ngp))
    end if


    do ix = 1, fft_cf%nx
      ixx(1) = pad_feq(ix, db(1), .true.)
      gx = temp(1)*ixx(1)

      if( mesh%sb%periodic_dim == 0 ) then
        call spline_init(cylinder_cutoff_f)
        xmax = sqrt((temp(2)*db(2)/2)**2 + (temp(3)*db(3)/2)**2)
        do k = 1, ngp
          x(k) = (k-1)*(xmax/(ngp-1))
          y(k) = poisson_cutoff_3D_1D_finite(gx, x(k), M_TWO*mesh%sb%xsize, M_TWO*mesh%sb%rsize)
        end do
        call spline_fit(ngp, x, y, cylinder_cutoff_f)
      end if

      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, db(3)
          ixx(3) = pad_feq(iz, db(3), .true.)

          gg(:) = temp(:)*ixx(:)

          gg(:) = matmul(gg, mesh%sb%klattice_primitive)
          do idim = 1, mesh%sb%dim
            gg(idim) = gg(idim) / lalg_nrm2(mesh%sb%dim, mesh%sb%klattice_primitive(:, idim))
          end do

          modg2 = sum(gg(:)**2)

          if(abs(modg2) > M_EPSILON) then
            gperp = hypot(gg(2), gg(3))
            if (mesh%sb%periodic_dim==1) then
              fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_1D(abs(gx), gperp, r_c)/modg2
            else if (mesh%sb%periodic_dim==0) then
              gy = gg(2)
              gz = gg(3)
              if ((gz >= M_ZERO) .and. (gy >= M_ZERO)) then
                fft_Coulb_FS(ix, iy, iz) = spline_eval(cylinder_cutoff_f, gperp)
              end if
              if ((gz >= M_ZERO) .and. (gy < M_ZERO)) then
                fft_Coulb_FS(ix, iy, iz) = fft_Coulb_FS(ix, -ixx(2) + 1, iz)
              end if
              if ((gz < M_ZERO) .and. (gy >= M_ZERO)) then
                fft_Coulb_FS(ix, iy, iz) = fft_Coulb_FS(ix, iy, -ixx(3) + 1)
              end if
              if ((gz < M_ZERO) .and. (gy < M_ZERO) ) then
                fft_Coulb_FS(ix, iy, iz) = fft_Coulb_FS(ix, -ixx(2) + 1, -ixx(3) + 1)
              end if
            end if

          else
            if (mesh%sb%periodic_dim == 1) then
              fft_Coulb_FS(ix, iy, iz) = -(M_HALF*log(r_c) - M_FOURTH)*r_c**2
            else if (mesh%sb%periodic_dim == 0) then
              fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_1D_finite(M_ZERO, M_ZERO, &
                M_TWO*mesh%sb%xsize, M_TWO*mesh%sb%rsize)
            end if

          end if
        end do
      end do
 
      if( mesh%sb%periodic_dim == 0 ) call spline_end(cylinder_cutoff_f)
    end do

    forall(iz=1:fft_cf%n(3), iy=1:fft_cf%n(2), ix=1:fft_cf%nx)
      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
    end forall

    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)
    POP_SUB(poisson_fft_build_3d_1d)
  end subroutine poisson_fft_build_3d_1d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_3d_0d(mesh, poisson_solver)
    type(mesh_t), intent(inout) :: mesh
    integer,      intent(in)    :: poisson_solver

    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), idim
    FLOAT :: temp(MAX_DIM), modg2
    FLOAT :: gx, r_c, gg(MAX_DIM)
    FLOAT :: DELTA_R = CNST(1.0e-12)
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_0d)

    select case(poisson_solver)
    case(POISSON_FFT_SPH)
      call mesh_double_box(mesh%sb, mesh, db)
      db(:) = maxval(db)
    case(POISSON_FFT_CORRECTED)
      db(:) = mesh%idx%ll(:)
    end select

    ! allocate cube function where we will perform the ffts
    call cube_function_init(fft_cf, db)

    ! initialize the fft
    call dcube_function_fft_init(fft_cf, mesh%sb)

    ! dimensions may have been optimized
    db(1:3) = fft_cf%n(1:3)

    if (poisson_solver .ne. POISSON_FFT_CORRECTED) then
      call parse_float(datasets_check('PoissonCutoffRadius'),&
        maxval(db(:)*mesh%spacing(:)/M_TWO), r_c, units_inp%length)

      write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
        trim(units_abbrev(units_out%length)), '] = ',       &
        units_from_atomic(units_out%length, r_c)
      call messages_info(1)
      if ( r_c > maxval(db(:)*mesh%spacing(:)/M_TWO) + DELTA_R) then
        message(1) = 'Poisson cutoff radius is larger than cell size.'
        message(2) = 'You can see electrons in next cell(s).'
        call messages_warning(2)
      end if
    end if

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
    fft_Coulb_FS = M_ZERO

    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))


    do ix = 1, fft_cf%nx
      ixx(1) = pad_feq(ix, db(1), .true.)
      gx = temp(1)*ixx(1)
      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, db(3)
          ixx(3) = pad_feq(iz, db(3), .true.)

          gg(:) = temp(:)*ixx(:)
          gg(:) = matmul(gg, mesh%sb%klattice_primitive)
          do idim = 1, mesh%sb%dim
            gg(idim) = gg(idim) / lalg_nrm2(mesh%sb%dim, mesh%sb%klattice_primitive(:, idim))
          end do
          modg2 = sum(gg(:)**2)

          if(abs(modg2) > M_EPSILON) then
            select case(poisson_solver)
            case(POISSON_FFT_SPH)
              fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_0D(sqrt(modg2),r_c)/modg2
            case(POISSON_FFT_CORRECTED)
              fft_Coulb_FS(ix, iy, iz) = M_ONE/modg2
            end select
          else
            select case(poisson_solver)
            case(POISSON_FFT_SPH)
              fft_Coulb_FS(ix, iy, iz) = r_c**2/M_TWO
            case (POISSON_FFT_CORRECTED)
              fft_Coulb_FS(ix, iy, iz) = M_ZERO
            end select
          end if
        end do
      end do

    end do

    forall(iz=1:fft_cf%n(3), iy=1:fft_cf%n(2), ix=1:fft_cf%nx)
      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
    end forall

    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_3d_0d)
  end subroutine poisson_fft_build_3d_0d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_2d_0d(mesh)
    type(mesh_t), intent(in) :: mesh

    type(spline_t) :: besselintf
    integer :: i, ix, iy, ixx(MAX_DIM), db(MAX_DIM), npoints
    FLOAT :: temp(MAX_DIM), vec, r_c, maxf, dk
    FLOAT :: DELTA_R = CNST(1.0e-12)
    FLOAT, allocatable :: x(:), y(:)
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_2d_0d)

    ! double the box to perform the fourier transforms
    call mesh_double_box(mesh%sb, mesh, db)                 ! get dimensions of the double box
    db(1:2) = maxval(db)

    ! allocate cube function where we will perform the ffts
    call cube_function_init(fft_cf, db)      
    ! initialize fft
    call dcube_function_fft_init(fft_cf, mesh%sb)
    ! dimensions may have been optimized
    db(1:3) = fft_cf%n(1:3)

    call parse_float(datasets_check('PoissonCutoffRadius'),&
      maxval(db(1:2)*mesh%spacing(1:2)/M_TWO), r_c, units_inp%length)

    write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
      trim(units_abbrev(units_out%length)), '] = ',       &
      units_from_atomic(units_out%length, r_c)
    call messages_info(1)
    if ( r_c > maxval(db(1:2)*mesh%spacing(1:2)/M_TWO) + DELTA_R) then
      message(1) = 'Poisson cutoff radius is larger than cell size.'
      message(2) = 'You can see electrons in next cell(s).'
      call messages_warning(2)
    end if
    call spline_init(besselintf)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
    fft_Coulb_FS = M_ZERO
    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))

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

    do iy = 1, db(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      do ix = 1, fft_cf%nx
        ixx(1) = pad_feq(ix, db(1), .true.)
        vec = sqrt( (temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
        if (vec > M_ZERO) then
           fft_coulb_fs(ix, iy, 1) = (M_TWO * M_PI / vec) * spline_eval(besselintf, vec*r_c)
        else
           fft_coulb_fs(ix, iy, 1) = M_TWO * M_PI * r_c
        end if
      end do
    end do

    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)
    call spline_end(besselintf)
    POP_SUB(poisson_fft_build_2d_0d)
  end subroutine poisson_fft_build_2d_0d
  !-----------------------------------------------------------------
    

  !-----------------------------------------------------------------
  subroutine poisson_fft_build_2d_1d(mesh)
    type(mesh_t), intent(in) :: mesh

    integer :: ix, iy, ixx(MAX_DIM), db(MAX_DIM)
    FLOAT :: temp(MAX_DIM), vec, r_c, gx, gy
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_2d_1d)

    ! double the box to perform the fourier transforms
    call mesh_double_box(mesh%sb, mesh, db)                 ! get dimensions of the double box

    ! allocate cube function where we will perform the ffts
    call cube_function_init(fft_cf, db)      
    ! initialize fft
    call dcube_function_fft_init(fft_cf, mesh%sb)
    ! dimensions may have been optimized
    db(1:3) = fft_cf%n(1:3)

    r_c = M_TWO * mesh%sb%lsize(2)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
    fft_Coulb_FS = M_ZERO
    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))

    ! First, the term ix = 0 => gx = 0.
    fft_coulb_fs(1, 1, 1) = -M_FOUR * r_c * (log(r_c)-M_ONE)
    do iy = 2, db(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      gy = temp(2)*ixx(2)
      fft_coulb_fs(1, iy, 1) = -M_FOUR * poisson_cutoff_intcoslog(r_c, gy, M_ONE )
    end do

    do ix = 2, fft_cf%nx
      ixx(1) = pad_feq(ix, db(1), .true.)
      gx = temp(1)*ixx(1)
      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        gy = temp(2)*ixx(2)
        vec = sqrt( (temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
        fft_coulb_fs(ix, iy, 1) = poisson_cutoff_2d_1d(gy, gx, r_c)
      end do
    end do

    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)

    POP_SUB(poisson_fft_build_2d_1d)
  end subroutine poisson_fft_build_2d_1d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_2d_2d(mesh)
    type(mesh_t), intent(in) :: mesh

    integer :: ix, iy, ixx(MAX_DIM), db(MAX_DIM)
    FLOAT :: temp(MAX_DIM), vec
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_2d_2d)

    db(:) = mesh%idx%ll(:)
    ! allocate cube function where we will perform the ffts
    call cube_function_init(fft_cf, db)      
    ! initialize fft
    call dcube_function_fft_init(fft_cf, mesh%sb)
    ! dimensions may have been optimized
    db(1:3) = fft_cf%n(1:3)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
    fft_Coulb_FS = M_ZERO
    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))

    do iy = 1, db(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      do ix = 1, fft_cf%nx
        ixx(1) = pad_feq(ix, db(1), .true.)
        vec = sqrt( (temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
        if (vec > M_ZERO) fft_coulb_fs(ix, iy, 1) = M_TWO * M_PI / vec
      end do
    end do

    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_2d_2d)
  end subroutine poisson_fft_build_2d_2d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_1d_1d(mesh, poisson_soft_coulomb_param)
    type(mesh_t), intent(in) :: mesh
    FLOAT,        intent(in) :: poisson_soft_coulomb_param

    integer            :: box(MAX_DIM), ix
    FLOAT              :: g
    FLOAT, allocatable :: fft_coulb_fs(:, :, :)

    PUSH_SUB(poisson_fft_build_1d_1d)

    box = mesh%idx%ll
    call cube_function_init(fft_cf, box)
    call dcube_function_fft_init(fft_cf, mesh%sb)
    box(1:3) = fft_cf%n(1:3)

    SAFE_ALLOCATE(fft_coulb_fs(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
    fft_coulb_fs = M_ZERO

    ! Fourier transform of Soft Coulomb interaction.
    do ix = 1, fft_cf%nx
      g = ix*M_PI/mesh%sb%lsize(1) ! note that g is always positive with this definition
      fft_coulb_fs(ix, 1, 1) = M_TWO * loct_bessel_k0(poisson_soft_coulomb_param*g)
    end do

    call dfourier_space_op_init(coulb, fft_cf, fft_coulb_fs)
    SAFE_DEALLOCATE_A(fft_coulb_fs)
    
    POP_SUB(poisson_fft_build_1d_1d)
  end subroutine poisson_fft_build_1d_1d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_1d_0d(mesh, poisson_soft_coulomb_param)
    type(mesh_t), intent(in) :: mesh
    FLOAT,        intent(in) :: poisson_soft_coulomb_param

    integer            :: box(MAX_DIM), ixx(MAX_DIM), ix
    FLOAT              :: temp(MAX_DIM), g, r_c
    FLOAT, allocatable :: fft_coulb_fs(:, :, :)

    PUSH_SUB(poisson_fft_build_1d_0d)

    call mesh_double_box(mesh%sb, mesh, box)
    
    call cube_function_init(fft_cf, box)
    call dcube_function_fft_init(fft_cf, mesh%sb)
    box(1:3) = fft_cf%n(1:3)

    r_c = box(1)*mesh%spacing(1)/M_TWO

    SAFE_ALLOCATE(fft_coulb_fs(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
    fft_coulb_fs = M_ZERO
    temp(:) = M_TWO*M_PI/(box(:)*mesh%spacing(:))

    ! Fourier transform of Soft Coulomb interaction.
    do ix = 1, fft_cf%nx
      ixx(1) = pad_feq(ix, box(1), .true.)
      g = temp(1)*ixx(1)
      fft_coulb_fs(ix, 1, 1) = poisson_cutoff_1D_0D(g, poisson_soft_coulomb_param, r_c)
    end do

    call dfourier_space_op_init(coulb, fft_cf, fft_coulb_fs)
    SAFE_DEALLOCATE_A(fft_coulb_fs)
    
    POP_SUB(poisson_fft_build_1d_0d)
  end subroutine poisson_fft_build_1d_0d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_end()
    PUSH_SUB(poisson_fft.end)

    call fourier_space_op_end(coulb)
    call cube_function_fft_end(fft_cf)
    call cube_function_end(fft_cf)

    POP_SUB(poisson_fft.end)
  end subroutine poisson_fft_end

  !-----------------------------------------------------------------

  subroutine poisson_fft(mesh, pot, rho, average_to_zero)
    type(mesh_t),      intent(in)  :: mesh
    FLOAT,             intent(out) :: pot(:)
    FLOAT,             intent(in)  :: rho(:)
    logical, optional, intent(in)  :: average_to_zero

    FLOAT, allocatable :: rho_global(:), pot_global(:)

    FLOAT :: average
    type(profile_t), save :: prof_bcast, prof_sct
    integer :: default_fft_library, fft_library, ii, last, first
    
    PUSH_SUB(poisson_fft)
    
    average = M_ZERO !this avoids a non-initialized warning
    
    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(rho_global(1:mesh%np_global))
      SAFE_ALLOCATE(pot_global(1:mesh%np_global))
    end if
    !Safe memory if PFFT is used  
    if (fft_cf%fft_library /= PFFT_LIB) then
      call dcube_function_alloc_RS(fft_cf)          ! allocate the cube in real space
    end if
    
    ! put the density in the cube
    if(mesh%parallel_in_domains) then
      if (fft_cf%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
        !all the data has to be distributed for the PFFT library
        call dvec_allgather(mesh%vp, rho_global, rho)
        call dmesh_to_cube_parallel(mesh, rho_global, fft_cf)
#else
        write(message(1),'(a)')'You have selected the PFFT for FFT, but it is not compiled.'
        call messages_fatal(1)
#endif
      else 
#if defined HAVE_MPI
        call dvec_gather(mesh%vp, mesh%vp%root, rho_global, rho)
        call dmesh_to_cube(mesh, rho_global, fft_cf)   
      end if
#endif
    else !not parallel in domains
      if (fft_cf%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
        ! Serial execution with PFFT
        call dmesh_to_cube_parallel(mesh, rho, fft_cf)
#endif 
      else
        call dmesh_to_cube(mesh, rho, fft_cf)
      end if
    end if

    ! apply the Couloumb term in Fourier space
    call dfourier_space_op_apply(coulb, fft_cf)

    !now the cube has the potential
    if(present(average_to_zero)) then
      if(average_to_zero) average = cube_function_surface_average(fft_cf)
#if defined HAVE_MPI
      ! Only root has the right average.
      if(mesh%parallel_in_domains) call MPI_Bcast(average, 1, MPI_FLOAT, 0, mesh%mpi_grp%comm, mpi_err)
#endif
    end if

    ! move the potential back to the mesh
    if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)   
      if (fft_cf%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
        call dcube_to_mesh_parallel(mesh, fft_cf, pot_global)
        call profiling_in(prof_sct, 'SCT')
        !skip the scatter, because all the nodes have each part of pot_global
        first = mesh%vp%xlocal(mesh%vp%partno)
        last = mesh%vp%np_local(mesh%vp%partno)   
        do ii = 1, last
          pot(ii) = pot_global(mesh%vp%local(ii+first-1))
        end do
        call profiling_out(prof_sct)
#else
        write(message(1),'(a)')'You have selected the PFFT for FFT, but it is not compiled.'
        call messages_fatal(1)
#endif
      else
        call dcube_to_mesh(mesh, fft_cf, pot_global)
        call dvec_scatter(mesh%vp, mesh%vp%root, pot_global, pot)
      end if
#endif
    else !not parallel in domains
      if (fft_cf%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
        ! Serial execution with PFFT
        call dcube_to_mesh_parallel(mesh, fft_cf, pot)!_global)
        !pot = pot_global
#endif 
      else
        call dcube_to_mesh(mesh, fft_cf, pot)
      end if
    end if

    if(present(average_to_zero)) then
      if(average_to_zero) pot(1:mesh%np) = pot(1:mesh%np) - average
    end if
    
    if (fft_cf%fft_library /= PFFT_LIB) then
      call dcube_function_free_RS(fft_cf)           ! memory is no longer needed
    end if
    
    if(mesh%parallel_in_domains) then
      SAFE_DEALLOCATE_A(rho_global)
      SAFE_DEALLOCATE_A(pot_global)
    end if
    
    POP_SUB(poisson_fft)
  end subroutine poisson_fft

end module poisson_fft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
