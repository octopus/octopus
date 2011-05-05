!! Copyright (C) 2011 J. Alberdi, P. Garcia Risueño
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
!! $Id: poisson_pfft.F90 2660 2011-02-08 15:11:54Z pablogr $

#include "global.h"

module poisson_pfft_m
  use cube_function_m
  use datasets_m
  use fft_m
  use pfft_m
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
    poisson_pfft_build_3d_3d, &
    poisson_pfft_end,         &
    poisson_pfft

  integer, public, parameter :: &
       POISSON_PFFT_SPH       =  0,      &
       POISSON_PFFT_CYL       =  1,      &
       POISSON_PFFT_PLA       =  2,      &
       POISSON_PFFT_NOCUT     =  3,      &
       POISSON_PFFT_CORRECTED =  4

  type(cube_function_t), public      :: fft_cf
  type(fourier_space_op_t) :: coulb

contains

  !-----------------------------------------------------------------
  subroutine poisson_pfft_build_3d_3d(mesh)
    type(mesh_t), intent(inout) :: mesh

#ifdef HAVE_PFFT
    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), idim
    FLOAT :: temp(MAX_DIM), modg2
    FLOAT :: gg(MAX_DIM)
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_pfft_build_3d_3d)
    
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
            gg(idim) = gg(idim) / lalg_nrm2(MAX_DIM, mesh%sb%klattice_primitive(:, idim))
          end do

          modg2 = sum(gg(:)**2)

          if(modg2 /= M_ZERO) then
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
    POP_SUB(poisson_pfft_build_3d_3d)
#endif
  end subroutine poisson_pfft_build_3d_3d
  !-----------------------------------------------------------------

  ! \todo I`m not going to implement then, for the moment. Joseba
  !-----------------------------------------------------------------
#if 0
!!$  subroutine poisson_pfft_build_3d_2d(mesh)
!!$    type(mesh_t), intent(inout) :: mesh
!!$
!!$    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), idim
!!$    FLOAT :: temp(MAX_DIM), modg2
!!$    FLOAT :: gpar, gx, gz, r_c, gg(MAX_DIM)
!!$    FLOAT :: DELTA_R = CNST(1.0e-12)
!!$    FLOAT, allocatable :: fft_coulb_FS(:,:,:)
!!$
!!$    PUSH_SUB(poisson_pfft_build_3d_2d)
!!$    
!!$    ! double the box to perform the fourier transforms
!!$    call mesh_double_box(mesh%sb, mesh, db)                 ! get dimensions of the double box
!!$
!!$    ! allocate cube function where we will perform the ffts
!!$    call dcf_new(db, fft_cf)
!!$
!!$    ! initialize the fft
!!$    call dcf_fft_init(fft_cf, mesh%sb)
!!$
!!$    ! dimensions may have been optimized
!!$    db(1:3) = fft_cf%n(1:3)
!!$
!!$    call parse_float(datasets_check('PoissonCutoffRadius'),&
!!$      maxval(db(:)*mesh%spacing(:)/M_TWO), r_c, units_inp%length)
!!$
!!$    write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
!!$      trim(units_abbrev(units_out%length)), '] = ',       &
!!$      units_from_atomic(units_out%length, r_c)
!!$    call messages_info(1)
!!$    if ( r_c > maxval(db(:)*mesh%spacing(:)/M_TWO) + DELTA_R) then
!!$      message(1) = 'Poisson cutoff radius is larger than cell size.'
!!$      message(2) = 'You can see electrons in next cell(s).'
!!$      call messages_warning(2)
!!$    end if
!!$
!!$    ! store the fourier transform of the Coulomb interaction
!!$    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
!!$    fft_Coulb_FS = M_ZERO
!!$
!!$    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))
!!$
!!$    do ix = 1, fft_cf%nx
!!$      ixx(1) = pad_feq(ix, db(1), .true.)
!!$      gx = temp(1)*ixx(1)
!!$      do iy = 1, db(2)
!!$        ixx(2) = pad_feq(iy, db(2), .true.)
!!$        do iz = 1, db(3)
!!$          ixx(3) = pad_feq(iz, db(3), .true.)
!!$
!!$          gg(:) = temp(:)*ixx(:)
!!$
!!$          gg(:) = matmul(gg, mesh%sb%klattice_primitive)
!!$          do idim = 1, mesh%sb%dim
!!$            gg(idim) = gg(idim) / lalg_nrm2(MAX_DIM, mesh%sb%klattice_primitive(:, idim))
!!$          end do
!!$
!!$          modg2 = sum(gg(:)**2)
!!$
!!$          if(modg2 /= M_ZERO) then
!!$            gz = abs(gg(3))
!!$            gpar = hypot(gg(1), gg(2))
!!$            fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_2D(gpar,gz,r_c)/modg2
!!$          else
!!$            fft_Coulb_FS(ix, iy, iz) = -M_HALF*r_c**2
!!$          end if
!!$        end do
!!$      end do
!!$
!!$    end do
!!$
!!$    forall(iz=1:fft_cf%n(3), iy=1:fft_cf%n(2), ix=1:fft_cf%nx)
!!$      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
!!$    end forall
!!$
!!$    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)
!!$
!!$    SAFE_DEALLOCATE_A(fft_Coulb_FS)
!!$    POP_SUB(poisson_pfft_build_3d_2d)
!!$  end subroutine poisson_pfft_build_3d_2d
!!$  !-----------------------------------------------------------------
!!$
!!$
!!$  !-----------------------------------------------------------------
!!$  subroutine poisson_pfft_build_3d_1d(mesh)
!!$    type(mesh_t), intent(inout) :: mesh
!!$
!!$    type(spline_t)     :: cylinder_cutoff_f
!!$    FLOAT, allocatable :: x(:), y(:)
!!$    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), k, ngp, idim
!!$    FLOAT :: temp(MAX_DIM), modg2, xmax
!!$    FLOAT :: gperp, gx, gy, gz, r_c, gg(MAX_DIM)
!!$    FLOAT :: DELTA_R = CNST(1.0e-12)
!!$    FLOAT, allocatable :: fft_coulb_FS(:,:,:)
!!$
!!$    PUSH_SUB(poisson_pfft_build_3d_1d)
!!$    
!!$    call mesh_double_box(mesh%sb, mesh, db)
!!$
!!$    ! allocate cube function where we will perform the ffts
!!$    call dcf_new(db, fft_cf)
!!$
!!$    ! initialize the fft
!!$    call dcf_fft_init(fft_cf, mesh%sb)
!!$
!!$    ! dimensions may have been optimized
!!$    db(1:3) = fft_cf%n(1:3)
!!$
!!$    call parse_float(datasets_check('PoissonCutoffRadius'),&
!!$      maxval(db(:)*mesh%spacing(:)/M_TWO), r_c, units_inp%length)
!!$
!!$    write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
!!$      trim(units_abbrev(units_out%length)), '] = ',       &
!!$      units_from_atomic(units_out%length, r_c)
!!$    call messages_info(1)
!!$    if ( r_c > maxval(db(:)*mesh%spacing(:)/M_TWO) + DELTA_R) then
!!$      message(1) = 'Poisson cutoff radius is larger than cell size.'
!!$      message(2) = 'You can see electrons in next cell(s).'
!!$      call messages_warning(2)
!!$    end if
!!$
!!$    ! store the fourier transform of the Coulomb interaction
!!$    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
!!$    fft_Coulb_FS = M_ZERO
!!$
!!$    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))
!!$
!!$    if( mesh%sb%periodic_dim == 0 ) then
!!$      ngp = 8*db(2)
!!$      SAFE_ALLOCATE(x(1:ngp))
!!$      SAFE_ALLOCATE(y(1:ngp))
!!$    end if
!!$
!!$
!!$    do ix = 1, fft_cf%nx
!!$      ixx(1) = pad_feq(ix, db(1), .true.)
!!$      gx = temp(1)*ixx(1)
!!$
!!$      if( mesh%sb%periodic_dim == 0 ) then
!!$        call spline_init(cylinder_cutoff_f)
!!$        xmax = sqrt((temp(2)*db(2)/2)**2 + (temp(3)*db(3)/2)**2)
!!$        do k = 1, ngp
!!$          x(k) = (k-1)*(xmax/(ngp-1))
!!$          y(k) = poisson_cutoff_3D_1D_finite(gx, x(k), M_TWO*mesh%sb%xsize, M_TWO*mesh%sb%rsize)
!!$        end do
!!$        call spline_fit(ngp, x, y, cylinder_cutoff_f)
!!$      end if
!!$
!!$      do iy = 1, db(2)
!!$        ixx(2) = pad_feq(iy, db(2), .true.)
!!$        do iz = 1, db(3)
!!$          ixx(3) = pad_feq(iz, db(3), .true.)
!!$
!!$          gg(:) = temp(:)*ixx(:)
!!$
!!$          gg(:) = matmul(gg, mesh%sb%klattice_primitive)
!!$          do idim = 1, mesh%sb%dim
!!$            gg(idim) = gg(idim) / lalg_nrm2(MAX_DIM, mesh%sb%klattice_primitive(:, idim))
!!$          end do
!!$
!!$          modg2 = sum(gg(:)**2)
!!$
!!$          if(modg2 /= M_ZERO) then
!!$            gperp = hypot(gg(2), gg(3))
!!$            if (mesh%sb%periodic_dim==1) then
!!$              fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_1D(abs(gx), gperp, r_c)/modg2
!!$            else if (mesh%sb%periodic_dim==0) then
!!$              gy = gg(2)
!!$              gz = gg(3)
!!$              if ((gz >= M_ZERO) .and. (gy >= M_ZERO)) then
!!$                fft_Coulb_FS(ix, iy, iz) = spline_eval(cylinder_cutoff_f, gperp)
!!$              end if
!!$              if ((gz >= M_ZERO) .and. (gy < M_ZERO)) then
!!$                fft_Coulb_FS(ix, iy, iz) = fft_Coulb_FS(ix, -ixx(2) + 1, iz)
!!$              end if
!!$              if ((gz < M_ZERO) .and. (gy >= M_ZERO)) then
!!$                fft_Coulb_FS(ix, iy, iz) = fft_Coulb_FS(ix, iy, -ixx(3) + 1)
!!$              end if
!!$              if ((gz < M_ZERO) .and. (gy < M_ZERO) ) then
!!$                fft_Coulb_FS(ix, iy, iz) = fft_Coulb_FS(ix, -ixx(2) + 1, -ixx(3) + 1)
!!$              end if
!!$            end if
!!$
!!$          else
!!$            if (mesh%sb%periodic_dim == 1) then
!!$              fft_Coulb_FS(ix, iy, iz) = -(M_HALF*log(r_c) - M_FOURTH)*r_c**2
!!$            else if (mesh%sb%periodic_dim == 0) then
!!$              fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_1D_finite(M_ZERO, M_ZERO, &
!!$                M_TWO*mesh%sb%xsize, M_TWO*mesh%sb%rsize)
!!$            end if
!!$
!!$          end if
!!$        end do
!!$      end do
!!$ 
!!$      if( mesh%sb%periodic_dim == 0 ) call spline_end(cylinder_cutoff_f)
!!$    end do
!!$
!!$    forall(iz=1:fft_cf%n(3), iy=1:fft_cf%n(2), ix=1:fft_cf%nx)
!!$      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
!!$    end forall
!!$
!!$    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)
!!$
!!$    SAFE_DEALLOCATE_A(fft_Coulb_FS)
!!$    SAFE_DEALLOCATE_A(x)
!!$    SAFE_DEALLOCATE_A(y)
!!$    POP_SUB(poisson_pfft_build_3d_1d)
!!$  end subroutine poisson_pfft_build_3d_1d
!!$  !-----------------------------------------------------------------
!!$
!!$
  !-----------------------------------------------------------------
  subroutine poisson_pfft_build_3d_0d(mesh, poisson_solver)
    type(mesh_t), intent(inout) :: mesh
    integer,      intent(in)    :: poisson_solver

    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), idim
    FLOAT :: temp(MAX_DIM), modg2
    FLOAT :: gx, r_c, gg(MAX_DIM)
    FLOAT :: DELTA_R = CNST(1.0e-12)
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_pfft_build_3d_0d)
    
    select case(poisson_solver)
    case(POISSON_PFFT_SPH)
      call mesh_double_box(mesh%sb, mesh, db)
      db(:) = maxval(db)
    case(POISSON_PFFT_CORRECTED)
      db(:) = mesh%idx%ll(:)
    end select

    ! allocate cube function where we will perform the ffts
    call dcf_new(db, fft_cf)

    ! initialize the fft
    call dcf_fft_init(fft_cf, mesh%sb)

    ! dimensions may have been optimized
    db(1:3) = fft_cf%n(1:3)

    if (poisson_solver .ne. POISSON_PFFT_CORRECTED) then
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
            gg(idim) = gg(idim) / lalg_nrm2(MAX_DIM, mesh%sb%klattice_primitive(:, idim))
          end do
          modg2 = sum(gg(:)**2)

          if(modg2 /= M_ZERO) then
            select case(poisson_solver)
            case(POISSON_PFFT_SPH)
              fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_0D(sqrt(modg2),r_c)/modg2
            case(POISSON_PFFT_CORRECTED)
              fft_Coulb_FS(ix, iy, iz) = M_ONE/modg2
            end select
          else
            select case(poisson_solver)
            case(POISSON_PFFT_SPH)
              fft_Coulb_FS(ix, iy, iz) = r_c**2/M_TWO
            case (POISSON_PFFT_CORRECTED)
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
    POP_SUB(poisson_pfft_build_3d_0d)
  end subroutine poisson_pfft_build_3d_0d
!!$  !-----------------------------------------------------------------
!!$  subroutine poisson_pfft_build_3d_0d(mesh, poisson_solver)
!!$    type(mesh_t), intent(inout) :: mesh
!!$    integer,      intent(in)    :: poisson_solver
!!$
!!$    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), idim
!!$    FLOAT :: temp(MAX_DIM), modg2
!!$    FLOAT :: gx, r_c, gg(MAX_DIM)
!!$    FLOAT :: DELTA_R = CNST(1.0e-12)
!!$    FLOAT, allocatable :: fft_coulb_FS(:,:,:)
!!$
!!$    PUSH_SUB(poisson_pfft_build_3d_0d)
!!$    
!!$    select case(poisson_solver)
!!$    case(POISSON_PFFT_SPH)
!!$      call mesh_double_box(mesh%sb, mesh, db)
!!$      db(:) = maxval(db)
!!$    case(POISSON_PFFT_CORRECTED)
!!$      db(:) = mesh%idx%ll(:)
!!$    end select
!!$
!!$    ! allocate cube function where we will perform the ffts
!!$    call dcf_new(db, fft_cf)
!!$
!!$    ! initialize the fft
!!$    call dcf_fft_init(fft_cf, mesh%sb)
!!$
!!$    ! dimensions may have been optimized
!!$    db(1:3) = fft_cf%n(1:3)
!!$
!!$    if (poisson_solver .ne. POISSON_PFFT_CORRECTED) then
!!$      call parse_float(datasets_check('PoissonCutoffRadius'),&
!!$        maxval(db(:)*mesh%spacing(:)/M_TWO), r_c, units_inp%length)
!!$
!!$      write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
!!$        trim(units_abbrev(units_out%length)), '] = ',       &
!!$        units_from_atomic(units_out%length, r_c)
!!$      call messages_info(1)
!!$      if ( r_c > maxval(db(:)*mesh%spacing(:)/M_TWO) + DELTA_R) then
!!$        message(1) = 'Poisson cutoff radius is larger than cell size.'
!!$        message(2) = 'You can see electrons in next cell(s).'
!!$        call messages_warning(2)
!!$      end if
!!$    end if
!!$
!!$    ! store the fourier transform of the Coulomb interaction
!!$    SAFE_ALLOCATE(fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), 1:fft_cf%n(3)))
!!$    fft_Coulb_FS = M_ZERO
!!$
!!$    temp(:) = M_TWO*M_PI/(db(:)*mesh%spacing(:))
!!$
!!$
!!$    do ix = 1, fft_cf%nx
!!$      ixx(1) = pad_feq(ix, db(1), .true.)
!!$      gx = temp(1)*ixx(1)
!!$      do iy = 1, db(2)
!!$        ixx(2) = pad_feq(iy, db(2), .true.)
!!$        do iz = 1, db(3)
!!$          ixx(3) = pad_feq(iz, db(3), .true.)
!!$
!!$          gg(:) = temp(:)*ixx(:)
!!$          gg(:) = matmul(gg, mesh%sb%klattice_primitive)
!!$          do idim = 1, mesh%sb%dim
!!$            gg(idim) = gg(idim) / lalg_nrm2(MAX_DIM, mesh%sb%klattice_primitive(:, idim))
!!$          end do
!!$          modg2 = sum(gg(:)**2)
!!$
!!$          if(modg2 /= M_ZERO) then
!!$            select case(poisson_solver)
!!$            case(POISSON_PFFT_SPH)
!!$              fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_0D(sqrt(modg2),r_c)/modg2
!!$            case(POISSON_PFFT_CORRECTED)
!!$              fft_Coulb_FS(ix, iy, iz) = M_ONE/modg2
!!$            end select
!!$          else
!!$            select case(poisson_solver)
!!$            case(POISSON_PFFT_SPH)
!!$              fft_Coulb_FS(ix, iy, iz) = r_c**2/M_TWO
!!$            case (POISSON_PFFT_CORRECTED)
!!$              fft_Coulb_FS(ix, iy, iz) = M_ZERO
!!$            end select
!!$          end if
!!$        end do
!!$      end do
!!$
!!$    end do
!!$
!!$    forall(iz=1:fft_cf%n(3), iy=1:fft_cf%n(2), ix=1:fft_cf%nx)
!!$      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
!!$    end forall
!!$
!!$    call dfourier_space_op_init(coulb, fft_cf, fft_Coulb_FS)
!!$
!!$    SAFE_DEALLOCATE_A(fft_Coulb_FS)
!!$    POP_SUB(poisson_pfft_build_3d_0d)
!!$  end subroutine poisson_pfft_build_3d_0d
!!$  !-----------------------------------------------------------------
#endif
  !-----------------------------------------------------------------
  subroutine poisson_pfft_end()
    PUSH_SUB(poisson_pfft.end)

    call dcube_function_free(fft_cf)
    call dfourier_space_op_end(coulb)

    POP_SUB(poisson_pfft.end)
  end subroutine poisson_pfft_end

  !-----------------------------------------------------------------

  subroutine poisson_pfft(mesh, pot, rho, average_to_zero)
    type(mesh_t),      intent(in)  :: mesh
    FLOAT,             intent(out) :: pot(:)
    FLOAT,             intent(in)  :: rho(:)
    logical, optional, intent(in)  :: average_to_zero

    FLOAT, allocatable :: rho_global(:), pot_global(:)

    FLOAT :: average

    PUSH_SUB(poisson_pfft)

    average=M_ZERO !this avoids a non-initialized warning
    
    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(rho_global(1:mesh%np_global))
      SAFE_ALLOCATE(pot_global(1:mesh%np_global))
    end if

    call dcube_function_alloc_RS(fft_cf)          ! allocate the cube in real space

    ! put the density in the cube
    if(mesh%parallel_in_domains) then
#if defined HAVE_MPI
      call dvec_gather(mesh%vp, mesh%vp%root, rho_global, rho)
      call dmesh_to_cube(mesh, rho_global, fft_cf) 
#endif
    else
      call dmesh_to_cube(mesh, rho, fft_cf)
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
      call dcube_to_mesh(mesh, fft_cf, pot_global)
      call dvec_scatter(mesh%vp, mesh%vp%root, pot_global, pot)
#endif
    else
      call dcube_to_mesh(mesh, fft_cf, pot)
    end if

    if(present(average_to_zero)) then
      if(average_to_zero) pot(1:mesh%np) = pot(1:mesh%np) - average
    end if

    call dcube_function_free_RS(fft_cf)           ! memory is no longer needed

    if(mesh%parallel_in_domains) then
      SAFE_DEALLOCATE_A(rho_global)
      SAFE_DEALLOCATE_A(pot_global)
    end if

    POP_SUB(poisson_pfft)
  end subroutine poisson_pfft

end module poisson_pfft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
