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
  use datasets_m
  use geometry_m
  use global_m
  use lib_oct_m
  use lib_oct_parser_m
  use mesh_m
  use messages_m
  use mpi_m
  use profiling_m
  use units_m
#ifdef HAVE_FFT
  use fft_m
  use cube_function_m
#endif
  use functions_m
  use grid_m
  use mesh_function_m
  use par_vec_m
  use poisson_cutoffs_m
  use lib_oct_gsl_spline_m

  implicit none

  private
  public :: & 
       poisson_fft_build_2d, &
       poisson_fft_build_3d, &
       poisson_fft_end,  &
       poisson_fft
  

  integer, public, parameter :: &
       FFT_SPH       =  0, &
       FFT_CYL       =  1, &
       FFT_PLA       =  2, &
       FFT_NOCUT     =  3, &
       FFT_CORRECTED =  4

#ifdef HAVE_FFT  
  type(dcf_t), public :: fft_cf
  FLOAT, pointer :: fft_coulb_FS(:,:,:)
#endif

contains

  !-----------------------------------------------------------------
  subroutine poisson_fft_build_3d(gr, poisson_solver)
    type(grid_t), intent(inout) :: gr
    integer, intent(in) :: poisson_solver
#if defined(HAVE_FFT)
    type(loct_spline_t) :: cylinder_cutoff_f
    FLOAT, allocatable :: x(:), y(:)
    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), k, ngp
    FLOAT :: temp(MAX_DIM), modg2, xmax
    FLOAT :: gpar, gperp, gx, gy, gz, r_c
    FLOAT :: DELTA_R = CNST(1.0e-12)

    ! double the box to perform the fourier transforms
    if(poisson_solver.ne.FFT_CORRECTED) then
       call mesh_double_box(gr%sb, gr%m, db)                 ! get dimensions of the double box
       if (poisson_solver == FFT_SPH) db(:) = maxval(db)
    else
       db(:) = gr%m%l(:)
    end if

    call dcf_new(db, fft_cf)    ! allocate cube function where we will perform
    call dcf_fft_init(fft_cf, gr%sb)   ! the ffts
    db = fft_cf%n               ! dimensions may have been optimized

    if (poisson_solver <= FFT_PLA .and. poisson_solver .ne. FFT_CORRECTED) then
       call loct_parse_float(check_inp('PoissonCutoffRadius'),&
            maxval(db(:)*gr%m%h(:)/M_TWO)/units_inp%length%factor , r_c)
       r_c = r_c*units_inp%length%factor
       write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
            trim(units_out%length%abbrev), '] = ',       &
            r_c/units_out%length%factor
       call write_info(1)
       if ( r_c > maxval(db(:)*gr%m%h(:)/M_TWO) + DELTA_R) then
          message(1) = 'Poisson cutoff radius is larger than cell size.'
          message(2) = 'You can see electrons in next cell(s).'
          call write_warning(2)
       end if
    end if

    ! store the fourier transform of the Coulomb interaction
    ALLOCATE(fft_Coulb_FS(fft_cf%nx, fft_cf%n(2), fft_cf%n(3)), fft_cf%nx*fft_cf%n(2)*fft_cf%n(3))
    fft_Coulb_FS = M_ZERO

    temp(:) = M_TWO*M_PI/(db(:)*gr%m%h(:))

    if( (poisson_solver .eq. FFT_CYL)  .and. (gr%sb%periodic_dim == 0) ) then
      ngp = 8*db(2)
      ALLOCATE(x(ngp), ngp)
      ALLOCATE(y(ngp), ngp)
    end if


    do ix = 1, fft_cf%nx
      ixx(1) = pad_feq(ix, db(1), .true.)
      gx = temp(1)*ixx(1)

      if( (poisson_solver .eq. FFT_CYL)  .and. (gr%sb%periodic_dim == 0) ) then
        call loct_spline_init(cylinder_cutoff_f)
        xmax = sqrt((temp(2)*db(2)/2)**2 + (temp(3)*db(3)/2)**2)
        do k = 1, ngp
          x(k) = (k-1)*(xmax/(ngp-1))
          y(k) = poisson_cutoff_finite_cylinder(gx, x(k), M_TWO*gr%m%sb%xsize, M_TWO*gr%m%sb%rsize)
        end do
        call loct_spline_fit(ngp, x, y, cylinder_cutoff_f)
      end if

      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, db(3)
          ixx(3) = pad_feq(iz, db(3), .true.)

             modg2 = sum((temp(:)*ixx(:))**2)
             if(modg2 /= M_ZERO) then
                select case(poisson_solver)
                case(FFT_SPH)
                   fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_sphere(sqrt(modg2),r_c)/modg2

                case(FFT_CYL)
                   gperp = sqrt((temp(2)*ixx(2))**2+(temp(3)*ixx(3))**2)
                   if (gr%sb%periodic_dim==1) then
                     fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_infinite_cylinder(abs(gx), gperp, r_c)/modg2

                   else if (gr%sb%periodic_dim==0) then
                     gy = temp(2)*ixx(2)
                     gz = temp(3)*ixx(3)
                     if ((gz >= M_ZERO) .and. (gy >= M_ZERO)) then
                       fft_Coulb_FS(ix, iy, iz) = loct_splint(cylinder_cutoff_f, gperp)
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

                 case(FFT_PLA)
                   gz = abs(temp(3)*ixx(3))
                   gpar = sqrt((temp(1)*ixx(1))**2+(temp(2)*ixx(2))**2)
                   fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_slab(gpar,gz,r_c)/modg2

                case(FFT_NOCUT, FFT_CORRECTED)
                   fft_Coulb_FS(ix, iy, iz) = M_ONE/modg2
                end select

             else
                select case(poisson_solver)
                case(FFT_SPH)
                  fft_Coulb_FS(ix, iy, iz) = r_c**2/M_TWO

                case (FFT_CYL)
                  if (gr%sb%periodic_dim == 1) then
                    fft_Coulb_FS(ix, iy, iz) = -(M_HALF*log(r_c) - M_FOURTH)*r_c**2

                  else if (gr%sb%periodic_dim == 0) then
                    fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_finite_cylinder(M_ZERO, M_ZERO, &
                         M_TWO*gr%m%sb%xsize, M_TWO*gr%m%sb%rsize)
                  end if

                case(FFT_PLA)
                  fft_Coulb_FS(ix, iy, iz) = -M_HALF*r_c**2

                case (FFT_NOCUT, FFT_CORRECTED)
                  fft_Coulb_FS(ix, iy, iz) = M_ZERO

                end select
             end if
          end do
       end do

      if( (poisson_solver .eq. FFT_CYL) .and. (gr%sb%periodic_dim == 0) ) then
        call loct_spline_end(cylinder_cutoff_f)
      end if
    end do

    do k=1, fft_cf%n(3)
      fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), k) = M_FOUR*M_PI*fft_Coulb_FS(1:fft_cf%nx, 1:fft_cf%n(2), k)
    end do

    if( (poisson_solver .eq. FFT_CYL) .and. (gr%sb%periodic_dim == 0) ) then
      deallocate(x, y)
    end if
#endif
  end subroutine poisson_fft_build_3d
    

  subroutine poisson_fft_build_2d(gr, poisson_solver)
    type(grid_t), intent(in) :: gr
    integer, intent(in) :: poisson_solver
#if defined(HAVE_FFT)

    type(loct_spline_t) :: besselintf
    integer :: i, ix, iy, ixx(MAX_DIM), db(MAX_DIM), npoints
    FLOAT :: temp(MAX_DIM), vec, r_c, maxf, dk
    FLOAT :: DELTA_R = CNST(1.0e-12)
    FLOAT, allocatable :: x(:), y(:)


    ! double the box to perform the fourier transforms
    if(poisson_solver.ne.FFT_CORRECTED) then
      call mesh_double_box(gr%sb, gr%m, db)                 ! get dimensions of the double box
      if (poisson_solver == FFT_SPH) db(1:2) = maxval(db)
    else
      db(:) = gr%m%l(:)
    end if

    call dcf_new(db, fft_cf)      ! allocate cube function where we will perform
    call dcf_fft_init(fft_cf, gr%sb) ! the ffts
    db = fft_cf%n                 ! dimensions may have been optimized

    call loct_parse_float(check_inp('PoissonCutoffRadius'),&
      maxval(db(:)*gr%m%h(:)/M_TWO)/units_inp%length%factor , r_c)
    r_c = r_c*units_inp%length%factor
    write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
      trim(units_out%length%abbrev), '] = ',       &
      r_c/units_out%length%factor
    call write_info(1)
    if ( r_c > maxval(db(:)*gr%m%h(:)/M_TWO) + DELTA_R) then
      message(1) = 'Poisson cutoff radius is larger than cell size.'
      message(2) = 'You can see electrons in next cell(s).'
      call write_warning(2)
    end if

    call loct_spline_init(besselintf)

    ! store the fourier transform of the Coulomb interaction
    ALLOCATE(fft_Coulb_FS(fft_cf%nx, fft_cf%n(2), fft_cf%n(3)), fft_cf%nx*fft_cf%n(2)*fft_cf%n(3))
    fft_Coulb_FS = M_ZERO

    temp(:) = M_TWO*M_PI/(db(:)*gr%m%h(:))

    maxf = r_c * sqrt((temp(1)*db(1)/2)**2 + (temp(2)*db(2)/2)**2)
    dk = CNST(0.3) ! This seems to be reasonable.
    npoints = nint(maxf/dk)
    ALLOCATE(x(npoints), npoints)
    ALLOCATE(y(npoints), npoints)
    do i = 1, npoints
       x(i) = (i-1) * maxf / (npoints-1)
       y(i) = besselint(x(i))
    end do
    call loct_spline_fit(npoints, x, y, besselintf)

    do iy = 1, db(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      do ix = 1, fft_cf%nx
        ixx(1) = pad_feq(ix, db(1), .true.)
        vec = sqrt( (temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
        fft_coulb_fs(ix, iy, 1) = M_TWO * M_PI * r_c * loct_splint(besselintf, vec*r_c)
      end do
    end do

    deallocate(x, y)
    call loct_spline_end(besselintf)
#endif
  end subroutine poisson_fft_build_2d

  !-----------------------------------------------------------------
  subroutine poisson_fft_end()
#if defined(HAVE_FFT)
      call dcf_free(fft_cf)
      deallocate(fft_coulb_FS); nullify(fft_coulb_FS)
#endif
  end subroutine poisson_fft_end

  !-----------------------------------------------------------------
  subroutine poisson_fft(m, pot, rho, average_to_zero)
    type(mesh_t), intent(in) :: m
    FLOAT, intent(out) :: pot(:) ! pot(m%np)
    FLOAT, intent(in)  :: rho(:) ! rho(m%np)
    logical, intent(in), optional :: average_to_zero

#if defined(HAVE_FFT)

    FLOAT, allocatable :: rho_global(:), pot_global(:)

    integer :: k, j, i
    FLOAT :: average

    call push_sub('poisson.poisson_fft')

    average=M_ZERO !this avoids a non-initialized warning
    
    if(m%parallel_in_domains) then
      ALLOCATE(rho_global(m%np_global), m%np_global)
      ALLOCATE(pot_global(m%np_global), m%np_global)
    end if

    call dcf_alloc_RS(fft_cf)          ! allocate the cube in real space
    call dcf_alloc_FS(fft_cf)          ! allocate the cube in Fourier space

    if(m%parallel_in_domains) then
#if defined HAVE_MPI
      call dvec_gather(m%vp, rho_global, rho)
      call dmf2cf(m, rho_global, fft_cf)        ! put the density in a cube
#endif
    else
      call dmf2cf(m, rho, fft_cf)        ! put the density in a cube
    end if
    call dcf_RS2FS(fft_cf)             ! Fourier transform


    ! multiply by the FS of the Coulomb interaction
    do k = 1, fft_cf%n(3)
      do j = 1, fft_cf%n(2)
        do i = 1, fft_cf%nx
          fft_cf%FS(i, j, k) = fft_cf%FS(i, j, k)*fft_Coulb_FS(i, j, k)
        end do
      end do
    end do

    call dcf_FS2RS(fft_cf)             ! Fourier transform back
    if(present(average_to_zero)) then
      if(average_to_zero) average = cf_surface_average(fft_cf)
#if defined HAVE_MPI
      ! Only root has the right average.
      if(m%parallel_in_domains) call MPI_Bcast(average, 1, MPI_FLOAT, 0, m%mpi_grp%comm, k)
#endif
    end if

    if(m%parallel_in_domains) then
#if defined(HAVE_MPI)
      call dcf2mf(m, fft_cf, pot_global)        ! put the density in a mesh
      call dvec_scatter(m%vp, pot_global, pot)
#endif
    else
      call dcf2mf(m, fft_cf, pot)        ! put the density in a mesh
    end if

    if(present(average_to_zero)) then
      if(average_to_zero) pot = pot - average
    end if

    call dcf_free_RS(fft_cf)           ! memory is no longer needed
    call dcf_free_FS(fft_cf)

    if(m%parallel_in_domains) then
      deallocate(rho_global, pot_global)
    end if
#endif
    call pop_sub()
  end subroutine poisson_fft



end module poisson_fft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
