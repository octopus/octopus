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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! ---------------------------------------------------------
subroutine poisson3D_init(gr, geo)
  type(grid_t), intent(inout) :: gr
  type(geometry_t), intent(in) :: geo

  call push_sub('poisson3D.poisson3D_init')

  ASSERT(poisson_solver >= FFT_SPH .or. poisson_solver <= MULTIGRILLA)

  !%Variable PoissonSolverMaxMultipole
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Description
  !% Order of the multipolar expansion for boundary
  !% corrections. Default is 4 for cg_corrected and multigrid and 2
  !% for fft_corrected.
  !%End


  !%Variable PoissonSolverThreshold
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Description
  !% The tolerance for the poisson solution, used by the cg and
  !% multigrid solvers. Default is <math>10^{-5}</math>.
  !%End

  select case(poisson_solver)
  case(CG)
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 4, maxl)
     write(message(1),'(a,i2)')'Info: Boundary conditions fixed up to L =',  maxl
     call write_info(1)
     call loct_parse_float(check_inp('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     call poisson_cg1_init(gr%m, maxl, threshold)

  case(CG_CORRECTED)
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 4, maxl)
     call loct_parse_float(check_inp('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call write_info(1)
     call poisson_cg2_init(gr%m, maxl, threshold)

  case(MULTIGRILLA)
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 4, maxl)
     call loct_parse_float(check_inp('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call write_info(1)

     call poisson_multigrid_init(gr%m, maxl, threshold)

     call grid_create_multigrid(gr, geo)
  end select

#ifdef HAVE_FFT
  if (poisson_solver <= FFT_CORRECTED) call init_fft(gr%m)

  if (poisson_solver == FFT_CORRECTED) then
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 2, maxl)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call write_info(1)
     call build_aux(gr%m)
     call build_phi(gr%m)
  end if
#endif

  call pop_sub()

contains

#ifdef HAVE_FFT
  ! ---------------------------------------------------------
  subroutine init_fft(m)
    type(mesh_t), intent(in) :: m

    type(loct_spline_t) :: cylinder_cutoff_f
    FLOAT, allocatable :: x(:), y(:)
    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), k, ngp
    FLOAT :: temp(MAX_DIM), modg2, xmax
    FLOAT :: gpar, gperp, gx, gy, gz, r_c
    FLOAT :: DELTA_R = CNST(1.0e-12)

    ! double the box to perform the fourier transforms
    if(poisson_solver.ne.FFT_CORRECTED) then
       call mesh_double_box(gr%sb, m, db)                 ! get dimensions of the double box
       if (poisson_solver == FFT_SPH) db(:) = maxval(db)
    else
       db(:) = m%l(:)
    end if

    call dcf_new(db, fft_cf)    ! allocate cube function where we will perform
    call dcf_fft_init(fft_cf, gr%sb)   ! the ffts
    db = fft_cf%n               ! dimensions may have been optimized

    if (poisson_solver <= FFT_PLA .and. poisson_solver .ne. FFT_CORRECTED) then
       call loct_parse_float(check_inp('PoissonCutoffRadius'),&
            maxval(db(:)*m%h(:)/M_TWO)/units_inp%length%factor , r_c)
       r_c = r_c*units_inp%length%factor
       write(message(1),'(3a,f12.6)')'Info: Poisson Cutoff Radius [',  &
            trim(units_out%length%abbrev), '] = ',       &
            r_c/units_out%length%factor
       call write_info(1)
       if ( r_c > maxval(db(:)*m%h(:)/M_TWO) + DELTA_R) then
          message(1) = 'Poisson cutoff radius is larger than cell size.'
          message(2) = 'You can see electrons in next cell(s).'
          call write_warning(2)
       end if
    end if

    ! store the fourier transform of the Coulomb interaction
    ALLOCATE(fft_Coulb_FS(fft_cf%nx, fft_cf%n(2), fft_cf%n(3)), fft_cf%nx*fft_cf%n(2)*fft_cf%n(3))
    fft_Coulb_FS = M_ZERO

    temp(:) = M_TWO*M_PI/(db(:)*m%h(:))

    if(poisson_solver .eq. FFT_CYL) then
      ngp = 8*db(2)
      ALLOCATE(x(ngp), ngp)
      ALLOCATE(y(ngp), ngp)
    end if

    do ix = 1, fft_cf%nx
      ixx(1) = pad_feq(ix, db(1), .true.)
      if(poisson_solver .eq. FFT_CYL) then
        call loct_spline_init(cylinder_cutoff_f)
        gx = temp(1)*ixx(1)
        xmax = sqrt((temp(2)*db(2)/2)**2 + (temp(3)*db(3)/2)**2)
        do k = 1, ngp
          x(k) = (k-1)*(xmax/(ngp-1))
          y(k) = loct_poisson_finite_cylinder(gx, x(k), M_TWO*m%sb%xsize, M_TWO*m%sb%rsize)
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
                   fft_Coulb_FS(ix, iy, iz) = cutoff0(sqrt(modg2),r_c)/modg2

                case(FFT_CYL)
                   gperp = sqrt((temp(2)*ixx(2))**2+(temp(3)*ixx(3))**2)
                   if (gr%sb%periodic_dim==1) then
                     fft_Coulb_FS(ix, iy, iz) = cutoff1(abs(gx), gperp, r_c)/modg2
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
                   fft_Coulb_FS(ix, iy, iz) = cutoff2(gpar,gz,r_c)/modg2

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
                    fft_Coulb_FS(ix, iy, iz) = loct_poisson_finite_cylinder(M_ZERO, M_ZERO, M_TWO*m%sb%xsize, M_TWO*m%sb%rsize)
                  end if

                case(FFT_PLA)
                  fft_Coulb_FS(ix, iy, iz) = -M_HALF*r_c**2

                case (FFT_NOCUT, FFT_CORRECTED)
                  fft_Coulb_FS(ix, iy, iz) = M_ZERO
                end select
             end if
          end do
       end do

      if(poisson_solver .eq. FFT_CYL) then
        call loct_spline_end(cylinder_cutoff_f)
      end if
    end do

    fft_Coulb_FS(:,:,:) = M_FOUR*M_PI*fft_Coulb_FS(:,:,:)

    if(poisson_solver .eq. FFT_CYL) then
      deallocate(x, y)
    end if

  end subroutine init_fft
#endif

end subroutine poisson3D_init



