!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id$

! ---------------------------------------------------------
subroutine poisson3D_init(gr)
  type(grid_t), intent(inout) :: gr

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
     call loct_parse_float(check_inp('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     call poisson_cg1_init(gr%m, maxl, threshold)

  case(CG_CORRECTED)
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 4, maxl)
     call loct_parse_float(check_inp('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     call poisson_cg2_init(gr%m, maxl, threshold)

  case(MULTIGRILLA)
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 4, maxl)
     call loct_parse_float(check_inp('PoissonSolverThreshold'), CNST(1.0e-6), threshold)

     call poisson_multigrid_init(gr%m, maxl, threshold)

     call grid_create_multigrid(gr)
  end select

#ifdef HAVE_FFT
  if (poisson_solver <= FFT_CORRECTED) call init_fft(gr%m)

  if (poisson_solver == FFT_CORRECTED) then
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 2, maxl)
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

    integer :: ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM)
    FLOAT :: temp(MAX_DIM), modg2
    FLOAT :: gpar,gperp,gx,gz,r_c
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

    do iz = 1, db(3)
       ixx(3) = pad_feq(iz, db(3), .true.)
       do iy = 1, db(2)
          ixx(2) = pad_feq(iy, db(2), .true.)
          do ix = 1, fft_cf%nx
             ixx(1) = pad_feq(ix, db(1), .true.)

             modg2 = sum((temp(:)*ixx(:))**2)
             if(modg2 /= M_ZERO) then
                select case(poisson_solver)
                case(FFT_SPH)
                   fft_Coulb_FS(ix, iy, iz) = cutoff0(sqrt(modg2),r_c)/modg2
                case(FFT_CYL)
                   gx = abs(temp(1)*ixx(1))
                   gperp = sqrt((temp(2)*ixx(2))**2+(temp(3)*ixx(3))**2)
                   fft_Coulb_FS(ix, iy, iz) = cutoff1(gx,gperp,r_c)/modg2
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
                   fft_Coulb_FS(ix, iy, iz) = -(M_HALF*log(r_c)-M_FOURTH)*r_c**2
                case(FFT_PLA)
                   fft_Coulb_FS(ix, iy, iz) = -M_HALF*r_c**2
                case (FFT_NOCUT, FFT_CORRECTED)
                   fft_Coulb_FS(ix, iy, iz) = M_ZERO
                end select
             end if
             fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
          end do
       end do
    end do

  end subroutine init_fft
#endif

end subroutine poisson3D_init



