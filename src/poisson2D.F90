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

subroutine poisson2D_init(gr)
  type(grid_type), intent(inout) :: gr

  ASSERT(poisson_solver == FFT_SPH .or. poisson_solver == DIRECT_SUM_2D)

#ifdef HAVE_FFT
  if (poisson_solver == FFT_SPH) call init_fft(gr%m)
#endif

  call pop_sub()
contains

#ifdef HAVE_FFT
  subroutine init_fft(m)
    type(mesh_type), intent(in) :: m

    integer :: ix, iy, ixx(3), db(3)
    FLOAT :: temp(3), vec, r_c
    FLOAT :: DELTA_R = CNST(1.0e-12)

    ! double the box to perform the fourier transforms
    if(poisson_solver.ne.FFT_CORRECTED) then
       call mesh_double_box(gr%sb, gr%m, db)                 ! get dimensions of the double box
       if (poisson_solver == FFT_SPH) db(1:2) = maxval(db)
    else
       db(:) = m%l(:)
    endif

    call dcf_new(db, fft_cf)      ! allocate cube function where we will perform
    call dcf_fft_init(fft_cf, gr%sb) ! the ffts
    db = fft_cf%n                 ! dimensions may have been optimized


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

    ! store the fourier transform of the Coulomb interaction
    allocate(fft_Coulb_FS(fft_cf%nx, fft_cf%n(2), fft_cf%n(3)))
    fft_Coulb_FS = M_ZERO

    temp(:) = M_TWO*M_PI/(db(:)*m%h(:))

    do iy = 1, db(2)
       ixx(2) = pad_feq(iy, db(2), .true.)
       do ix = 1, fft_cf%nx
          ixx(1) = pad_feq(ix, db(1), .true.)
          vec = sqrt((temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
          fft_coulb_fs(ix, iy, 1) = r_c*besselint(vec*r_c)
       end do
    end do

  end subroutine init_fft
#endif

end subroutine poisson2D_init


subroutine poisson2D_solve(m, pot, rho)
  type(mesh_type), intent(in) :: m
  FLOAT, intent(out)          :: pot(m%np)
  FLOAT, intent(in)           :: rho(m%np)

  integer  :: i, j
  FLOAT    :: x(2), y(2)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
  FLOAT    :: tmp
  FLOAT, allocatable :: pvec(:)
#endif

  ASSERT(poisson_solver == -2)

  call push_sub('poisson2D.poisson2D_solve')

#if defined(HAVE_MPI) && defined(HAVE_METIS)

  allocate(pvec(1:m%np))

  pot = M_ZERO
  do i = 1, m%np_global
     x(:) = m%x_global(i,:)
     do j = 1, m%np
        if(m%vp%global(i, m%vp%partno) == j) then
           pvec(j) = M_TWO*sqrt(M_PI)*rho(j)/m%h(1)
        else
           y(:) = m%x(j,:)
           pvec(j) = rho(j)/sqrt(sum((x-y)**2))
        endif
     enddo
     tmp = dmf_integrate(m, pvec)
     if (m%part(i).eq.m%vp%partno) then
        pot(m%vp%global(i, m%vp%partno)) = tmp
     endif
  end do

  deallocate(pvec)

#else

  pot = M_ZERO
  do i = 1, m%np
     x(:) = m%x(i,:)
     do j = 1, m%np
        if(i == j) then
           pot(i) = pot(i) + M_TWO*sqrt(M_PI)*rho(i)/m%h(1)*m%vol_pp(j)
        else
           y(:) = m%x(j,:)
           pot(i) = pot(i) + rho(j)/sqrt(sum((x-y)**2))*m%vol_pp(j)
        end if
     end do
  end do

#endif

  call pop_sub()
end subroutine poisson2D_solve
