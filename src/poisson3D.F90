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

subroutine poisson3D_init(m, geo)
  type(mesh_type),     intent(inout) :: m
  type(geometry_type), intent(IN)    :: geo ! this is needed to call mesh_create_xyz


  ASSERT(poisson_solver >= 0 .or. poisson_solver <= 4)

  select case(poisson_solver)
#ifdef HAVE_FFT
    case(FFT_SPH)
      message(1) = 'Info: Using FFTs with spherical cutoff to solve Poisson equation .'
    case(FFT_CYL)
      message(1) = 'Info: Using FFTs with cylindrical cutoff to solve Poisson equation .'
    case(FFT_PLA)
      message(1) = 'Info: Using FFTs with planar cutoff to solve Poisson equation .'
    case(FFT_NOCUT)
      message(1) = 'Info: Using FFTs without cutoff to solve Poisson equation.'
#endif
    case(CG)
      message(1) = 'Info: Using conjugated gradients method to solve poisson equation.'
      call init_real()
  end select
  call write_info(1)

#ifdef HAVE_FFT
  if (poisson_solver <= 3) call init_fft()
#endif

  call pop_sub()

contains
  subroutine init_real()
 
    ! setup mesh including ghost points
    allocate(cg_m_aux)
    cg_m_aux = m
    nullify(cg_m_aux%lxyz, cg_m_aux%lxyz_inv, cg_m_aux%grad)
    call mesh_create_xyz(cg_m_aux, geo, cg_m_aux%der_order)
    call mf_add_der(cg_m_aux)

  end subroutine init_real

#ifdef HAVE_FFT
  subroutine init_fft()
    integer :: ix, iy, iz, ixx(3), db(3)
    FLOAT :: temp(3), vec
    FLOAT :: gpar,gperp,gx,gz,r_c
    FLOAT :: DELTA_R = CNST(1.0e-12)

    ! double the box to perform the fourier transforms
    call mesh_double_box(m, db)                 ! get dimensions of the double box
    if (poisson_solver == FFT_SPH) db(:) = maxval(db)

    call dcf_new(db, fft_cf)    ! allocate cube function where we will perform
    call dcf_fft_init(fft_cf)   ! the ffts
    db = fft_cf%n               ! dimensions may have been optimized

  if (poisson_solver <= 2) then
    call loct_parse_float('PoissonCutoffRadius',&
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
    allocate(fft_Coulb_FS(fft_cf%nx, fft_cf%n(2), fft_cf%n(3)))
    fft_Coulb_FS = M_ZERO

    temp(:) = M_TWO*M_PI/(db(:)*m%h(:))
      
    do iz = 1, db(3)
      ixx(3) = pad_feq(iz, db(3), .true.)
      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do ix = 1, fft_cf%nx
          ixx(1) = pad_feq(ix, db(1), .true.)

           vec = sum((temp(:)*ixx(:))**2)
           if(vec /= M_ZERO) then
             select case(poisson_solver)
               case(FFT_SPH)
                 fft_Coulb_FS(ix, iy, iz) = cutoff0(sqrt(vec)*r_c)/vec
               case(FFT_CYL)
                 gx = abs(temp(1)*ixx(1))
                 gperp = sqrt((temp(2)*ixx(2))**2+(temp(3)*ixx(3))**2)
                 fft_Coulb_FS(ix, iy, iz) = cutoff1(gx*r_c,gperp*r_c)/vec
               case(FFT_PLA)
                 gz = abs(temp(3)*ixx(3))
                 gpar = sqrt((temp(1)*ixx(1))**2+(temp(2)*ixx(2))**2)
                 fft_Coulb_FS(ix, iy, iz) = cutoff2(gpar*r_c,gz*r_c)/vec
               case(FFT_NOCUT)
                 fft_Coulb_FS(ix, iy, iz) = M_ONE/vec
              end select               
            else
              select case(poisson_solver)
              case(FFT_SPH)
                fft_Coulb_FS(ix, iy, iz) = r_c**2/M_TWO
              case (FFT_CYL,FFT_PLA,FFT_NOCUT)
                fft_Coulb_FS(ix, iy, iz) = M_ZERO
              end select
            endif
            fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
        end do
      end do
    end do

  end subroutine init_fft
#endif

end subroutine poisson3D_init

subroutine poisson_cg(m, pot, rho)
  type(mesh_type), intent(IN) :: m
  FLOAT, intent(inout) :: pot(m%np)
  FLOAT, intent(in)    :: rho(m%np)
  
  integer, parameter :: ml = 4 ! maximum multipole moment to include
  FLOAT, parameter :: fpi = M_FOUR*M_PI

  integer :: i, iter, add_lm, l, mm
  FLOAT :: sa, x(3), r, s1, s2, s3, ak, ck
  FLOAT, allocatable :: zk(:), tk(:), pk(:), rholm(:), wk(:), lwk(:)

  call push_sub('poisson_cg')

  allocate(rholm((ml+1)**2))             ! to store the multipole moments
  allocate(zk(m%np), tk(m%np), pk(m%np)) ! for the cg

  ! calculate multipole moments until ml
  rholm = M_ZERO
  do i = 1, m%np
    call mesh_r(m, i, r, x=x)
    
    zk(i) = -fpi*rho(i) ! this will be needed later

    add_lm = 1
    do l = 0, ml
      if(l == 0) then 
        s1 = r**l*rho(i)
      else
        s1 = rho(i)
      end if

      do mm = -l, l
        sa = loct_ylm(x(1), x(2), x(3), l, mm)
        rholm(add_lm) = rholm(add_lm) + s1*sa
        add_lm = add_lm+1
      end do
    end do
  end do
  rholm = rholm*m%vol_pp

  ! build initial guess for the potential
  allocate(wk(cg_m_aux%np), lwk(cg_m_aux%np))
  wk(1:m%np) = pot(1:m%np)
  do i = cg_m_aux%np_in+1, cg_m_aux%np ! boundary conditions
    call mesh_r(cg_m_aux, i, r, x=x)
    
    add_lm = 1
    do l = 0, ml
      s1 = fpi/((M_TWO*l + M_ONE)*r**(l + 1))
      do mm = -l, l
        sa = loct_ylm(x(1), x(2), x(3), l, mm)
        wk(i) = wk(i) +  sa * rholm(add_lm) * s1
        add_lm = add_lm+1
      end do
    end do
  end do
  deallocate(rholm)

  call dmf_laplacian(cg_m_aux, wk, lwk)

  zk(:) = zk(:) - lwk(1:m%np)
  deallocate(wk, lwk) ! they are no longer needed

  pk = zk
  s1 = dmf_nrm2(m, zk)**2

  ! now we start the conjugate gradient cycles
  iter = 0
  do 
    iter = iter + 1
    call dmf_laplacian(m, pk, tk)

    s2 = dmf_dotp(m, zk, tk)
    ak = s1/s2
    pot = pot + ak*pk
    zk = zk - ak*tk
    s3 = dmf_nrm2(m, zk)**2

    if(iter >= 400 .or. abs(s3) < CNST(1.0e-5)) exit

    ck = s3/s1
    s1 = s3
    pk = zk + ck*pk
  end do

  if(iter >= 400) then
    message(1) = "Poisson using conjugated gradients: Not converged!"
    write(message(2),*) "iter = ", iter, " s3 = ", s3
    call write_warning(2)
  endif

  deallocate(zk, tk, pk)

  call pop_sub()
  return
end subroutine poisson_cg

#if defined(HAVE_FFT)
subroutine poisson_fft(m, pot, rho)
  type(mesh_type), intent(IN) :: m
  FLOAT, intent(out) :: pot(m%np)
  FLOAT, intent(in)  :: rho(m%np)

  call push_sub('poisson_fft')

  call dcf_alloc_RS(fft_cf)          ! allocate the cube in real space
  call dcf_alloc_FS(fft_cf)          ! allocate the cube in Fourier space

  call dmf2cf(m, rho, fft_cf)        ! put the density in a cube
  call dcf_RS2FS(fft_cf)             ! Fourier transform
  fft_cf%FS = fft_cf%FS*fft_Coulb_FS ! multiply by the FS of the Coulomb interaction

  call dcf_FS2RS(fft_cf)             ! Fourier transform back
  call dcf2mf(m, fft_cf, pot)        ! put the density in a cube
  
  call dcf_free_RS(fft_cf)           ! memory is no longer needed
  call dcf_free_FS(fft_cf)

  call pop_sub()
  return
end subroutine poisson_fft
#endif
