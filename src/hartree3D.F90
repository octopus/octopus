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

subroutine hartree3D_init(h, m)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(inout) :: m

  select case(h%solver)
    case(1)
      message(1) = 'Info: Using conjugated gradients method to solve poisson equation.'
      call write_info(1)
      call init_real()
#ifdef HAVE_FFT
    case(2)
      message(1) = 'Info: Using FFTs to solve poisson equation.'
      call write_info(1)
      call init_fft()
    case(3)
      message(1) = 'Info: Using FFTs to solve poisson equation with spherical cutoff.'
      call write_info(1)
      call init_fft()
#endif
    case default
      message(1) = "Internal error in hartree3D_init"
      call write_fatal(1)
  end select

  call pop_sub()

contains
  subroutine init_real()
 
    ! setup mesh including ghost points
    allocate(h%m_aux)
    h%m_aux = m
    nullify(h%m_aux%Lxyz, h%m_aux%Lxyz_inv,  h%m_aux%der_lookup)
    call mesh_create_xyz(h%m_aux, m%d%norder)
    
  end subroutine init_real

#ifdef HAVE_FFT
  subroutine init_fft()
    integer :: ix, iy, iz, ixx(3), db(3), dbc(3)
    real(r8) :: r_0, temp(3), vec

    ! double the box to perform the fourier transforms
    call mesh_double_box(m, db)
    if(h%solver == 3) db(:) = maxval(db)

    ! initialize the ffts
    call fft_init(db, fft_real, h%fft)
    call fft_getdim_real   (h%fft, db)
    call fft_getdim_complex(h%fft, dbc)

    ! store the fourier transform of the Coulomb interaction
    allocate(h%ff(dbc(1), dbc(2), dbc(3)))
    h%ff = M_ZERO
    r_0 = maxval(db(:)*m%h(:))/M_TWO
    temp(:) = 2.0_r8*M_PI/(db(:)*m%h(:))
      
    do ix = 1, dbc(1)
      ixx(1) = pad_feq(ix, db(1), .true.)
      do iy = 1, dbc(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, dbc(3)
          ixx(3) = pad_feq(iz, db(3), .true.)

          vec = sum((temp(:)*ixx(:))**2) 
          if(vec.ne.M_ZERO) then
            if(h%solver==2) then
              h%ff(ix, iy, iz) = 4.0_r8*M_PI / vec
            else
              h%ff(ix, iy, iz) = 4.0_r8*M_Pi*(1.0_8 - cos(sqrt(vec)*r_0)) / vec 
            endif
          else
            if(h%solver==2) then
              h%ff(ix, iy, iz) = M_ZERO
            else
              h%ff(ix, iy, iz) = 2.0_r8*M_Pi*r_0**2 
            endif
          end if
        end do
      end do
    end do

  end subroutine init_fft
#endif

end subroutine hartree3D_init

subroutine hartree_cg(h, m, pot, dist)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(IN) :: m
  real(r8), dimension(:), intent(inout) :: pot
  real(r8), dimension(:, :), intent(IN) :: dist
  
  integer, parameter :: ml = 4 ! maximum multipole moment to include
  real(r8), parameter :: fpi = M_FOUR*M_PI

  integer :: i, iter, add_lm, l, mm
  real(r8) :: sa, x(3), r, s1, s2, s3, ak, ck, temp
  real(r8), allocatable :: zk(:), tk(:), pk(:), rholm(:), wk(:), lwk(:)

  call push_sub('hartree_cg')

  allocate(rholm((ml+1)**2))             ! to store the multipole moments
  allocate(zk(m%np), tk(m%np), pk(m%np)) ! for the cg

  ! calculate multipole moments until ml
  rholm = M_ZERO
  do i = 1, m%np
    call mesh_r(m, i, r, x=x)
    
    temp   = sum(dist(i, :))
    zk(i) = -fpi*temp ! this will be needed later

    add_lm = 1
    do l = 0, ml
      s1 = r**l*temp
      do mm = -l, l
        sa = oct_ylm(x(1), x(2), x(3), l, mm)
        rholm(add_lm) = rholm(add_lm) + s1*sa
        add_lm = add_lm+1
      end do
    end do
  end do
  rholm = rholm*m%vol_pp

  ! build initial guess for the potential
  allocate(wk(h%m_aux%np), lwk(h%m_aux%np))
  wk(1:m%np) = pot(1:m%np)
  do i = h%m_aux%np_in+1, h%m_aux%np ! boundary conditions
    call mesh_r(h%m_aux, i, r, x=x)
    
    add_lm = 1
    do l = 0, ml
      s1 = fpi/((M_TWO*l + M_ONE)*r**(l + 1))
      do mm = -l, l
        sa = oct_ylm(x(1), x(2), x(3), l, mm)
        wk(i) = wk(i) +  sa * rholm(add_lm) * s1
        add_lm = add_lm+1
      end do
    end do
  end do
  deallocate(rholm)

  call dmf_laplacian(h%m_aux, wk, lwk)

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

    if(iter >= 400 .or. abs(s3) < 1.0e-5_r8) exit

    ck = s3/s1
    s1 = s3
    pk = zk + ck*pk
  end do

  if(iter >= 400) then
    message(1) = "Hartree using conjugated gradients: Not converged!"
    write(message(2),*) "iter = ", iter, " s3 = ", s3
    call write_warning(2)
  endif

  h%ncgiter = iter

  deallocate(zk, tk, pk)

  call pop_sub()
  return
end subroutine hartree_cg

#if defined(HAVE_FFT)
subroutine hartree_fft(h, m, pot, dist)
  type(hartree_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  real(r8), dimension(:), intent(out) :: pot
  real(r8), dimension(:, :), intent(IN) :: dist

  integer :: i, db(3), dbc(3), c(3)
  real(r8), allocatable :: f(:,:,:)
  complex(r8), allocatable :: gg(:,:,:)

  call push_sub('hartree_fft')

  call fft_getdim_real   (h%fft, db)  ! get dimensions of real array
  call fft_getdim_complex(h%fft, dbc) ! get dimensions of complex array
  c(:) = db(:)/2 + 1                  ! get center of double box
  
  allocate(f(db(1), db(2), db(3)), gg(dbc(1), dbc(2), dbc(3)))

  f = M_ZERO
  do i = 1, m%np
    f(m%Lxyz(1, i)+c(1), m%Lxyz(2, i)+c(2), m%Lxyz(3, i)+c(3)) = sum(dist(i, :))
  end do
  
  ! Fourier transform the charge density
  
  call rfft_forward (h%fft, f, gg)
  gg = gg*h%ff
  call rfft_backward(h%fft, gg, f)
  
  do i = 1, m%np
    pot(i) = f(m%Lxyz(1, i) + c(1), m%Lxyz(2, i) + c(2), m%Lxyz(3, i) + c(3))
  end do

  deallocate(f, gg)

  call pop_sub()
  return
end subroutine hartree_fft
#endif
