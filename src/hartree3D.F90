module hartree
use global
use liboct
use mesh
use math
#ifdef HAVE_FFTW
use fft
#endif
implicit none

private

type hartree_type
  integer :: solver

  ! used by conjugated gradients (method 1)
  integer :: ncgiter

#ifdef HAVE_FFTW
  ! used in method 3
  real(r8), pointer :: ff(:, :, :)   
#endif
end type hartree_type

public :: hartree_type, hartree_init, hartree_solve, hartree_end

contains

subroutine hartree_init(h, m)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(inout) :: m

  integer :: ix, iy, iz, ixx(3)
  real(r8) :: l, r_0, temp(3), vec

  sub_name = 'hartree_init'; call push_sub()

  call oct_parse_int(C_string('PoissonSolver'), 3, h%solver)
  if(h%solver<1 .or. h%solver>3 ) then
    write(message(1), '(a,i2,a)') "Input: '", h%solver, &
         "' is not a valid PoissonSolver"
    message(2) = 'PoissonSolver = 1(cg) | 2(fft) | 3(fft spherical cutoff)'
    call write_fatal(2)
  end if
  if(h%solver.eq.3 .and. m%box_shape.eq.4 ) then
    write(message(1), '(a,i2,a)') "Input: '", h%solver, "'"
    message(2) = "  is not a valid PoissonSolver when parallelpiped box-shape is used."
    message(3) = 'PoissonSolver = 1(cg)'
    call write_fatal(3)
  end if
#ifdef POLYMERS
  if(h%solver.ne.2) then
    message(1) = "Sorry, but for polymers only the fft method (method 2)"
    message(2) = "is available for solving the poisson equation"
    call write_fatal(2)
  end if
#endif


  select case(h%solver)
  case(1)
    message(1) = 'Info: Using conjugated gradients method to solve poisson equation.'
    call write_info(1)
    
#ifdef HAVE_FFTW
  case(2)
    message(1) = 'Info: Using FFTs to solve poisson equation.'
    call write_info(1)
    allocate(h%ff(m%hfft_n2, m%fft_n2(2), m%fft_n2(3)))

    temp(:) = 2.0_r8*M_PI/(m%fft_n2(:)*m%h(:))
      
    do ix = 1, m%hfft_n2
      ixx(1) = pad_feq(ix, m%fft_n2(1), .true.)
      do iy = 1, m%fft_n2(2)
        ixx(2) = pad_feq(iy, m%fft_n2(2), .true.)
        do iz = 1, m%fft_n2(3)
          ixx(3) = pad_feq(iz, m%fft_n2(3), .true.)
          vec = sum((temp(:)*ixx(:))**2)

          if(vec.ne.0.0_r8) then
            h%ff(ix, iy, iz) = 4.0_r8*M_PI / vec
          else
            h%ff(ix, iy, iz) = 0._r8 
          end if          
        end do
      end do
    end do

    call mesh_alloc_ffts(m, 2)    
  case(3) ! setup spherical ffts
    message(1) = 'Info: Using FFTs to solve poisson equation with spherical cutoff.'
    call write_info(1)

    allocate(h%ff(m%hfft_n2, m%fft_n2(1), m%fft_n2(1)))
    h%ff = 0.0_r8

    l = m%fft_n2(1)*m%h(1)
    r_0 = l/2.0_r8
    temp(1) = (2.0_r8*M_PI/l)
      
    do ix = 1, m%hfft_n2
      ixx(1) = pad_feq(ix, m%fft_n2(1), .true.)
      do iy = 1, m%fft_n2(1)
        ixx(2) = pad_feq(iy, m%fft_n2(2), .true.)
        do iz = 1, m%fft_n2(1)
          ixx(3) = pad_feq(iz, m%fft_n2(3), .true.)
          vec = temp(1)*sqrt(real(ixx(1)**2 + ixx(2)**2 + ixx(3)**2, r8))
          if(vec.ne.0.0_r8) then
            h%ff(ix, iy, iz) = 4.0_r8*M_Pi*(1.0_8 - cos(vec*r_0)) / vec**2 
          else
            h%ff(ix, iy, iz) = 2.0_r8*M_Pi*r_0**2 
          end if
        end do
      end do
    end do

    call mesh_alloc_ffts(m, 2)
#endif
  end select

  call pop_sub()
  return
end subroutine hartree_init

subroutine hartree_end(h)
  type(hartree_type), intent(inout) :: h

  sub_name = 'hartree_end'; call push_sub()

  select case(h%solver)
  case(1)
#ifdef HAVE_FFTW
  case(2,3)
    if(associated(h%ff)) then ! has been allocated => destroy
      deallocate(h%ff); nullify(h%ff)
    end if
#endif
  end select

  call pop_sub()
  return
end subroutine hartree_end

subroutine hartree_solve(h, m, pot, dist)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(IN) :: m
  real(r8), dimension(:), intent(inout) :: pot
  real(r8), dimension(:,:), intent(IN) :: dist

  sub_name = 'hartree_solve'; call push_sub()

  select case(h%solver)
  case(1)
    call hartree_cg(h, m, pot, dist)

#ifdef HAVE_FFTW
  case(2,3)
    call hartree_fft(h, m, pot, dist)
#endif

  case default
    message(1) = "Hartree structure not initialized"
    call write_fatal(1)
  end select

  call pop_sub()
  return
end subroutine hartree_solve

subroutine hartree_cg(h, m, pot, dist)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(IN) :: m
  real(r8), dimension(:), intent(inout) :: pot
  real(r8), dimension(:,:), intent(IN) :: dist
  
  integer :: i, ix, iy, iz, iter, a, j, imax, add_lm, l, mm, ml
  real(r8) :: sa, x(3), r, s1, s2, s3, ak, ck, xx, yy, zz, temp
  real(r8), allocatable :: zk(:), tk(:), pk(:), rholm(:), wk(:,:,:)

  sub_name = 'hartree_cg'; call push_sub()

  ml = 4
  allocate(zk(m%np), tk(m%np), pk(m%np))
  allocate(rholm((ml+1)**2))
  rholm = 0.0_r8

  ! calculate multipole moments till ml
  do i = 1, m%np
    call mesh_r(m, i, r, x=x)

    temp   = sum(dist(i, :))
    add_lm = 1
    do l = 0, ml
      do mm = -l, l
        sa = oct_ylm(x(1), x(2), x(3), l, mm)
        if(l == 0) then ! this avoids 0**0
          rholm(add_lm) = rholm(add_lm) + temp*sa
        else
          rholm(add_lm) = rholm(add_lm) + r**l*temp*sa
        end if
        add_lm = add_lm+1
      end do
    end do
  end do
  rholm = rholm*m%vol_pp

  allocate(wk(-m%nx(1):m%nx(1), -m%nx(2):m%nx(2), -m%nx(3):m%nx(3)))

  wk = 0.0d0
  do i = 1, m%nk
    x(1) = m%kx(i)*m%H(1); x(2) = m%ky(i)*m%H(2); x(3) = m%kz(i)*m%H(3)
    r = sqrt(sum(x**2)) + 1e-50_r8

    add_lm = 1
    do l = 0, ml
      do mm = -l, l
        sa = oct_ylm(x(1), x(2), x(3), l, mm)

        wk(m%Kx(i), m%Ky(i), m%Kz(i)) = wk(m%Kx(i), m%Ky(i), m%Kz(i)) +   &
            sa * ((4.0_r8*M_Pi)/(2.0_r8*l + 1.0_r8)) *            &
            rholm(add_lm)/r**(l + 1)
        add_lm = add_lm+1
      enddo
    enddo
  enddo

  do i = 1, m%np
    zk(i) = -4.0d0*M_Pi*sum(dist(i,:))
    wk(m%Lx(i), m%Ly(i), m%Lz(i)) = pot(i)
  enddo 

  ! I think this just calculates the laplacian of pot(i) -> zk?
  ! should not this call derivatives?
  do i = 1, m%np
    ix = m%Lx(i); iy = m%Ly(i); iz = m%Lz(i)
    !zk(i) = zk(i) - 3.0_r8*m%d%dlidfj(0)*wk(ix, iy, iz)/m%h**2
    zk(i) = zk(i) - m%d%dlidfj(0)*wk(ix,iy,iz)*(1/m%h(1)**2+1/m%h(2)**2+1/m%h(3)**2)
    do j = 1, m%d%norder
      zk(i) = zk(i) -  m%d%dlidfj(j)*(  &
          (wk(ix+j, iy, iz) + wk(ix-j, iy, iz))/m%h(1)**2 + &
          (wk(ix, iy+j, iz) + wk(ix, iy-j, iz))/m%h(2)**2 + &
          (wk(ix, iy, iz+j) + wk(ix, iy, iz-j))/m%h(3)**2 )
    end do
  end do

  wk = 0.d0
  pk = zk
  s1 = dmesh_nrm2(m, zk)**2

  iter = 0
  do 
    iter = iter + 1
    do i = 1, m%np
      wk(m%Lx(i), m%Ly(i), m%Lz(i)) = pk(i)
    enddo
    ! again the laplacian of wk -> tk
    do i = 1, m%np
      ix = m%Lx(i); iy = m%Ly(i); iz = m%Lz(i)
      !tk(i) = 3.0_r8*m%d%dlidfj(0)*wk(ix, iy, iz)
      tk(i) = m%d%dlidfj(0)*wk(ix,iy,iz)*(1/m%h(1)**2+1/m%h(2)**2+1/m%h(3)**2)
      do j = 1, m%d%norder
        tk(i) = tk(i) + m%d%dlidfj(j)*( &
            (wk(ix+j, iy, iz) + wk(ix-j, iy, iz))/m%h(1)**2 + &
            (wk(ix, iy+j, iz) + wk(ix, iy-j, iz))/m%h(2)**2 + &
            (wk(ix, iy, iz+j) + wk(ix, iy, iz-j))/m%h(3)**2 )
      enddo
      !tk(i) = tk(i) /m%h**2
    enddo

    s2 = dmesh_dotp(m, zk, tk)
    ak = s1/s2
    pot = pot + ak*pk
    zk = zk - ak*tk
    s3 = dmesh_nrm2(m, zk)**2

    if(iter >= 400 .or. abs(s3) < 1.0E-5_r8) exit

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

  deallocate(wk, zk, tk, pk, rholm)

  call pop_sub()
  return
end subroutine hartree_cg

#ifdef HAVE_FFTW
subroutine hartree_fft(h, m, pot, dist)
  type(hartree_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  real(r8), dimension(:), intent(out) :: pot
  real(r8), dimension(:,:), intent(IN) :: dist

  integer :: i, ix, iy, iz, s(3), c(3)
  real(r8), allocatable :: f(:,:,:)
  complex(r8), allocatable :: gg(:,:,:)

  sub_name = 'hartree_fft'; call push_sub()

  if(h%solver == 3) then
    s(:) = m%fft_n2(1)
  else
    s(:) = m%fft_n2(:)
  end if
  c(:) = s(:)/2 + 1
  
  allocate(f(s(1), s(2), s(3)), gg(s(1)/2+1, s(2), s(3)))
  f = 0.0_8
  
  do i = 1, m%np
    ix = m%Lx(i) + c(1)
    iy = m%Ly(i) + c(2)
    iz = m%Lz(i) + c(3)
    f(ix, iy, iz) = sum(dist(i,:))
  end do
  
  ! Fourier transform the charge density
  call rfft3d(f, gg, s(1), s(2), s(3), m%dplanf2, fftw_forward)
  gg = gg*h%ff
  call rfft3d(f, gg, s(1), s(2), s(3), m%dplanb2, fftw_backward)
  
  do i = 1, m%np
    ix = m%Lx(i) + c(1)
    iy = m%Ly(i) + c(2)
    iz = m%Lz(i) + c(3)
    pot(i) = f(ix, iy, iz)
  end do

  deallocate(f, gg)

  call pop_sub()
  return
end subroutine hartree_fft
#endif

end module hartree
