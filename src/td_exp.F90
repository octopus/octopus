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

subroutine td_dtexp(h, sys, td, ik, zpsi, timestep, t, &
                    order)
  type(hamiltonian_type), intent(in) :: h
  type(system_type), intent(in)      :: sys
  type(td_type), intent(in)          :: td
  integer, intent(in) :: ik
  complex(r8), intent(inout) :: zpsi(sys%m%np, sys%st%dim)
  real(r8), intent(in) :: timestep, t
  integer, optional, intent(out) :: order ! For the methods that rely on Hamiltonian-vector
                                          ! multipolication, the number of these.

  call push_sub('td_dtexp')

  select case(td%exp_method)
   case(FOURTH_ORDER);      call fourth
   case(LANCZOS_EXPANSION); call lanczos
#ifdef HAVE_FFT
   case(SPLIT_OPERATOR);    call split
   case(SUZUKI_TROTTER);    call suzuki
#endif
   case(CHEBYSHEV);         call cheby
  end select

  call pop_sub()
contains

  subroutine fourth
    complex(r8) :: zfact
    complex(r8), allocatable :: zpsi1(:,:), hzpsi1(:,:)
    integer :: i

    call push_sub('fourth')

    allocate(zpsi1(sys%m%np, sys%st%dim), hzpsi1(sys%m%np, sys%st%dim))
    zfact = M_z1
    call zcopy(sys%m%np*sys%st%dim, zpsi, 1, zpsi1, 1)
    do i = 1, td%exp_order
      zfact = zfact*(-M_zI*timestep)/i
      call zHpsi(h, sys%m, sys%st, sys, ik, zpsi1, hzpsi1, t)
      call zaxpy(sys%m%np*sys%st%dim, zfact, hzpsi1, 1, zpsi, 1)
      if(i .ne. td%exp_order) call zcopy(sys%m%np*sys%st%dim, hzpsi1, 1, zpsi1, 1)
    end do
    deallocate(zpsi1, hzpsi1)

    if(present(order)) order = td%exp_order
    call pop_sub()
  end subroutine fourth

  subroutine cheby
    ! /* Calculates the exponential of the hamiltonian through a expansion in Chebyshev's polynomials.
    ! For that purposes it uses the closed form of the coefficients[1] and Clenshaw-Gordons's[2]
    ! recursive algorithm.
    ! [1] H. Tal-Ezer and R. Kosloff, J. Chem. Phys 81, 3967 (1984).
    ! [2] C. W. Clenshaw, MTAC 9, 118 (1955).
    ! Since I don't have access to MTAC, I copied next Maple algorithm from Dr. F. G. Lether's 
    ! (University of Georgia) homepage: (www.math.uga.edu/~fglether):
    ! {twot := t + t; u0 := 0; u1 := 0;
    !  for k from n to 0 by -1 do
    !   u2 := u1; u1 := u0;
    !   u0 := twot*u1 - u2 + c[k];
    !  od;
    !  ChebySum := 0.5*(u0 - u2);} */
    integer :: j, n
    complex(r8) :: zfact
    complex(r8), allocatable :: zpsi1(:,:,:)

    call push_sub('cheby')
    
    allocate(zpsi1(sys%m%np, sys%st%dim, 0:2))
    zpsi1 = M_z0
    n = sys%m%np*sys%st%dim
    do j = td%exp_order, 0, -1
       call zcopy(n, zpsi1(1, 1, 1), 1, zpsi1(1, 1, 2), 1)
       call zcopy(n, zpsi1(1, 1, 0), 1, zpsi1(1, 1, 1), 1)
       call zhpsi(h, sys%m, sys%st, sys, ik, zpsi1(:, :, 1), zpsi1(1:, :, 0), t)
            zfact = 2*(-M_zI)**j*oct_bessel(j, h%spectral_half_span*timestep)
       call zaxpy(n, cmplx(-h%spectral_middle_point, 0.0_r8, r8), &
                                                                     zpsi1(1, 1, 1), 1, zpsi1(1, 1, 0), 1)
       call zscal(n, cmplx(1./h%spectral_half_span,0._r8, r8),    zpsi1(1, 1, 0), 1)
       call zscal(n, cmplx(2._r8,0._r8, r8),                      zpsi1(1, 1, 0), 1)
       call zaxpy(n, zfact, zpsi(1, 1),                      1,   zpsi1(1, 1, 0), 1)
       call zaxpy(n, cmplx(-1._r8,0._r8,r8), zpsi1(1, 1, 2), 1,   zpsi1(1, 1, 0), 1)
    end do
    zpsi(:, :) = M_HALF*(zpsi1(:, :, 0) - zpsi1(:, :, 2))
    call zscal(n, exp(-M_zI*h%spectral_middle_point*timestep), zpsi, 1)
    deallocate(zpsi1)

    if(present(order)) order = td%exp_order
    call pop_sub()
  end subroutine cheby

  subroutine lanczos
    integer ::  korder, is, n, nn, i, info, j
    complex(r8), allocatable :: hm(:, :), v(:, :, :), w(:, :), f(:, :), hh(:), expo(:, :)
    real(r8) :: alpha, beta, res, tol, nrm

    call push_sub('lanczos')

    korder = td%exp_order
    tol = td%lanczos_tol
    allocate(v(sys%m%np, sys%st%dim, korder), &
             w(sys%m%np, sys%st%dim),         &
             f(sys%m%np, sys%st%dim),         &
             hm(korder, korder),         &
             expo(korder, korder),      &
             hh(korder))

    j = 1    
    ! Normalize input vector, and put it into v(:, :, 1)
    nrm = zstates_nrm2(sys%m, sys%st%dim, zpsi(1:sys%m%np, 1:sys%st%dim))
    v(:, :, 1) = zpsi(:, :)/nrm
    ! Operate on v(:, :, 1) and place it onto w.
    call zhpsi(h, sys%m, sys%st, sys, ik, v(:, 1:sys%st%dim, 1), w(1:sys%m%np, 1:sys%st%dim), t)
    alpha = zstates_dotp(sys%m, sys%st%dim, v(1:sys%m%np, 1:sys%st%dim, 1), w(:, :))
    f(:, :) = w(:, :) - alpha*v(1:sys%m%np, 1:sys%st%dim, 1)
    hm = M_z0; hm(1, 1) = alpha
    beta = zstates_nrm2(sys%m, sys%st%dim, f)
    do n = 1, korder - 1
       if(h%ab.ne.1) j = n
       v(1:sys%m%np, 1:sys%st%dim, n + 1) = f(1:sys%m%np, 1:sys%st%dim)/beta
       hm(n+1, n) = beta
       call zhpsi(h, sys%m, sys%st, sys, ik, &
            v(:, 1:sys%st%dim, n+1), w(1:sys%m%np, 1:sys%st%dim), t)
       hh = M_z0
       do nn = j, n + 1
          hh(nn) = zstates_dotp(sys%m, sys%st%dim, v(1:, :, nn), w)
       enddo
       f = w
       do nn = j, n + 1
          f = f - v(1:, :, nn)*hh(nn)
       enddo
       do nn = j, n  + 1
          hm(nn, n + 1) = hh(nn)
       enddo
       if(h%ab .ne. 1) then
         call zmatexp(n+1, hm(1:n+1, 1:n+1), expo(1:n+1, 1:n+1), -M_zI*timestep, method = 2)
       else
         call zmatexp(n+1, hm(1:n+1, 1:n+1), expo(1:n+1, 1:n+1), -M_zI*timestep, method = 1)
       endif
       res = abs(beta*abs(expo(1, n+1)))
       beta = zstates_nrm2(sys%m, sys%st%dim, f)
       if(beta < 1.0e-12_r8) exit
       if(n>1 .and. res<tol) exit
    enddo
    korder = n + 1
    if(res > tol .and. beta > 1.0e-12_r8) then
      write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
      call write_warning(1)
    endif

    zpsi = M_z0
    do nn = 1, korder
       call zaxpy(sys%m%np*sys%st%dim, nrm*expo(nn, 1), v(1, 1, nn), 1, zpsi, 1)
    enddo

    if(present(order)) order = korder
    deallocate(v, w, f, hm, expo, hh)
    call pop_sub()
  end subroutine lanczos

#ifdef HAVE_FFT
  subroutine split
    call push_sub('split')

    if(h%gauge == 2) then
      message(1) = 'Split operator does not work well if velocity gauge is used.'
      call write_fatal(1)
    endif

    call zexp_vlpsi (h, sys%m, sys%st, ik, zpsi, t, -M_zI*timestep/M_TWO)
    if(sys%nlpp) call zexp_vnlpsi (ik, sys%m, sys%st, sys, zpsi, -M_zI*timestep/M_TWO, .true.)
    call zexp_kinetic(h, ik, sys%m, sys%st, zpsi, -M_zI*timestep)
    if(sys%nlpp) call zexp_vnlpsi (ik, sys%m, sys%st, sys, zpsi, -M_zI*timestep/M_TWO, .false.)
    call zexp_vlpsi (h, sys%m, sys%st, ik, zpsi, t, -M_zI*timestep/M_TWO)

    if(present(order)) order = 0
    call pop_sub()
  end subroutine split

  subroutine suzuki
    real(r8) :: tim(5), tt, dt(5)
    real(r8) :: p, pp(5)
    integer :: ist, k

    call push_sub('suzuki')

    if(h%gauge == 2) then
      message(1) = 'Suzuki-Trotter operator does not work well if velocity gauge is used.'
      call write_fatal(1)
    endif

    p = M_ONE/(M_FOUR - M_FOUR**(M_THIRD))
    pp = (/ p, p, M_ONE-M_FOUR*p, p, p /)
    dt(1:5) = pp(1:5)*timestep

    do k = 1, 5
       call zexp_vlpsi (h, sys%m, sys%st, ik, zpsi, t, -M_zI*dt(k)/M_TWO)
       if(sys%nlpp) call zexp_vnlpsi (ik, sys%m, sys%st, sys, zpsi, -M_zI*dt(k)/M_TWO, .true.)
       call zexp_kinetic(h, ik, sys%m, sys%st, zpsi, -M_zI*dt(k))
       if(sys%nlpp) call zexp_vnlpsi (ik, sys%m, sys%st, sys, zpsi, -M_zI*dt(k)/M_TWO, .false.)
       call zexp_vlpsi (h, sys%m, sys%st, ik, zpsi, t, -M_zI*dt(k)/M_TWO)
    end do

    if(present(order)) order = 0
    call pop_sub()
  end subroutine suzuki
#endif

end subroutine td_dtexp
