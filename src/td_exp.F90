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

subroutine td_dtexp(h, sys, td, ik, zpsi, timestep, t)
  type(hamiltonian_type), intent(in) :: h
  type(system_type), intent(in)      :: sys
  type(td_type), intent(in)          :: td
  integer, intent(in) :: ik
  complex(r8), intent(inout) :: zpsi(0:sys%m%np, sys%st%dim)
  real(r8), intent(in) :: timestep, t

  sub_name = 'td_dtexp'; call push_sub()

  select case(td%exp_method)
   case(FOURTH_ORDER);      call fourth
   case(LANCZOS_EXPANSION); call lanczos
   case(SPLIT_OPERATOR);    call split
   case(SUZUKI_TROTTER);    call suzuki
   case(CHEBYSHEV);         call cheby
  end select

  call pop_sub(); return
contains

  subroutine fourth
    integer :: order
    complex(r8) :: zfact
    complex(r8), allocatable :: zpsi1(:,:), hzpsi1(:,:)
    integer :: i

    sub_name = 'fourth'; call push_sub()

    order = td%exp_order
    allocate(zpsi1(0:sys%m%np, sys%st%dim), hzpsi1(sys%m%np, sys%st%dim))
    zfact = 1._r8
    zpsi1 = zpsi
    do i = 1, order
      zfact = zfact*(-M_zI*timestep)/i
      call zHpsi(h, sys%m, sys%st, sys, ik, zpsi1, hzpsi1, t)
      zpsi(1:,:) = zpsi(1:,:) + zfact*hzpsi1(:,:)
      if(i .ne. order) zpsi1(1:,:) = hzpsi1(1:,:)
    end do
    deallocate(zpsi1, hzpsi1)

    call pop_sub(); return
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
    integer :: order = 4, j
    complex(r8) :: zfact
    complex(r8), allocatable :: zpsi1(:,:,:)
    sub_name = 'cheby'; call push_sub()
    
    order = td%exp_order
    allocate(zpsi1(0:sys%m%np, sys%st%dim, 0:2))
    zpsi1 = M_z0
    do j = order, 0, -1
       call zcopy((sys%m%np+1)*sys%st%dim, zpsi1(0, 1, 1), 1, zpsi1(0, 1, 2), 1)
       call zcopy((sys%m%np+1)*sys%st%dim, zpsi1(0, 1, 0), 1, zpsi1(0, 1, 1), 1)
       call zhpsi(h, sys%m, sys%st, sys, ik, zpsi1(:, :, 1), zpsi1(1:, :, 0), t)
            zfact = 2*(-M_zI)**j*oct_bessel(j, h%spectral_half_span*timestep)
       call zaxpy((sys%m%np+1)*sys%st%dim, cmplx(-h%spectral_middle_point, 0.0_r8, r8), &
                                                                     zpsi1(0, 1, 1), 1, zpsi1(0, 1, 0), 1)
       call zscal((sys%m%np+1)*sys%st%dim, cmplx(1./h%spectral_half_span,0._r8, r8),    zpsi1(0, 1, 0), 1)
       call zscal((sys%m%np+1)*sys%st%dim, cmplx(2._r8,0._r8, r8),                      zpsi1(0, 1, 0), 1)
       call zaxpy((sys%m%np+1)*sys%st%dim, zfact, zpsi(0, 1),                      1,   zpsi1(0, 1, 0), 1)
       call zaxpy((sys%m%np+1)*sys%st%dim, cmplx(-1._r8,0._r8,r8), zpsi1(0, 1, 2), 1,   zpsi1(0, 1, 0), 1)
    enddo
    zpsi(:, :) = 0.5_r8*(zpsi1(:, :, 0) - zpsi1(:, :, 2))
    call zscal((sys%m%np+1)*sys%st%dim, exp(-M_zI*h%spectral_middle_point*timestep), zpsi, 1)
    deallocate(zpsi1)

    call pop_sub(); return
  end subroutine cheby

  subroutine lanczos
    integer ::  korder, is, n, nn, i, info, order, j
    complex(r8), allocatable :: hm(:, :), v(:, :, :), w(:, :), f(:, :), hh(:), expo(:, :)
    real(r8) :: alpha, beta, res, tol, nrm

    sub_name = 'lanczos'; call push_sub()

    korder = td%exp_order
    tol = td%lanczos_tol
    allocate(v(0:sys%m%np, sys%st%dim, korder), &
             w(  sys%m%np, sys%st%dim),         &
             f(  sys%m%np, sys%st%dim),         &
             hm(korder, korder),         &
             expo(korder, korder),      &
             hh(korder))
    v(0, :, :) = M_z0

    j = 1    
    ! Normalize input vector, and put it into v(:, :, 1)
    nrm = zstates_nrm2(sys%m, sys%st%dim, zpsi(1:sys%m%np, 1:sys%st%dim))
    v(:, :, 1) = zpsi(:, :)/nrm
    ! Operate on v(:, :, 1) and place it onto w.
    call zhpsi(h, sys%m, sys%st, sys, ik, v(0:sys%m%np, 1:sys%st%dim, 1), w(1:sys%m%np, 1:sys%st%dim), t)
    alpha = zstates_dotp(sys%m, sys%st%dim, v(1:sys%m%np, 1:sys%st%dim, 1), w(:, :))
    f(:, :) = w(:, :) - alpha*v(1:sys%m%np, 1:sys%st%dim, 1)
    hm = M_z0; hm(1, 1) = alpha
    beta = zstates_nrm2(sys%m, sys%st%dim, f)
    do n = 1, korder - 1
       if(h%ab.ne.1) j = n
       v(1:sys%m%np, 1:sys%st%dim, n + 1) = f(1:sys%m%np, 1:sys%st%dim)/beta
       hm(n+1, n) = beta
       call zhpsi(h, sys%m, sys%st, sys, ik, &
            v(0:sys%m%np, 1:sys%st%dim, n+1), w(1:sys%m%np, 1:sys%st%dim), t)
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
         call mat_exp(n+1, hm(1:n+1, 1:n+1), expo(1:n+1, 1:n+1), timestep, hermitian = .true.)
       else
         call mat_exp(n+1, hm(1:n+1, 1:n+1), expo(1:n+1, 1:n+1), timestep)
       endif
       res = abs(beta*abs(expo(1, n+1)))
       beta = zstates_nrm2(sys%m, sys%st%dim, f)
       if(beta < 1.0e-12_r8) exit
       if(n>1 .and. res<tol) exit
    enddo
    order = min(korder, n + 1)
    !if(present(lanczos_order)) lanczos_order = order
    if(res > tol .and. beta > 1.0e-12_r8) then
      write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
      call write_warning(1)
    endif

    zpsi = M_z0
    do nn = 1, order
       call zaxpy((sys%m%np+1)*sys%st%dim, nrm*expo(nn, 1), v(0, 1, nn), 1, zpsi, 1)
    enddo

    deallocate(v, w, f, hm, expo, hh)
    call pop_sub(); return
  end subroutine lanczos

  subroutine split
    sub_name = 'split'; call push_sub()

    if(h%gauge == 2) then
      message(1) = 'Split operator does not work well if velocity gauge is used.'
      call write_fatal(1)
    endif

    call local_part(sys%m, t, timestep/M_TWO)
    if(sys%nlpp) call non_local_pp(sys%m, timestep/M_TWO, .true.)
    call kinetic(sys%m, timestep)
    if(sys%nlpp) call non_local_pp(sys%m, timestep/M_TWO, .false.)
    call local_part(sys%m, t, timestep/M_TWO)

    call pop_sub(); return
  end subroutine split

  subroutine suzuki
    real(r8) :: p, tim(5), tt, dt(5), pp(5)
    integer :: ist, k

    sub_name = 'suzuki'; call push_sub()

    if(h%gauge == 2) then
      message(1) = 'Suzuki-Trotter operator does not work well if velocity gauge is used.'
      call write_fatal(1)
    endif

    p = 1.0_r8/(4.0_r8 - 4.0_r8**(1.0_r8/3.0_r8))
    pp(1) = p
    pp(2) = p
    pp(3) = 1 - 4*p
    pp(4) = p
    pp(5) = p
    dt(1:5) = pp(1:5)*timestep

    do k = 1, 5
       call local_part(sys%m, t, dt(k)/M_TWO)
       if(sys%nlpp) call non_local_pp(sys%m, dt(k)/M_TWO, .true.)
       call kinetic(sys%m, dt(k))
       if(sys%nlpp) call non_local_pp(sys%m, dt(k)/M_TWO, .false.)
       call local_part(sys%m, t, dt(k)/M_TWO)
    end do

    call pop_sub(); return
  end subroutine suzuki

  subroutine kinetic(m, dt)
    type(mesh_type), intent(IN) :: m
    real(r8), intent(in) :: dt

    integer :: ix, iy, iz, ixx(3), idim
    real(r8) :: vec, temp(3)
    complex(r8), allocatable :: wf_r(:,:,:), wf_k(:,:,:)

    sub_name = 'kinetic'; call push_sub()

    temp = M_ZERO
    temp(1:conf%dim) = 2.0_r8*M_PI/(sys%m%fft_n(1:conf%dim)*sys%m%h(1:conf%dim))

    allocate(&
         wf_r(sys%m%fft_n(1), sys%m%fft_n(2), sys%m%fft_n(3)), &
         wf_k(sys%m%fft_n(1), sys%m%fft_n(2), sys%m%fft_n(3))  &
         )

    do idim = 1, sys%st%dim

       wf_r = M_z0
       call zmesh_to_cube(sys%m, zpsi(1:,idim), wf_r)
       call fftwnd_f77_one(sys%m%zplanf, wf_r, wf_k)

       do ix = 1, m%fft_n(1)
          ixx(1) = pad_feq(ix, m%fft_n(1), .true.)
          do iy = 1, m%fft_n(2)
             ixx(2) = pad_feq(iy, m%fft_n(2), .true.)
             do iz = 1, m%fft_n(3)
                ixx(3) = pad_feq(iz, m%fft_n(3), .true.)
                vec = sum((temp(:)*ixx(:))**2)/M_TWO
                wf_k(ix, iy, iz) = exp(- M_zI*dt*vec) * wf_k(ix, iy, iz)
             end do
          end do
       end do

       call fftwnd_f77_one(sys%m%zplanb, wf_k, wf_r)
       wf_r = wf_r/real(sys%m%fft_n(1)*sys%m%fft_n(2)*sys%m%fft_n(3), r8)

       zpsi(1:,idim)= M_z0
       call zcube_to_mesh(sys%m, wf_r, zpsi(1:,idim))

    enddo

    deallocate(wf_r, wf_k)
    call pop_sub(); return
  end subroutine kinetic

  subroutine non_local_pp (m, dt, order)
    type(mesh_type), intent(IN) :: m
    real(r8), intent(in) :: dt
    logical, intent(in) :: order

    integer :: is, idim, ia, ikbc, jkbc, l, lm, add_lm, &
               ia_start, ia_end, step, l_start, l_end, kbc_start, kbc_end
    complex(r8) :: uvpsi, p2, ctemp
    complex(r8), allocatable :: lpsi(:), lhpsi(:), initzpsi(:, :)
    type(atom_type), pointer :: atm
    type(specie_type), pointer :: spec
    complex(r8), external :: zdotc

    sub_name = 'vnlpsi'; call push_sub()

    if(order) then
      step = 1;  ia_start = 1; ia_end = sys%natoms
    else
      step = -1; ia_start = sys%natoms; ia_end = 1
    end if

    allocate(initzpsi(0:m%np, 1:sys%st%dim))
    initzpsi = zpsi
    ! Ionic pseudopotential
    do_atm: do ia = ia_start, ia_end, step
      atm => sys%atom(ia)
      spec => atm%spec
      ! do we have a pseudopotential, or a local pot?
      if(spec%local) cycle do_atm

      do_dim: do idim = 1, sys%st%dim
        allocate(lpsi(atm%mps), lHpsi(atm%mps))
        lpsi(:) = initzpsi(atm%jxyz(:), idim)
        lHpsi(:) = M_z0
        if(order) then
          l_start   = 0; l_end   = spec%ps%L_max
          kbc_start = 1; kbc_end = spec%ps%kbc
          add_lm = 1
        else
          l_start   = spec%ps%L_max; l_end = 0
          kbc_start = spec%ps%kbc; kbc_end = 1
          add_lm = (spec%ps%L_max + 1)**2
        end if
        do_l: do l = l_start, l_end, step
          if (l == spec%ps%L_loc) then
            add_lm = add_lm + (2*l + 1)*step
            cycle do_l
          end if

          do_m: do lm = -l*step, l*step, step
            do ikbc = kbc_start, kbc_end, step
              do jkbc = kbc_start, kbc_end, step
                 p2 = zdotc(atm%mps, atm%zuv(:, add_lm, ikbc), 1, atm%zuv(:, add_lm, ikbc), 1)*m%vol_pp
                 ctemp = atm%zuvu(add_lm, ikbc, jkbc)*p2*(-M_zI*dt)
                 uvpsi = zdotc(atm%mps, atm%zuv(:, add_lm, ikbc), 1, lpsi(:), 1) * m%vol_pp* &
                       (exp(ctemp) - M_z1)/p2
                 call zaxpy (atm%mps, uvpsi, atm%zuv(:, add_lm, jkbc), 1, lHpsi(:), 1)
              end do
            end do
            add_lm = add_lm + step
          end do do_m
        end do do_l
        zpsi(atm%jxyz(:), idim) = zpsi(atm%jxyz(:), idim) + lhpsi(:)
        deallocate(lpsi, lHpsi)
      end do do_dim
    end do do_atm

    deallocate(initzpsi)
    call pop_sub(); return
  end subroutine non_local_pp

  subroutine local_part(m, t, dt)
    type(mesh_type), intent(IN) :: m
    real(r8), intent(in) :: t, dt

    integer :: j, idim
    real(r8) :: r, x(3), las(3)

    sub_name = 'local_part'; call push_sub()

    if(h%no_lasers > 0 .and. h%gauge == 1) then
      call laser_field(h%no_lasers, h%lasers, t, las)
    end if

    ! Propagation of the local part and external field
    do j = 1, m%np
      r = h%Vpsl(j)
      if(.not.h%ip_app) then
        r = r + h%Vhartree(j)
        select case(sys%st%ispin)
        case(1) ! dim = 1
          r = r + h%Vxc(j, 1)
        case(2) ! dim = 1
          if(modulo(ik, 2) == 0) then ! we have a spin down
            r = r + h%Vxc(j, 1)
          else ! spin down
            r = r + h%Vxc(j, 2)
          end if
        case(3) ! dim = 2
          message(1) = "Suzuki-Trotter propagation not implemented for ispin=3"
          call write_fatal(1)
        end select
      end if
      
      if(h%ab .eq. 1) then
        r = r + h%ab_pot(j)
      end if
      
      if(h%no_lasers > 0 .and. h%gauge == 1) then
        call mesh_xyz(m, j, x)
        r = r + sum(x(1:conf%dim)*las(1:conf%dim))
      end if
      
      ! Warning: this does not work for spinors...
      do idim = 1, sys%st%dim
        zpsi(j, idim) = zpsi(j, idim) * exp(-M_zI*dt*r)
      enddo
    end do

    call pop_sub(); return
  end subroutine local_part


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! /* Calculates exp(-M_zI*dt*in), and places it into out, where in and out
  ! are order x order complex matrices. Optional argument hermitian should be
  ! passed if in is hermitian.
  !
  ! Two ways to calculate the exponential of a (small) matrix: if it is
  ! hermitian, by decomposing it into UDU(+), where D is diagonal and
  ! U is unitary; U(+) is its hermitian conjugate. If it is not, I use the 
  ! definition: power expansion (I think 10th order should be enough). Smarter 
  ! methods could be used in this case: a Sch\"{u}r decomposition could be 
  ! performed, (through LAPACK's zhseqr since input matrix is upper Heisenberg) 
  ! and then apply Parlett's algorithm (B. N. Parlett, Lin. Alg. Appl. 14,
  ! 117 (1976)) to calculate the exponential of a upper triangular matrix.
  ! Probably this is not needed. */
  subroutine mat_exp(order, in, out, dt, hermitian)
    implicit none
    integer, intent(in)           :: order
    complex(r8), intent(in)       :: in(order, order)
    complex(r8), intent(out)      :: out(order, order)
    real(r8), intent(in)          :: dt
    logical, intent(in), optional :: hermitian

    complex(r8) :: aux(order, order), dd(order, order), zfact
    integer ::  n, info, lwork
    real(r8), allocatable :: w(:)
    complex(r8), allocatable :: work(:), rwork(:)

    if(present(hermitian)) then
      aux = in
      lwork = 4*order; allocate(work(lwork), rwork(lwork), w(order))
      call zheev('v', 'u', order, aux, order, w, work, lwork, rwork, info)
      if(info .ne. 0) then
         write(message(1),'(a,i4)') 'Error: "info" parameter of z returned', info
         call write_fatal(1)
      endif
      dd = M_z0
      do n = 1, order
         dd(n, n) = exp(-M_zI*dt*w(n))
      enddo
      out = matmul(dd, transpose(aux))
      out = matmul(aux, out)
      deallocate(work, w, rwork)
    else
      out = M_z0
      do n = 1, order
         out(n, n) = M_z1
      enddo
      zfact = M_z1
      aux   = out
      do n = 1, max(10, order)
         aux = matmul(aux, in)
         zfact = zfact*(-M_zI*dt)/n
         out = out + zfact*aux
      enddo
    endif

  end subroutine mat_exp

end subroutine td_dtexp
