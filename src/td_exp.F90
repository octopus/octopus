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
   case(FOURTH_ORDER);      call fourth_order
   case(LANCZOS_EXPANSION); call lanczos_expansion
   case(SPLIT_OPERATOR);    call split_operator
  end select

  call pop_sub(); return
contains

  subroutine fourth_order
    integer, parameter :: order = 4
    complex(r8) :: zfact
    complex(r8), allocatable :: zpsi1(:,:), hzpsi1(:,:)
    integer i, k, idim

    allocate(zpsi1(0:sys%m%np, sys%st%dim), hzpsi1(sys%m%np, sys%st%dim))
    zfact = 1._r8
    zpsi1 = zpsi
    do i = 1, order ! forth order method
      zfact = zfact*(-M_zI*timestep)/i
      call zHpsi(h, sys, ik, zpsi1, hzpsi1, t)
      zpsi(1:,:) = zpsi(1:,:) + zfact*hzpsi1(:,:)
      if(i .ne. order) zpsi1(1:,:) = hzpsi1(:,:)
    end do
    deallocate(zpsi1, hzpsi1)

  end subroutine fourth_order

  subroutine lanczos_expansion
    integer ::  korder, is, n, nn, i, info, order
    complex(r8), allocatable :: hm(:, :), v(:, :, :), w(:, :), f(:, :), hh(:), expo(:, :)
    real(r8) :: alpha, beta, res, tol

    sub_name = 'td_lanczos'; call push_sub()

    korder = td%lanczos_max
    tol    = td%lanczos_tol

    allocate(v(0:sys%m%np, sys%st%dim, korder), &
             w(  sys%m%np, sys%st%dim),         &
             f(  sys%m%np, sys%st%dim),         &
             hm(korder, korder),         &
             expo(korder, korder),      &
             hh(korder))

    ! Lanczos Procedure.

    ! Normalize input vector, and put it into v(:, :, 1)
    call zcopy(sys%m%np*sys%st%dim+2, zpsi(:, :), 1, v(:, :, 1), 1)
    beta   = zstates_nrm2(sys%m, sys%st%dim, v(1:sys%m%np, 1:sys%st%dim, 1))
    call zscal(sys%m%np*sys%st%dim+2, cmplx(1/beta, 0.0_r8, r8), v(:, :, 1), 1)
    ! Operate on v(:, :, 1) and place it onto w.
    call zhpsi(h, sys, ik, v(:, :, 1), w(:, :), t)
    alpha = zstates_dotp(sys%m, sys%st%dim, v(1:sys%m%np, 1:sys%st%dim, 1), w(:, :))
    f(:, :) = w(:, :) - alpha*v(1:sys%m%np, 1:sys%st%dim, 1)
    hm(1, 1) = alpha
    do n = 1, korder - 1
       beta = zstates_nrm2(sys%m, sys%st%dim, f)
       v(1:sys%m%np, 1:sys%st%dim, n + 1) = f(:, :)/beta
       ! Check orthogonality
       !call zstates_gram_schmidt(n+1, sys%m, sys%st%dim, v, start = n)
       hm(n+1, n) = beta
       call zhpsi(h, sys, ik, v(:, :, n+1), w(:, :), t)
       do nn = 1, n + 1
          hh(nn) = zstates_dotp(sys%m, sys%st%dim, v(1:sys%m%np, 1:sys%st%dim, nn), w(:, :))
       enddo
       do is = 1, sys%st%dim
          do i = 1, sys%m%np
             f(i, is) = w(i, is) - sum(v(i, is, 1:n+1)*hh(1:n+1))
          enddo
       enddo
       do nn = 1, n + 1
          hm(nn, n+1) = hh(nn)
       enddo
       call mat_exp(n+1, hm(1:n+1, 1:n+1), expo(1:n+1, 1:n+1), timestep)
       res = abs(beta*abs(expo(1, n+1)))
       !write(*, *) res
       if(n>1 .and. res<tol) exit
    enddo
    order = min(korder, n + 1)
    if(res > tol) then
      write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
      call write_warning(1)
    endif

    hh(1:order) = expo(1, 1:order) * &
                  zstates_dotp(sys%m, sys%st%dim, v(1:sys%m%np, 1:sys%st%dim, 1), &
                                                  zpsi(1:sys%m%np, 1:sys%st%dim))
    do is = 1, sys%st%dim
       do i = 1, sys%m%np
          zpsi(i, is) = sum(v(i, is, 1:order)*hh(1:order))
       enddo
    enddo

    deallocate(v, w, f, hm, expo, hh)
    call pop_sub(); return
  end subroutine lanczos_expansion

  subroutine split_operator
    integer :: idim
    real(r8) :: tt

    sub_name = 'split_operator'; call push_sub()

    if(h%gauge == 2) then
      message(1) = 'Split operator does not work well if velocity gauge is used.'
      call write_fatal(1)
    endif

    ! get reference time right
    tt = t - td%dt
    ! propagate with T/2
    call kinetic(sys%m, td%dt/2._r8)
    call non_local_pp(sys%m, td%dt/2._r8, .true.)
    ! propagate local part of the potential
    call local_part(sys%m, tt, td%dt)
    ! propagate with T/2
    call non_local_pp(sys%m, td%dt/2._r8, .false.)
    call kinetic(sys%m, td%dt/2._r8)

    call pop_sub(); return
  end subroutine split_operator

  subroutine kinetic(m, dt)
    type(mesh_type), intent(IN) :: m
    real(r8), intent(in) :: dt

    integer :: ix, iy, iz, ixx(3), idim
    real(r8) :: vec, temp(3)
    complex(r8), allocatable :: wf_r(:,:,:), wf_k(:,:,:)

    sub_name = 'kinetic'; call push_sub()

    temp(:) = 2.0_r8*M_PI/(sys%m%fft_n(:)*sys%m%h(:))

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
          vec = sum((temp(:)*ixx(:))**2)/2._r8

          ! TODO add vector potential
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

  subroutine non_local_pp(m, dt, order)
    type(mesh_type), intent(IN) :: m
    real(r8), intent(in) :: dt
    logical, intent(in) :: order

    integer :: step, ia, ia_start, ia_end, l, l_start, l_end, lm, add_lm, idim
    integer :: ikbc, jkbc, kbc_start, kbc_end
    complex(r8) :: uVpsi, ctemp
    type(atom_type), pointer :: atm
    type(specie_type), pointer :: spec

    sub_name = 'non_local_pp'; call push_sub()

    if(order) then
      step = 1;  ia_start = 1; ia_end = sys%natoms
    else
      step = -1; ia_start = sys%natoms; ia_end = 1
    end if

    do_atm: do ia = ia_start, ia_end, step
      atm => sys%atom(ia)
      spec => atm%spec

      if(spec%local) cycle do_atm

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
          add_lm = add_lm + step*(2*l + 1)
          cycle do_l
        end if

        do_lm: do lm = -l*step, l*step, step
          do ikbc = kbc_start, kbc_end, step
            do jkbc = kbc_start, kbc_end, step
              do idim = 1, sys%st%dim
                 uVpsi = sum(atm%zuV(:, add_lm, ikbc)*zpsi(atm%Jxyz(:), idim))*sys%m%vol_pp
                 ctemp = uVpsi * (exp(-M_zI*dt*atm%zuVu(add_lm, ikbc, jkbc)) - 1.0_r8)
                 zpsi(atm%Jxyz(:), idim) = zpsi(atm%Jxyz(:), idim) + &
                     ctemp * atm%zuV(:, add_lm, jkbc)
              enddo
            end do
          end do

          add_lm = add_lm + step
        end do do_lm
      end do do_l

    end do do_atm

    call pop_sub(); return
  end subroutine non_local_pp

  subroutine local_part(m, t, dt)
    type(mesh_type), intent(IN) :: m
    real(r8), intent(in) :: t, dt

    integer :: j, ik, ix, iy, iz, idim
    real(r8) :: r, x(3), las(3)

    sub_name = 'local_part'; call push_sub()

    if(h%no_lasers > 0 .and. h%gauge == 1) then
      call laser_field(h%no_lasers, h%lasers, t + dt/2._r8, las)
    end if

    ! Propagation of the local part and external field
        do j = 1, m%np
          ix = m%Lx(j) + (m%fft_n(1) - 1)/2 + 1
          iy = m%Ly(j) + (m%fft_n(2) - 1)/2 + 1
          iz = m%Lz(j) + (m%fft_n(3) - 1)/2 + 1

          r = h%Vpsl(j) + h%Vhartree(j)
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

          if(h%ab .eq. 1) then
            r = r + h%ab_pot(j)
          end if

          if(h%no_lasers > 0 .and. h%gauge == 1) then
            call mesh_xyz(m, j, x)
            r = r + sum(x(:)*las(:))
          end if

          ! Warning: this does not work for spinors...
          do idim = 1, sys%st%dim
             zpsi(j, 1) = zpsi(j, 1) * exp(- M_zI*td%dt*r)
          enddo
        end do

    call pop_sub(); return
  end subroutine local_part

  subroutine mat_exp(order, in, out, dt)
    implicit none
    integer, intent(in)      :: order
    complex(r8), intent(in)     :: in(order, order)
    complex(r8), intent(out) :: out(order, order)
    real(r8), intent(in)     :: dt

    complex(r8) :: aux(order, order), b(order, order)
    complex(r8) :: q(order, order), dd(order, order)
    integer ::  n, nn, i, info, lwork
    real(r8), allocatable :: dsygv_w(:)
    complex(r8), allocatable :: work(:), rwork(:)

    aux = in
    b = 0.0_r8
    do n = 1, order
       b(n, n) = 1._r8
    enddo
    lwork = 3*order; allocate(work(lwork), rwork(lwork), dsygv_w(order))
    call zhegv(1, 'v', 'u', order, aux, order, b, order, dsygv_w, work, lwork, rwork, info)
    if(info .ne. 0) then
       write(message(1),'(a,i4)') 'Error: "info" parameter of zpttrf returned', info
       call write_fatal(1)
    endif
    q = aux
    dd = M_z0
    do n = 1, order
       dd(n, n) = exp(-M_zI*dt*dsygv_w(n))
    enddo
    out = matmul(dd, transpose(q))
    out = matmul(q, out)
    deallocate(work, dsygv_w)

  end subroutine mat_exp

end subroutine td_dtexp
