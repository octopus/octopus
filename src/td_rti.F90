subroutine td_rti(sys, h, td, t)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  type(td_type), intent(inout) :: td
  real(r8), intent(in) :: t
  
  sub_name = 'td_rti'; call push_sub()
  
  select case(td%evolution_method)
  case(1)
    call td_rti1(sys%m, sys%st)
  case(2)
    call td_rti2(sys%m, sys%st)
  case(3) ! split operator
!    call td_rti3(sys, h, td, t)
  end select

  call pop_sub()
  return

contains

  subroutine td_dtexp(ik, zpsi, timestep, t)
    integer, intent(in) :: ik
    complex(r8), intent(inout) :: zpsi(0:sys%m%np)
    real(r8), intent(in) :: timestep, t

    integer, parameter :: order = 4

    complex(r8) :: zfact
    complex(r8), allocatable :: zpsi1(:), hzpsi1(:), grad(:,:)
    real(r8) :: x(3), f(3)
    integer i, k, l

    allocate(zpsi1(0:sys%m%np), hzpsi1(sys%m%np))

    zfact = 1._r8
    zpsi1 = zpsi
    do i = 1, order ! forth order method
      zfact = zfact*(-M_zI*timestep)/i
      call zHpsi(h, sys, ik, zpsi1, hzpsi1)
      
      ! apply lasers
      if(td%no_lasers > 0) then
        select case(td%gauge)
        case(1) ! length gauge
          call laser_field(td%no_lasers, td%lasers, t, f)
          
          do k = 1, sys%m%np
            call mesh_xyz(sys%m, k, x)
            hzpsi1(k) = hzpsi1(k) + sum(x*f) * zpsi1(k)
          end do
          
        case(2) ! velocity gauge
          call laser_vector_field(td%no_lasers, td%lasers, t, f)
          allocate(grad(3, sys%m%np))
          call zmesh_derivatives(sys%m, zpsi1, grad=grad)
          do k = 1, sys%m%np
            hzpsi1(k) = hzpsi1(k) - M_zI * sum(f(:)*grad(:, k)) + &
                 sum(f**2)/2._r8 * zpsi1(k)
          end do
          deallocate(grad)
        end select
      end if
      
      ! absorbing potential
      if(td%ab .eq. 1) then
        hzpsi1 = hzpsi1 + M_zI*td%ab_pot(:)*zpsi1(1:)
      end if
      
      zpsi(1:) = zpsi(1:) + zfact*hzpsi1(:)
      
      if(i .ne. order) zpsi1(1:) = hzpsi1(:)
    end do
    
    deallocate(zpsi1, hzpsi1)
    return
  end subroutine td_dtexp
  
  ! Warning: this subroutine should only be used with LDA/GGA functionals
  subroutine td_rti1(m, st)
    type(mesh_type), intent(IN) :: m
    type(states_type), intent(inout) :: st

    integer is, ik, ist
    real(r8), allocatable :: aux(:,:)
    complex(r8), allocatable :: zpsi1(:,:,:,:)
    
    allocate(aux(m%np, st%nspin))
    
    do is = 1, sys%st%nspin
      aux(:, is) = 1.875_r8*(h%VHartree(:) + h%Vxc(:, is)) &
           -1.25_r8*td%v_old1(:, is) + 0.375_r8*td%v_old2(:, is)
    end do
    
    td%v_old2 = td%v_old1
    do is = 1, st%nspin
      td%v_old1(:, is) = h%VHartree(:) + h%Vxc(:, is)
    end do
    
    h%VHartree = 0._r8; h%Vxc = aux
    allocate(zpsi1(0:m%np, st%dim, st%st_start:st%st_end, st%nik))
    zpsi1 = st%zpsi
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:, :, ist, ik), td%dt, t-td%dt)
      end do
    end do
    st%zpsi = zpsi1
    deallocate(zpsi1)
    
    call zcalcdens(st, m%np, aux, .true.)
    st%rho = (st%rho +  aux) / 2.0_r8
    deallocate(aux)
    
    call hamiltonian_setup(h, sys)
    
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, sys%st%zpsi(:, :, ist, ik), td%dt, t-td%dt)
      end do
    end do
    
    return
  end subroutine td_rti1

  subroutine td_rti2(m, st)
    type(mesh_type), intent(inout) :: m
    type(states_type), intent(inout) :: st
    
    real(r8), allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
    complex(r8), allocatable :: zpsi1(:,:,:,:)
    integer is, ik, ist
#if defined(HAVE_MPI) && defined(MPI_TD)
    real(r8), allocatable :: reduce_rho(:,:) ! temporary to do MPI_reduce
    integer :: ierr
#endif

    allocate(zpsi1(0:m%np, st%dim, st%st_start:st%st_end, st%nik))
    zpsi1 = st%zpsi ! store zpsi
    
    allocate(vhxc_t1(m%np, st%nspin))
    do is = 1, st%nspin ! store Vhxc
      Vhxc_t1(:, is) = h%VHartree(:) + h%Vxc(:, is)
    end do
    
    ! propagate dt with H(t-dt)
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:,:, ist, ik), td%dt, t-td%dt)
      end do
    end do
    
    call zcalcdens(st, m%np, st%rho, .true.)
    call hamiltonian_setup(h, sys)
    
    st%zpsi = zpsi1
    deallocate(zpsi1)
    
    ! store Vhxc at t
    allocate(vhxc_t2(m%np, st%nspin))
    do is = 1, st%nspin
      Vhxc_t2(:, is) = h%VHartree(:) + h%Vxc(:, is)
    end do
    
    ! propagate dt/2 with H(t-dt)
    h%Vhartree = 0._r8
    h%Vxc = Vhxc_t1
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t-td%dt)
      end do
    end do
    
    deallocate(vhxc_t1)
    
    ! propagate dt/2 with H(t)
    h%Vxc = Vhxc_t2
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t)
      end do
    end do
    
    deallocate(vhxc_t2)
    
    return
  end subroutine td_rti2

end subroutine td_rti

! TODO
!!$subroutine systm_td_rti3(sys, td, t)
!!$  type(sys_type), intent(inout) :: sys
!!$  type(td_type), intent(inout) :: td
!!$  real(r8), intent(in) :: t
!!$
!!$  integer :: ia, is, p, j, ix, iy, iz, ixx, iyy, izz, l, lm
!!$  real(r8) :: f(3), temp, vec
!!$  complex(r8) :: ctemp
!!$  complex(r8), allocatable :: wf(:,:,:)
!!$
!!$  if(td%no_lasers > 0) then 
!!$    select case(td%gauge)
!!$    case(1)
!!$      call laser_field(td%no_lasers, td%lasers, t, f)
!!$    case(2)
!!$      call laser_vector_field(td%no_lasers, td%lasers, t, f)
!!$    end select
!!$  end if
!!$
!!$  allocate(wf(td%fft_n, td%fft_n, td%fft_n))
!!$  wf = 0._r8
!!$  temp = 2.0_r8*M_PI/(td%fft_n*sys%m%h)
!!$
!!$  do is = 1, sys%nspin
!!$    do p = sys%st_start, sys%st_end
!!$      ! setup wf in cubic grid
!!$      do j = 1, sys%m%NL
!!$        ix = sys%m%Lx(j) + (td%fft_n-1)/2 + 1
!!$        iy = sys%m%Ly(j) + (td%fft_n-1)/2 + 1
!!$        iz = sys%m%Lz(j) + (td%fft_n-1)/2 + 1
!!$        
!!$        wf(ix, iy, iz) = td%zpsi(j, p, is)
!!$      end do
!!$
!!$      ! Fourier transform
!!$      call fftwnd_f77_one(td%planf, wf, 0)
!!$      
!!$      ! propagate with T/2
!!$      wf = wf * td%kin_2
!!$      if(td%no_lasers > 0 .and. td%gauge == 2) then
!!$        do ix = 1, td%fft_n
!!$          do iy = 1, td%fft_n
!!$            do iz = 1, td%fft_n
!!$              ixx = pad_feq(ix, td%fft_n, .true.)
!!$              iyy = pad_feq(iy, td%fft_n, .true.)
!!$              izz = pad_feq(iz, td%fft_n, .true.)
!!$              
!!$              wf(ix, iy, iz) = wf(ix, iy, iz) * exp(- M_zI*P_H2M*td%dt/2._r8*( &
!!$                   2._r8*temp*(f(1)*real(ixx,r8) + f(2)*real(iyy,r8) + f(3)*real(izz,r8)) &
!!$                   + sum(f**2)))
!!$            end do
!!$          end do
!!$        end do
!!$      end if
!!$
!!$      ! Fourier transform back
!!$      call fftwnd_f77_one(td%planb, wf, 0)
!!$      wf = wf/real(td%fft_n**3, r8)
!!$
!!$      ! now propagate with non-local part of the pseudopotentials
!!$      do ia = 1, sys%nions 
!!$        do l = 0 , sys%spec(sys%ion(ia)%spec)%ps_lmax 
!!$          if (l /= sys%spec(sys%ion(ia)%spec)%ps_lloc) then
!!$            do lm = 1, (sys%spec(sys%ion(ia)%spec)%ps_Lmax+1)**2
!!$              ctemp = (0.0_r8,0.0_r8)
!!$              do j = 1, sys%ion(ia)%Mps
!!$                ix = sys%m%Lx(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                iy = sys%m%Ly(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                iz = sys%m%Lz(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                ctemp = ctemp +  sys%ion(ia)%uV(j, lm)*wf(ix, iy, iz)
!!$              end do
!!$              ctemp = ctemp * sys%m%vol_pp * (exp(-M_zI*td%dt/2._r8*              &
!!$                   sys%spec(sys%ion(ia)%spec)%ps_kb(-1,l)/             &
!!$                   sys%spec(sys%ion(ia)%spec)%ps_kb(-2,l)**2)-1.0_r8)* &
!!$                   sys%spec(sys%ion(ia)%spec)%ps_kb(-2,l)**2 
!!$              
!!$              do j = 1, sys%ion(ia)%Mps
!!$                ix = sys%m%Lx(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                iy = sys%m%Ly(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                iz = sys%m%Lz(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                wf(ix, iy, iz) = wf(ix, iy, iz) +  sys%ion(ia)%uV(j, lm) * ctemp
!!$              end do
!!$              
!!$            end do
!!$          end if
!!$        end do
!!$      end do
!!$
!!$      !Propagation of the local part and external field
!!$      do j = 1, sys%m%NL
!!$        ix = sys%m%Lx(j) + (td%fft_n-1)/2 + 1
!!$        iy = sys%m%Ly(j) + (td%fft_n-1)/2 + 1
!!$        iz = sys%m%Lz(j) + (td%fft_n-1)/2 + 1
!!$        
!!$        wf(ix, iy, iz) = wf(ix, iy, iz) * exp( - M_zI*td%dt*( &
!!$             sys%Vpsl(j) + sys%Vh(j) + sys%Vx(j, is) + sys%Vc(j, is)))
!!$
!!$        if(td%ab .eq. 1) then
!!$          wf(ix, iy, iz) = wf(ix, iy, iz) * exp( td%dt*td%ab_pot(j) )               
!!$        end if
!!$      end do
!!$
!!$      if(td%no_lasers > 0 .and. td%gauge == 1) then
!!$        do j = 1, sys%m%nl
!!$          ix = sys%m%Lx(j) + (td%fft_n-1)/2 + 1
!!$          iy = sys%m%Ly(j) + (td%fft_n-1)/2 + 1
!!$          iz = sys%m%Lz(j) + (td%fft_n-1)/2 + 1
!!$          
!!$          wf(ix, iy, iz) = wf(ix, iy, iz) * exp( - M_zI*td%dt* &
!!$               sys%m%h*(sys%m%Lx(j)*f(1) + sys%m%Ly(j)*f(2) + sys%m%Lz(j)*f(3)))       
!!$        end do
!!$      end if
!!$
!!$      ! now propagate (inverse-order) non-local part of the pseudopotentials
!!$      do ia = sys%nions , 1 , -1
!!$        do l =  sys%spec(sys%ion(ia)%spec)%ps_lmax , 0 , -1
!!$          if (l /= sys%spec(sys%ion(ia)%spec)%ps_lloc) then
!!$            do lm = (sys%spec(sys%ion(ia)%spec)%ps_Lmax+1)**2, 1, -1
!!$              ctemp = (0.0_r8,0.0_r8)
!!$              do j = 1, sys%ion(ia)%Mps
!!$                ix = sys%m%Lx(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                iy = sys%m%Ly(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                iz = sys%m%Lz(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                ctemp = ctemp +  sys%ion(ia)%uV(j, lm)*wf(ix, iy, iz)
!!$              end do
!!$              ctemp = ctemp * sys%m%vol_pp * (exp(-M_zI*td%dt/2._r8*              &
!!$                   sys%spec(sys%ion(ia)%spec)%ps_kb(-1,l)/             &
!!$                   sys%spec(sys%ion(ia)%spec)%ps_kb(-2,l)**2)-1.0_r8)* &
!!$                   sys%spec(sys%ion(ia)%spec)%ps_kb(-2,l)**2 
!!$              
!!$              do j = 1, sys%ion(ia)%Mps
!!$                ix = sys%m%Lx(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                iy = sys%m%Ly(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                iz = sys%m%Lz(sys%ion(ia)%Jxyz(j)) + (td%fft_n-1)/2 + 1
!!$                wf(ix, iy, iz) = wf(ix, iy, iz) +  sys%ion(ia)%uV(j, lm) * ctemp
!!$              end do
!!$              
!!$            end do
!!$          end if
!!$        end do
!!$      end do
!!$
!!$      call fftwnd_f77_one(td%planf, wf, 0)
!!$
!!$      ! propagate with T/2
!!$      wf = wf * td%kin_2
!!$      if(td%no_lasers > 0 .and. td%gauge == 2) then
!!$        do ix = 1, td%fft_n
!!$          do iy = 1, td%fft_n
!!$            do iz = 1, td%fft_n
!!$              ixx = pad_feq(ix, td%fft_n, .true.)
!!$              iyy = pad_feq(iy, td%fft_n, .true.)
!!$              izz = pad_feq(iz, td%fft_n, .true.)
!!$              
!!$              wf(ix, iy, iz) = wf(ix, iy, iz) * exp(- M_zI*P_H2M*td%dt/2._r8*( &
!!$                   2._r8*temp*(f(1)*real(ixx,r8) + f(2)*real(iyy,r8) + f(3)*real(izz,r8)) &
!!$                   + sum(f**2)))
!!$            end do
!!$          end do
!!$        end do
!!$      end if
!!$
!!$      ! Fourier transform back
!!$      call fftwnd_f77_one(td%planb, wf, 0)
!!$      wf = wf/real(td%fft_n**3, r8)
!!$
!!$      do j = 1, sys%m%NL
!!$        ix = sys%m%Lx(j) + (td%fft_n-1)/2 + 1
!!$        iy = sys%m%Ly(j) + (td%fft_n-1)/2 + 1
!!$        iz = sys%m%Lz(j) + (td%fft_n-1)/2 + 1
!!$
!!$        td%zpsi(j, p, is) = wf(ix, iy, iz)
!!$      end do
!!$
!!$    end do
!!$  end do
!!$
!!$  deallocate(wf)
!!$end subroutine systm_td_rti3
