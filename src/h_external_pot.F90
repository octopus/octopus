! This subroutine should be in specie.F90, but due to the limitations
! of f90 to handle circular dependences it had to come here!

subroutine specie_local_fourier_init(ns, s, m)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)
  type(mesh_type), intent(IN) :: m

  integer :: i, j, ix, iy, iz, n, ixx(3)
  real(r8) :: x(3), r, vl, r1
  complex(r8) :: c
  real(r8), allocatable :: fr(:,:,:), dfr(:,:,:,:)

  sub_name = 'specie_local_fourier_init'; call push_sub()

  allocate(fr(m%fft_n2(1),  m%fft_n2(2), m%fft_n2(3)))

  do i = 1, ns
    allocate( &
         s(i)%local_fw(m%hfft_n2, m%fft_n2(2), m%fft_n2(3)), &
         s(i)%rhocore_fw(m%hfft_n2, m%fft_n2(2), m%fft_n2(3)))

    do ix = 1, m%fft_n2(1)
      ixx(1) = ix - (m%fft_n2(1)/2 + 1)
      do iy = 1, m%fft_n2(2)
        ixx(2) = iy - (m%fft_n2(2)/2 + 1)
        do iz = 1, m%fft_n2(3)
          ixx(3) = iz - (m%fft_n2(3)/2 + 1)
          r = sqrt(sum((m%h(:)*ixx(:))**2))
          vl  = splint(s(i)%ps%vlocal,  r)
          if(r >= r_small) then
            fr(ix, iy, iz) = (vl - s(i)%Z_val)/r
          else
            fr(ix, iy, iz) = s(i)%ps%vlocal_origin
          endif
        end do
      end do
    end do

    call rfftwnd_f77_one_real_to_complex(m%dplanf2, fr, s(i)%local_fw)

    n = m%fft_n2(2)*m%fft_n2(3)*m%hfft_n2
    c  = cmplx(1.0_r8/(m%fft_n2(1)*m%fft_n2(2)*m%fft_n2(3)), 0.0_r8, r8)
    call zscal(n, c, s(i)%local_fw,  1)

    if(s(i)%ps%icore /= 'nc  ') then
      fr = 0.0_r8
      do ix = 1, m%fft_n2(1)
        do iy = 1, m%fft_n2(2)
          do iz = 1, m%fft_n2(3)
             r = sqrt( ( m%h(1)*(ix - m%fft_n2(1)/2 - 1) )**2 + &
                       ( m%h(2)*(iy - m%fft_n2(2)/2 - 1) )**2 + &
                       ( m%h(3)*(iz - m%fft_n2(3)/2 - 1) )**2 )
            vl  = splint(s(i)%ps%core, r)
            fr(ix, iy, iz) = vl
          end do
        end do
      end do
      call rfftwnd_f77_one_real_to_complex(m%dplanf2, fr, s(i)%rhocore_fw)
      call zscal(n, c, s(i)%rhocore_fw, 1)
    end if
  end do
  
  deallocate(fr)
  call pop_sub()
end subroutine specie_local_fourier_init

subroutine specie_nl_fourier_init(ns, s, m, nextra)
  integer, intent(in) :: ns, nextra
  type(specie_type), pointer :: s(:)
  type(mesh_type), intent(IN) :: m

  integer :: n(3), hn, i, j, ix, iy, iz, ixx, iyy, izz, l, lm, add_lm, kx, ky, kz
  integer(POINTER_SIZE) :: nl_planf
  real(r8) :: r, x(3), vl, g(3)
  complex(r8) :: c
  real(r8), allocatable :: fr(:,:,:), dfr(:,:,:,:)
  complex(r8), allocatable :: fw(:,:,:), dfw(:,:,:,:)

  sub_name = 'specie_nl_fourier_init'; call push_sub()
  
  specie_loop: do i = 1, ns
    ! first get the dimensions of the thing
    n(1:3) = 2*nint(s(i)%ps%rc_max/m%h(1:3)) + 5 ! this should be enough
    hn = n(1)/2 + 1
    s(i)%nl_fft_n(1:3) = n(1:3)
    s(i)%nl_hfft_n = hn

    ! allocate memory and FFT plans
    allocate(s(i)%nl_fw(hn, n(2), n(3), (s(i)%ps%L_max + 1)**2), &
        s(i)%nl_dfw(hn, n(2), n(3), 3, (s(i)%ps%L_max + 1)**2))
    call rfftw3d_f77_create_plan(s(i)%nl_planb, n(1), n(2), n(3), &
         fftw_backward, fftw_measure + fftw_threadsafe)    

    n(1:3) = (n(1:3) - 1)*(1 + nextra) + 1
    hn = n(1)/2 + 1

    call rfftw3d_f77_create_plan(nl_planf, n(1), n(2), n(3), &
         fftw_forward, fftw_measure + fftw_threadsafe)

    ! fill in structure
    allocate(fr(n(1), n(2), n(3)), fw(hn, n(2), n(3)), &
         dfr(n(1), n(2), n(3), 3), dfw(hn, n(2), n(3), 3))

    ! we will recalculate this value on the mesh
    add_lm = 1
    l_loop: do l = 0, s(i)%ps%L_max
      if(l == s(i)%ps%L_loc) then
        add_lm = add_lm + (2*l + 1)
        cycle l_loop
      end if
      
      lm_loop: do lm = -l, l
        do ix = 1, n(1)
          x(1) = m%h(1) * real(ix - n(1)/2 - 1, r8) / real(1 + nextra, r8)
          do iy = 1, n(2)
            x(2) = m%h(2) * real(iy - n(2)/2 - 1, r8) / real(1 + nextra, r8)
            do iz = 1, n(3)
              x(3) = m%h(3) * real(iz - n(3)/2 - 1, r8) / real(1 + nextra, r8)
              r = sqrt(sum(x**2))
              
              call get_nl_part(s(i)%ps, x, l, lm, fr(ix, iy, iz), dfr(ix, iy, iz, :))
            end do
          end do
        end do
        c = cmplx(1._r8/real(n(1)*n(2)*n(3)), 0.0_r8, r8)
        call rfftwnd_f77_one_real_to_complex(nl_planf, fr, fw)
        call zscal(hn*n(2)*n(3), c, fw(1, 1, 1), 1)

        do j = 1, 3
          call rfftwnd_f77_one_real_to_complex(nl_planf, dfr(1, 1, 1, j), dfw(1, 1, 1, j))
          call zscal(hn*n(2)*n(3), c, dfw(1, 1, 1, j), 1)
        end do

        ! have to cut the high frequencies
        do ix = 1, s(i)%nl_hfft_n
          kx  = pad_feq(ix, s(i)%nl_fft_n(1), .true.)
          ixx = pad_feq(kx, n(1), .false.)
          do iy = 1, s(i)%nl_fft_n(2)
            ky  = pad_feq(iy, s(i)%nl_fft_n(2), .true.)
            iyy = pad_feq(ky, n(2), .false.)
            do iz = 1, s(i)%nl_fft_n(3)
              kz  = pad_feq(iz, s(i)%nl_fft_n(3), .true.)
              izz = pad_feq(kz, n(3), .false.)
              c =  exp(-M_PI*M_zi* (&
                   kx*(1._r8/n(1)-1._r8/s(i)%nl_fft_n(1)) + &
                   ky*(1._r8/n(2)-1._r8/s(i)%nl_fft_n(2)) + &
                   kz*(1._r8/n(3)-1._r8/s(i)%nl_fft_n(3))))
              s(i)%nl_fw (ix, iy, iz, add_lm)   = c * fw (ixx, iyy, izz)
              s(i)%nl_dfw(ix, iy, iz,:, add_lm) = c * dfw(ixx, iyy, izz, :)
            end do
          end do
        end do
        
        add_lm = add_lm + 1
      end do lm_loop
    end do l_loop

    call fftw_f77_destroy_plan(nl_planf)
    deallocate(fr, fw, dfr, dfw)
  end do specie_loop

  call pop_sub()
end subroutine specie_nl_fourier_init

subroutine generate_external_pot(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys

  integer :: ia
  type(specie_type), pointer :: s
  type(atom_type),   pointer :: a
  complex(r8), allocatable :: fw(:,:,:), fwc(:,:,:)
  real(r8), allocatable :: fr(:,:,:)

  ! WARNING DEBUG
!!$  integer :: j
!!$  real(r8), allocatable :: f(:,:,:)

  integer :: i

  sub_name = 'generate_external_pot'; call push_sub()

  ! first we assume that we need to recalculate the ion_ion energy
  call ion_ion_energy(sys)

  if(h%vpsl_space == 1) then
    allocate(fw(sys%m%hfft_n2, sys%m%fft_n2(2), sys%m%fft_n2(3)), &
         fwc(sys%m%hfft_n2, sys%m%fft_n2(2), sys%m%fft_n2(3)))
    fw = M_z0; fwc = M_z0
  end if

  h%Vpsl  = 0._r8
  do ia = 1, sys%natoms
    a => sys%atom(ia) ! shortcuts
    s => a%spec

    select case(s%label(1:5))
    case('jelli', 'point')
      call from_jellium(sys%m)
    case default
      call from_pseudopotential(sys%m)
    end select
  end do

  if(h%vpsl_space == 1) then
    allocate(fr(sys%m%fft_n2(1), sys%m%fft_n2(2), sys%m%fft_n2(3)))

    ! first the potential
    call rfftwnd_f77_one_complex_to_real(sys%m%dplanb2, fw, fr)
    call dcube_to_mesh(sys%m, fr, h%Vpsl, t=2)

    ! and the non-local core corrections
    call rfftwnd_f77_one_complex_to_real(sys%m%dplanb2, fwc, fr)
    call dcube_to_mesh(sys%m, fr, h%rho_core, 2)
    
    deallocate(fw, fwc, fr)
  end if

  if (h%classic_pot) then
    h%Vpsl = h%Vpsl + h%Vclassic
  end if

  ! WARNING DEBUG
!!$  allocate(f(sys%m%fft_n, sys%m%fft_n, sys%m%fft_n))
!!$  f = 0._r8
!!$  call dmesh_to_cube(sys%m, h%Vpsl, f)
!!$  do i = 1, sys%m%fft_n
!!$    do j = 1, sys%m%fft_n
!!$      print *, i, j, f(i, j, sys%m%hfft_n)
!!$    end do
!!$    print *
!!$  end do
!!$  stop
  
  call pop_sub()
contains
  !***************************************************
  !  jellium stuff
  !***************************************************
  subroutine from_jellium(m)
    type(mesh_type), intent(in) :: m

    integer :: i
    real(r8) :: a1, a2, Rb2, x(3), r

    a1 = s%Z/(2._r8*s%jradius**3)
    a2 = s%Z/s%jradius
    Rb2= s%jradius**2
    do i = 1, h%np
      call mesh_r(m, i, r, x=x, a=a%x)
      if(r <= s%jradius) then
        h%Vpsl(i)     = h%Vpsl(i) + (a1*(r*r - Rb2) - a2)
      else
        h%Vpsl(i)     = h%Vpsl(i) - s%Z/r
      end if
    end do
    
  end subroutine from_jellium

  !***************************************************
  !  pseudopotential stuff
  !***************************************************
  subroutine from_pseudopotential(m)
    type(mesh_type), intent(in) :: m

    call build_local_part(m)
    call build_kb_sphere(m)
    call build_nl_part(m)

  end subroutine from_pseudopotential

  subroutine build_local_part(m)
    type(mesh_type), intent(in) :: m

    integer :: i
    real(r8) :: r, vl, dvl, x(3)

    if(h%vpsl_space == 0) then
      do i = 1, m%np
        call mesh_r(m, i, r, x=x, a=a%x)
        if(r >= r_small) then
          vl  = splint(s%ps%vlocal,  r)
          dvl = splint(s%ps%dvlocal, r)

          h%Vpsl(i)     = h%Vpsl(i) + (vl - s%Z_val)/r
        else
          h%Vpsl(i) = h%Vpsl(i) + s%ps%Vlocal_origin
        end if
        if(s%ps%icore /= 'nc  ' ) then
          h%rho_core(i) = h%rho_core(i) + splint(s%ps%core, r)
        end if
      end do
    else
      call phase_factor(m, m%fft_n2, a%x, s%local_fw, fw)
      if(s%ps%icore /= 'nc  ') then
        call phase_factor(m, m%fft_n2(1:3), a%x, s%rhocore_fw, fwc)
      end if
    end if

  end subroutine build_local_part
  
  subroutine build_kb_sphere(m)
    type(mesh_type), intent(in) :: m
    
    integer :: j, k
    real(r8) :: r

    ! This is for the ions movement; probably it is not too elegant, I
    ! will rethink it later.
    if(associated(a%jxyz)) deallocate(a%jxyz, a%uv, a%duv, a%uvu)
    
    j = 0
    do k = 1, h%np
      call mesh_r(m, k, r, a=a%x)
      if(r <= s%ps%rc_max + m%h(1)) j = j + 1
    end do
    a%Mps = j
    
    allocate(a%Jxyz(j), a%uV(j, (s%ps%L_max+1)**2), &
         a%duV(3, j, (s%ps%L_max+1)**2), a%uVu((s%ps%L_max+1)**2))
    
    a%uV  = 0.0_r8; a%duV = 0.0_r8; a%uVu = 0.0_r8
    
    j = 0
    do k = 1, h%np
      call mesh_r(m, k, r, a=a%x)
      ! we enlarge slightly the mesh (good for the interpolation scheme)
      if(r <= s%ps%rc_max + m%h(1)) then
        j = j + 1
        a%Jxyz(j) = k
      end if
    end do

  end subroutine build_kb_sphere

  subroutine build_nl_part(m)
    type(mesh_type), intent(IN) :: m

    integer :: j, l, lm, add_lm, p, ix, iy, iz, center(3)
    real(r8) :: r, x(3), ylm
    complex(r8), allocatable :: nl_fw(:,:,:)
    real(r8), allocatable :: nl_fr(:,:,:), nl_dfr(:,:,:,:)

    if(h%vnl_space == 0) then
      j_loop: do j = 1, a%Mps
        call mesh_xyz(m, a%Jxyz(j), x)
        x = x - a%x
        
        add_lm = 1
        l_loop: do l = 0, s%ps%L_max
          lm_loop: do lm = -l , l
            if(l .ne. s%ps%L_loc) then
              call get_nl_part(s%ps, x, l, lm, a%uV(j, add_lm), a%duV(:, j, add_lm))
            end if

            add_lm = add_lm + 1
          end do lm_loop
        end do l_loop
      end do j_loop
    else ! Fourier space
      center(:) = nint(a%x(:)/m%h(:))
      x(:)  = a%x(:) - center(:)*m%h(:)
      
      allocate(nl_fw (s%nl_hfft_n,   s%nl_fft_n(1), s%nl_fft_n(1)), &
               nl_fr (s%nl_fft_n(1), s%nl_fft_n(1), s%nl_fft_n(1)), &
               nl_dfr(s%nl_fft_n(1), s%nl_fft_n(1), s%nl_fft_n(1), 3))
      
      add_lm = 1
      l_loop2: do l = 0, s%ps%L_max
        if(l == s%ps%L_loc) then
          add_lm = add_lm + (2*l + 1)
          cycle l_loop2
        end if
        
        lm_loop2: do lm = -l , l
          nl_fw = M_z0
          call phase_factor(m, s%nl_fft_n(1:3), x, s%nl_fw(:,:,:, add_lm), nl_fw)
          call rfftwnd_f77_one_complex_to_real(s%nl_planb, nl_fw, nl_fr)

          ! now the gradient
          do j = 1, 3
            nl_fw = M_z0
            call phase_factor(m, s%nl_fft_n(1:3), x, s%nl_dfw(:,:,:, j, add_lm), nl_fw)
            call rfftwnd_f77_one_complex_to_real(s%nl_planb, nl_fw, nl_dfr(1, 1, 1, j))
          end do
          
          j_loop2: do j = 1, a%Mps
            p = a%Jxyz(j)
            ix = m%Lx(p) - center(1) + s%nl_fft_n(1)/2 + 1
            iy = m%Ly(p) - center(2) + s%nl_fft_n(2)/2 + 1
            iz = m%Lz(p) - center(3) + s%nl_fft_n(3)/2 + 1
            
            a%uV(j, add_lm)     = nl_fr (ix, iy, iz)
            a%duV(:, j, add_lm) = nl_dfr(ix, iy, iz, :)
          end do j_loop2

          add_lm = add_lm + 1
        end do lm_loop2
      end do l_loop2
      deallocate(nl_fw, nl_fr, nl_dfr)
    end if
    
    ! and here we calculate the uVu
    a%uVu = 0._r8
    add_lm = 1
    do l = 0, s%ps%L_max
      do lm = -l , l
        if(l .ne. s%ps%L_loc) then
          do j = 1, a%Mps
            call mesh_r(m, a%Jxyz(j), r, x=x, a=a%x)
            ylm = oct_ylm(x(1), x(2), x(3), l, lm)
            
            if(r > 0._r8) then ! 0**l crashes in osf
              a%uVu(add_lm) = a%uVu(add_lm) + a%uV(j, add_lm)* &
                   splint(s%ps%Ur(l), r) * ylm * (r**l)
            end if
          end do
          ! uVu can be calculated exactly, or numerically
          a%uVu(add_lm) = s%ps%dkbcos(l)
          !a%uVu(add_lm) = sum(a%uV(:, add_lm)**2)/(a%uVu(add_lm)*s%ps%dknrm(l))
          if(abs(a%uVu(add_lm) - s%ps%dkbcos(l))/s%ps%dkbcos(l) > 0.25_r8) then
            write(message(1), '(a,i4)') "Low precision in the calculation of the uVu for lm = ", &
                 add_lm
            write(message(2), '(f14.6,a,f14.6)') s%ps%dkbcos(l), ' .ne. ', a%uVu(add_lm)
            message(3) = "Please consider decreasing the spacing, or changing pseudopotential"
            call write_warning(3)

            a%uVu(add_lm) = s%ps%dkbcos(l)
          end if
        end if
          
        add_lm = add_lm + 1
      end do
    end do

  end subroutine build_nl_part

end subroutine generate_external_pot

subroutine get_nl_part(ps, x, l, lm, uV, duV)
  type(ps_type), intent(IN) :: ps
  real(r8), intent(in) :: x(3)
  integer, intent(in) :: l, lm
  real(r8), intent(out) :: uV, duV(3)

  real(r8) :: r, f, uVr0, duvr0, ylm, gylm(3)

  r = sqrt(sum(x**2))
  uVr0  = splint(ps%kb(l), r)
  duvr0 = splint(ps%dkb(l), r)
        
  call grylmr(x(1), x(2), x(3), l, lm, ylm, gylm)

  select case(l)
  case(0)
    if(r >= r_small) then
      f = ylm*duvr0/r
    else
      f = 0.0_r8
    end if
    Uv = uvr0*ylm
    dUv(:) = f*x(:)
  case(1)
    Uv = uvr0 * ylm * r
    dUv(:) = duvr0*x(:)*ylm
    select case(lm)
    case(1)
      dUv(2) = dUv(2) - 0.488602511903_r8*uvr0
    case(2)
      dUv(3) = dUv(3) + 0.488602511903_r8*uvr0
    case(3)
      dUv(1) = dUv(1) - 0.488602511903_r8*uvr0
    end select
  case default
    if(r >= r_small) then
      f = ylm * (duVr0 * r**(l-1) + uVr0 * l * r**(l-2))
    else
      f = 0._r8
    end if
    
    uV = uVr0 * Ylm * (r**l)
    duV(:) = f*x(:) + uVr0*gYlm(:)*(r**l)
  end select

end subroutine get_nl_part

! this actually adds to outp
subroutine phase_factor(m, n, vec, inp, outp)
  implicit none
  type(mesh_type), intent(IN) :: m
  integer, intent(in)         :: n(3)
  real(r8), intent(IN)        :: vec(3)
  complex(r8), intent(IN)     :: inp (n(1)/2+1, n(2), n(3))
  complex(r8), intent(inout)  :: outp(n(1)/2+1, n(2), n(3))
  
  complex(r8) :: k(3)
  integer     :: ix, iy, iz, ixx, iyy, izz
  
  k(1:3) = M_zI * ((2.0_r8*M_Pi)/(n(1:3)*m%h(1:3)))
  do iz = 1, n(3)
    izz = pad_feq(iz, n(3), .true.)
    do iy = 1, n(2)
      iyy = pad_feq(iy, n(2), .true.)
      do ix = 1, n(1)/2 + 1
        ixx = pad_feq(ix, n(1), .true.)
        outp(ix, iy, iz) = outp(ix, iy, iz) + &
             exp( -(k(1)*vec(1)*ixx + k(2)*vec(2)*iyy + k(3)*vec(3)*izz) ) * inp(ix, iy, iz)
      end do
    end do
  end do
end subroutine phase_factor

subroutine generate_classic_pot(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(IN) :: sys

  integer i, ia
  real(r8) :: r

  sub_name = 'generate_classic_pot'; call push_sub()

  h%Vclassic = 0._r8
  do ia = 1, sys%ncatoms
    do i = 1, sys%m%np
      call mesh_r(sys%m, i, r, a=sys%catom(ia)%x)
      if(r < r_small) r = r_small
      h%Vclassic(i) = h%Vclassic(i) - sys%catom(ia)%charge/r
    end do
  end do

  call pop_sub()
end subroutine generate_classic_pot
