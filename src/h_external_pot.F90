! This subroutine should be in specie.F90, but due to the limitations
! of f90 to handle circular dependences it had to come here!
subroutine specie_fourier_init(ns, s, m)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)
  type(mesh_type), intent(IN) :: m

  integer :: i, j, ix, iy, iz
  real(r8) :: r, vl
  real(r8), allocatable :: fr(:,:,:)

  sub_name = 'specie_init_fourier'; call push_sub()

  allocate(fr(m%fft_n,  m%fft_n, m%fft_n))
  do i = 1, ns

    allocate(s(i)%local_fw(m%hfft_n, m%fft_n, m%fft_n))
    fr = 0.0_r8
    do j = 1, m%np
      call mesh_r(m, j, r)
      ix = m%lx(j) + m%fft_n/2 + 1;
      iy = m%ly(j) + m%fft_n/2 + 1;
      iz = m%lz(j) + m%fft_n/2 + 1;
      vl  = splint(s(i)%ps%vlocal, r)
      if(r >= r_small) then
        fr(ix, iy, iz) = (vl - s(i)%Z_val)/r
      else
        fr(ix, iy, iz) = s(i)%ps%vlocal_origin
      endif
    enddo
    call rfftwnd_f77_one_real_to_complex(m%dplanf, fr, s(i)%local_fw)
    call zscal(m%fft_n**2*m%hfft_n, cmplx(1.0_r8/m%fft_n**3, 0.0_r8, r8), s(i)%local_fw, 1)

    if(s(i)%ps%icore /= 'nc  ') then
      fr = 0._r8
      do j = 1, m%np
        call mesh_r(m, j, r)
        ix = m%lx(j) + m%fft_n/2 + 1;
        iy = m%ly(j) + m%fft_n/2 + 1;
        iz = m%lz(j) + m%fft_n/2 + 1;
        vl  = splint(s(i)%ps%core, r)
        fr(ix, iy, iz) = vl
      enddo
      call rfftwnd_f77_one_real_to_complex(m%dplanf, fr, s(i)%rhocore_fw)
      call zscal(m%fft_n**2*m%hfft_n, cmplx(1.0_r8/m%fft_n**3, 0.0_r8, r8), s(i)%rhocore_fw, 1)
    end if
  end do
  
  call pop_sub()
end subroutine specie_fourier_init

subroutine generate_external_pot(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys

  integer :: ia
  type(specie_type), pointer :: s
  type(atom_type),   pointer :: a
  complex(r8), allocatable :: fw(:,:,:), fwc(:,:,:)
  real(r8), allocatable :: fr(:,:,:)

  ! WARNING DEBUG
!!$  integer :: i, j
!!$  real(r8), allocatable :: f(:,:,:)

  sub_name = 'generate_external_pot'; call push_sub()

  ! first we assume that we need to recalculate the ion_ion energy
  call ion_ion_energy(sys)

  if(h%vpsl_space == 1) then
    allocate(fw(sys%m%hfft_n, sys%m%fft_n, sys%m%fft_n), &
         fwc(sys%m%hfft_n, sys%m%fft_n, sys%m%fft_n))
    fw = M_z0; fwc = M_z0
  end if

  h%Vpsl = 0._r8
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
    allocate(fr(sys%m%fft_n, sys%m%fft_n, sys%m%fft_n))
    call rfftwnd_f77_one_complex_to_real(sys%m%dplanb, fw, fr)
    call dcube_to_mesh(sys%m, fr, h%vpsl)
    call rfftwnd_f77_one_complex_to_real(sys%m%dplanb, fwc, fr)
    call dcube_to_mesh(sys%m, fr, h%rho_core)

    deallocate(fw, fwc, fr)
  end if

  ! WARNING DEBUG
!!$  allocate(f(sys%m%fft_n, sys%m%fft_n, sys%m%fft_n))
!!$  call dmesh_to_cube(sys%m, h%Vpsl, f)
!!$  do i = 1, sys%m%fft_n
!!$    do j = 1, sys%m%fft_n
!!$      print *, i, j, f(i, j, sys%m%fft_n/2 + 1)
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
    real(r8) :: a1, a2, Rb2, r

    a1 = s%Z/(2._r8*s%jradius**3)
    a2 = s%Z/s%jradius
    Rb2= s%jradius**2
    do i = 1, h%np
      call mesh_r(m, i, r, a=a%x)
      if(r <= s%jradius) then
        h%Vpsl(i) = h%Vpsl(i) + (a1*(r*r - Rb2) - a2)
      else
        h%Vpsl(i) = h%Vpsl(i) - s%Z/r
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

  ! this actually adds to outp
  subroutine phase_factor(m, vec, inp, outp)
    implicit none
    type(mesh_type), intent(IN) :: m
    real(r8), intent(IN)        :: vec(3)
    complex(r8), intent(IN)     :: inp(m%hfft_n, m%fft_n, m%fft_n)
    complex(r8), intent(inout)  :: outp(m%hfft_n, m%fft_n, m%fft_n)
    
    complex(r8) :: k
    integer     :: ix, iy, iz, ixx, iyy, izz
    
    k = M_zI * ((2.0_r8*M_Pi)/(m%fft_n*m%h))
    do iz = 1, m%fft_n
      izz = pad_feq(iz, m%fft_n, .true.)
      do iy = 1, m%fft_n
        iyy = pad_feq(iy, m%fft_n, .true.)
        do ix = 1, m%hfft_n
          ixx = pad_feq(ix, m%fft_n, .true.)
          outp(ix, iy, iz) = outp(ix, iy, iz) + &
               exp( k * (vec(1)*ixx + vec(2)*iyy + vec(3)*izz) ) * inp(ix, iy, iz)
        end do
      end do
    end do
  end subroutine phase_factor

  subroutine build_local_part(m)
    type(mesh_type), intent(in) :: m

    integer :: i
    real(r8) :: r, vl

    if(h%vpsl_space == 0) then
      do i = 1, m%np
        call mesh_r(m, i, r, a=a%x)
        vl  = splint(s%ps%vlocal, r)
        if(r >= r_small) then
          h%Vpsl(i) = h%Vpsl(i) + (vl - s%Z_val)/r
        else
          h%Vpsl(i) = h%Vpsl(i) + s%ps%Vlocal_origin
        end if
        if(s%ps%icore /= 'nc  ' ) then
          h%rho_core(i) = h%rho_core(i) + splint(s%ps%core, r)
        end if
      end do
    else
      call phase_factor(m, a%x, s%local_fw, fw)
      if(s%ps%icore /= 'nc  ') then
        call phase_factor(m, a%x, s%rhocore_fw, fwc)
      end if
    end if

  end subroutine build_local_part
  
  subroutine build_kb_sphere(m)
    type(mesh_type), intent(in) :: m
    
    integer :: j, k
    real(r8) :: r
    
    j = 0
    do k = 1, h%np
      call mesh_r(m, k, r, a=a%x)
      if(r < s%ps%rc_max) j = j + 1
    end do
    a%Mps = j
    
    allocate(a%Jxyz(j), a%uV(j, (s%ps%L_max+1)**2), &
         a%duV(j, (s%ps%L_max+1)**2, 3), a%uVu((s%ps%L_max+1)**2))
    
    a%uV  = 0.0_r8; a%duV = 0.0_r8; a%uVu = 0.0_r8
    
    j = 0
    do k = 1, h%np
      call mesh_r(m, k, r, a=a%x)
      if(r < s%ps%rc_max) then
        j = j + 1
        a%Jxyz(j) = k
      end if
    end do

  end subroutine build_kb_sphere

  subroutine build_nl_part(m)
    type(mesh_type), intent(IN) :: m

    integer :: j, l, k, lm
    real(r8) :: r, f, uVr0, duvr0, x(3), ylm, gylm(3)

    j_loop: do j = 1, a%Mps
      call mesh_r(m, a%Jxyz(j), r, x=x, a=a%x)
     
      k = 1
      l_loop: do l = 0, s%ps%L_max
        if(l == s%ps%L_loc) then
          k = k + (2*l + 1)
          cycle l_loop
        end if

        uVr0  = splint(s%ps%kb(l), r)
        duvr0 = splint(s%ps%dkb(l), r)
        
        lm_loop: do lm = -l , l
          call grylmr(x(1), x(2), x(3), l, lm, ylm, gylm)

          select case(l)
          case(0)
            if(r >= r_small) then
              f = ylm*duvr0/r
            else
              f = 0.0_r8
            end if
            a%uv(j, k) = uvr0*ylm
            a%duv(j, k, :) = f*x(:)
          case(1)
            a%uv(j, k) = uvr0 * ylm * r
            a%duv(j, k, :) = duvr0*x(:)*ylm
            select case(lm)
            case(1)
              a%duv(j, k, 2) = a%duv(j, k, 2) - 0.488602511903_r8*uvr0
            case(2)
              a%duv(j, k, 3) = a%duv(j, k, 3) + 0.488602511903_r8*uvr0
            case(3)
              a%duv(j, k, 1) = a%duv(j, k, 1) - 0.488602511903_r8*uvr0
            end select
          case default
            if(r >= r_small) then
              f = ylm * (duVr0 * r**(l-1) + uVr0 * l * r**(l-2))
            else
              f = 0._r8
            end if
            
            a%uV(j, k) = uVr0 * Ylm * (r**l)
            a%duV(j, k, :) = f*x(:) + uVr0*gYlm(:)*(r**l)
          end select
            
          ! now for the uVus
          a%uVu(k) = a%uVu(k) + a%uV(j, k)* &
               splint(s%ps%Ur(l), r) * ylm * (r**l)

          k = k + 1
        end do lm_loop

      end do l_loop
    end do j_loop

    ! now we finish calculating the uVus
    
    k = 1
    do l = 0 , s%ps%L_max
      do lm = -l , l
        if(l.ne.s%ps%L_loc) then
          ! uVu can be calculated exactly, or numerically
          a%uVu(k) = s%ps%dkbcos(l)
          !a%uVu(k) = sum(a%uV(:, k)*a%uV(:, k))/(a%uVu(k)*s%ps%dknrm(l))
        end if
        k = k + 1
      end do
    end do
  end subroutine build_nl_part
  
end subroutine generate_external_pot
