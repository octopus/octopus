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

! This subroutine should be in specie.F90, but due to the limitations
! of f90 to handle circular dependences it had to come here!

subroutine specie_local_fourier_init(ns, s, m, nlcc)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)
  type(mesh_type), intent(IN) :: m
  logical, intent(in) :: nlcc

  integer :: i, j, ix, iy, iz, n, ixx(3)
  real(r8) :: x(3), g(3), g2
  real(r8), allocatable :: fr(:,:,:)
  complex(r8) :: c

  sub_name = 'specie_local_fourier_init'; call push_sub()

  allocate(fr(m%fft_n2(1),  m%fft_n2(2), m%fft_n2(3)))

  do i = 1, ns
    allocate(s(i)%local_fw(m%hfft_n2, m%fft_n2(2), m%fft_n2(3)))

    do ix = 1, m%fft_n2(1)
      ixx(1) = ix - (m%fft_n2(1)/2 + 1)
      do iy = 1, m%fft_n2(2)
        ixx(2) = iy - (m%fft_n2(2)/2 + 1)
        do iz = 1, m%fft_n2(3)
          ixx(3) = iz - (m%fft_n2(3)/2 + 1)

          x(:) = m%h(:)*ixx(:)
          fr(ix, iy, iz) = specie_get_local(s(i), x)
        end do
      end do
    end do

    call rfftwnd_f77_one_real_to_complex(m%dplanf2, fr, s(i)%local_fw)
    
    n = m%fft_n2(2)*m%fft_n2(3)*m%hfft_n2
    c  = cmplx(1.0_r8/(m%fft_n2(1)*m%fft_n2(2)*m%fft_n2(3)), 0.0_r8, r8)
    call zscal(n, c, s(i)%local_fw,  1)
    
    ! now we built the non-local core corrections in momentum space
    if(nlcc) then
      allocate(s(i)%rhocore_fw(m%hfft_n2, m%fft_n2(2), m%fft_n2(3)))

      fr = 0.0_r8
      do ix = 1, m%fft_n2(1)
        ixx(1) = ix - (m%fft_n2(1)/2 + 1)
        do iy = 1, m%fft_n2(2)
          ixx(2) = iy - (m%fft_n2(2)/2 + 1)
          do iz = 1, m%fft_n2(3)
            ixx(3) = iz - (m%fft_n2(3)/2 + 1)

            x(:) = m%h(:)*ixx(:)
            fr(ix, iy, iz) = specie_get_nlcc(s(i), x)
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

  integer :: n(3), hn, i, ii, j, ix, iy, iz, ixx, iyy, izz, l, lm, add_lm, kx, ky, kz
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
    allocate(s(i)%nl_fw(hn, n(2), n(3), (s(i)%ps%L_max + 1)**2, s(i)%ps%kbc), &
        s(i)%nl_dfw(hn, n(2), n(3), 3, (s(i)%ps%L_max + 1)**2, s(i)%ps%kbc))
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
      ii_loop : do ii = 1, s(i)%ps%kbc
        do ix = 1, n(1)
          x(1) = m%h(1) * real(ix - n(1)/2 - 1, r8) / real(1 + nextra, r8)
          do iy = 1, n(2)
            x(2) = m%h(2) * real(iy - n(2)/2 - 1, r8) / real(1 + nextra, r8)
            do iz = 1, n(3)
              x(3) = m%h(3) * real(iz - n(3)/2 - 1, r8) / real(1 + nextra, r8)
              r = sqrt(sum(x**2))
              
              call specie_get_nl_part(s(i), x, l, lm, ii, fr(ix, iy, iz), dfr(ix, iy, iz, :))
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
              s(i)%nl_fw (ix, iy, iz, add_lm, ii)   = c * fw (ixx, iyy, izz)
              s(i)%nl_dfw(ix, iy, iz,:, add_lm, ii) = c * dfw(ixx, iyy, izz, :)
            end do
          end do
        end do
      
      end do ii_loop  
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
  real(r8), allocatable :: fr(:,:,:)
  complex(r8), allocatable :: fw(:,:,:)
  complex(r8), allocatable :: fwc(:,:,:) ! for the nl core corrections
  sub_name = 'generate_external_pot'; call push_sub()

  ! first we assume that we need to recalculate the ion_ion energy
  sys%eii = ion_ion_energy(sys%natoms, sys%atom)

  if(h%vpsl_space == RECIPROCAL_SPACE) then
    allocate(fw(sys%m%hfft_n2, sys%m%fft_n2(2), sys%m%fft_n2(3)))
    fw = M_z0
  end if
  h%Vpsl  = 0._r8

  if(sys%nlcc) then
    sys%st%rho_core = 0._r8
    allocate(fwc(sys%m%hfft_n2, sys%m%fft_n2(2), sys%m%fft_n2(3)))
    fwc = M_z0
  end if

  do ia = 1, sys%natoms
    a => sys%atom(ia) ! shortcuts
    s => a%spec
    
    call build_local_part(sys%m)

    if(.not.s%local) then
      call build_kb_sphere(sys%m)
      call build_nl_part(sys%m)
    end if

  end do

  if(h%vpsl_space == 1) then
    allocate(fr(sys%m%fft_n2(1), sys%m%fft_n2(2), sys%m%fft_n2(3)))

    ! first the potential
    call rfftwnd_f77_one_complex_to_real(sys%m%dplanb2, fw, fr)
    call dcube_to_mesh(sys%m, fr, h%Vpsl, t=2)

    ! and the non-local core corrections
    if(sys%nlcc) then
      call rfftwnd_f77_one_complex_to_real(sys%m%dplanb2, fwc, fr)
      call dcube_to_mesh(sys%m, fr, sys%st%rho_core, 2)
      deallocate(fwc)
    end if
    
    deallocate(fw, fr)
  end if

  if (h%classic_pot > 0) then
    h%Vpsl = h%Vpsl + h%Vclassic
  end if

  call pop_sub()

contains
  subroutine build_local_part(m)
    type(mesh_type), intent(in) :: m
    
    integer :: i
    real(r8) :: x(3), r

    sub_name = 'build_local_part'; call push_sub()
    if(h%vpsl_space == REAL_SPACE) then ! real space
      do i = 1, h%np
        call mesh_xyz(m, i, x)
        x = x - a%x
        h%Vpsl(i) = h%Vpsl(i) + specie_get_local(s, x)
        
        if(s%nlcc) then
          sys%st%rho_core(i) = sys%st%rho_core(i) + specie_get_nlcc(s, x)
        end if
      end do

    else ! momentum space
      call phase_factor(m, m%fft_n2, a%x, s%local_fw, fw)
      if(s%nlcc) then
        call phase_factor(m, m%fft_n2(1:3), a%x, s%rhocore_fw, fwc)
      end if
    end if

    call pop_sub(); return
  end subroutine build_local_part

  subroutine build_kb_sphere(m)
    type(mesh_type), intent(in) :: m
    
    integer :: j, k
    real(r8) :: r

    sub_name = 'build_kb_sphere'; call push_sub()

    ! This is for the ions movement; probably it is not too elegant, I will rethink it later.
    if(associated(a%jxyz)) deallocate(a%jxyz, a%duv,  a%dduv,  a%duvu,     &
                                              a%zuv, a%zduv, a%zuvu,       &
                                              a%so_uv, a%so_duv, a%so_uvu, a%so_luv)
    
    j = 0
    do k = 1, h%np
      call mesh_r(m, k, r, a=a%x)
      if(r <= s%ps%rc_max + m%h(1)) j = j + 1
    end do
    a%Mps = j
    allocate(a%Jxyz(j), a%duV(j, (s%ps%L_max+1)**2, s%ps%kbc), &
         a%dduV(3, j, (s%ps%L_max+1)**2, s%ps%kbc), a%duVu((s%ps%L_max+1)**2, s%ps%kbc, s%ps%kbc))
    allocate(a%zuV(j, (s%ps%L_max+1)**2, s%ps%kbc), &
         a%zduV(3, j, (s%ps%L_max+1)**2, s%ps%kbc), a%zuVu((s%ps%L_max+1)**2, s%ps%kbc, s%ps%kbc))
    allocate(a%so_uV(j, (s%ps%L_max+1)**2, s%ps%kbc), &
         a%so_duV(3, j, (s%ps%L_max+1)**2, s%ps%kbc), a%so_uVu((s%ps%L_max+1)**2, s%ps%kbc, s%ps%kbc), &
         a%so_luv(j, (s%ps%L_max+1)**2, s%ps%kbc, 3))

    a%duv    = 0.0_r8; a%dduV   = 0.0_r8; a%duVu   = 0.0_r8
    a%zuv    = M_z0;   a%zduV   = M_z0;   a%zuVu   = M_z0
    a%so_uv  = M_z0;   a%so_duV = M_z0;   a%so_uvu = M_z0; a%so_luv = M_z0

    j = 0
    do k = 1, h%np
      call mesh_r(m, k, r, a=a%x)
      ! we enlarge slightly the mesh (good for the interpolation scheme)
      if(r <= s%ps%rc_max + m%h(1)) then
        j = j + 1
        a%Jxyz(j) = k
      end if
    end do

    call pop_sub(); return
  end subroutine build_kb_sphere

  subroutine build_nl_part(m)
    type(mesh_type), intent(IN) :: m

    integer :: i, j, l, lm, add_lm, p, ix(3), center(3)
    real(r8) :: r, x(3), ylm
    complex(r8), allocatable :: nl_fw(:,:,:)
    real(r8), allocatable :: nl_fr(:,:,:), nl_dfr(:,:,:,:)
    real(r8) :: so_uv, so_duv(3)

    sub_name = 'build_nl_part'; call push_sub()

    ! This loop is done always, for spin-orbit is, up to now, only done in real space.
    ! If we want the nl part also in real space, it is read.
    j_loop: do j = 1, a%Mps
        call mesh_xyz(m, a%Jxyz(j), x)
        x = x - a%x
        
        add_lm = 1
        l_loop: do l = 0, s%ps%L_max
          lm_loop: do lm = -l , l
            i_loop : do i = 1, s%ps%kbc
              if(l .ne. s%ps%L_loc .and. h%vnl_space == REAL_SPACE) then
                 call specie_get_nl_part(s, x, l, lm, i, a%duV(j, add_lm, i), a%dduV(:, j, add_lm, i))
              end if
              call specie_get_nl_part(s, x, l, lm, i, so_uv, so_duv(:), so=.true.)
              a%so_uv(j, add_lm, i) = so_uv
              a%so_duv(:, j, add_lm, i) = so_duv(:)

            end do i_loop
            add_lm = add_lm + 1
          end do lm_loop
        end do l_loop
    end do j_loop

    if(h%reltype == SPIN_ORBIT) then
      do j = 1, a%mps
        call mesh_xyz(m, a%Jxyz(j), x)
        x = x - a%x
        a%so_luv(j, 1:(a%spec%ps%l_max+1)**2, 1:a%spec%ps%kbc, 1) = &
                x(2)*a%so_duv(3, j, 1:(a%spec%ps%l_max+1)**2, 1:a%spec%ps%kbc) - &
                x(3)*a%so_duv(2, j, 1:(a%spec%ps%l_max+1)**2, 1:a%spec%ps%kbc)
        a%so_luv(j, 1:(a%spec%ps%l_max+1)**2, 1:a%spec%ps%kbc, 2) = &
                x(3)*a%so_duv(1, j, 1:(a%spec%ps%l_max+1)**2, 1:a%spec%ps%kbc) - &
                x(1)*a%so_duv(3, j, 1:(a%spec%ps%l_max+1)**2, 1:a%spec%ps%kbc)
        a%so_luv(j, 1:(a%spec%ps%l_max+1)**2, 1:a%spec%ps%kbc , 3) = &
                x(1)*a%so_duv(2, j, 1:(a%spec%ps%l_max+1)**2, 1:a%spec%ps%kbc ) - &
                x(2)*a%so_duv(1, j, 1:(a%spec%ps%l_max+1)**2, 1:a%spec%ps%kbc )
      enddo
      a%so_luv = -M_zI*a%so_luv
    endif

    if(h%vnl_space == RECIPROCAL_SPACE) then
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
        i_loop2 : do i = 1, s%ps%kbc
          nl_fw = M_z0
          call phase_factor(m, s%nl_fft_n(1:3), x, s%nl_fw(:,:,:, add_lm, i), nl_fw)
          call rfftwnd_f77_one_complex_to_real(s%nl_planb, nl_fw, nl_fr)

          ! now the gradient
          do j = 1, 3
            nl_fw = M_z0
            call phase_factor(m, s%nl_fft_n(1:3), x, s%nl_dfw(:,:,:, j, add_lm, i), nl_fw)
            call rfftwnd_f77_one_complex_to_real(s%nl_planb, nl_fw, nl_dfr(1, 1, 1, j))
          end do
          
          j_loop2: do j = 1, a%Mps
            p = a%Jxyz(j)
            ix(:) = m%Lxyz(:,p) - center(:) + s%nl_fft_n(:)/2 + 1
            
            a%duV(j, add_lm, i)     = nl_fr (ix(1), ix(2), ix(3))
            a%dduV(:, j, add_lm, i) = nl_dfr(ix(1), ix(2), ix(3), :)
          end do j_loop2

        end do i_loop2
        add_lm = add_lm + 1
        end do lm_loop2
      end do l_loop2
      deallocate(nl_fw, nl_fr, nl_dfr)
    end if
    
    ! and here we calculate the uVu
    if(s%ps%flavour(1:2) == 'tm') then
      a%duVu = 0._r8
      add_lm = 1
      do l = 0, s%ps%L_max
        do lm = -l , l
          if(l .ne. s%ps%L_loc) then
            do j = 1, a%Mps
              call mesh_r(m, a%Jxyz(j), r, x=x, a=a%x)
              ylm = oct_ylm(x(1), x(2), x(3), l, lm)
              
              if(r > 0._r8 .or. l>0) then ! 0**l crashes in osf
                a%duVu(add_lm, 1, 1) = a%duVu(add_lm, 1, 1) + a%duV(j, add_lm, 1)* &
                     splint(s%ps%ur(l, 1), r) * ylm * (r**l)
              else
                a%duvu(add_lm, 1, 1) = a%duvu(add_lm, 1, 1) + a%duv(j, add_lm, 1)* &
                     splint(s%ps%ur(l, 1), r)*ylm
              endif
            end do
            a%duVu(add_lm, 1, 1) = sum(a%duV(:, add_lm, 1)**2)/(a%duVu(add_lm, 1, 1)*s%ps%dknrm(l))
             if(abs((a%duVu(add_lm, 1, 1) - s%ps%h(l,1,1))/s%ps%h(l,1,1)) > 0.05_r8) then
              write(message(1), '(a,i4)') "Low precision in the calculation of the uVu for lm = ", &
                   add_lm
              write(message(2), '(f14.6,a,f14.6)') s%ps%h(l,1,1), ' .ne. ', a%duVu(add_lm, 1, 1)
              message(3) = "Please consider decreasing the spacing, or changing pseudopotential"
              call write_warning(3)
            end if
            ! uVu can be calculated exactly, or numerically
            a%duvu(add_lm, 1, 1) = s%ps%h(l, 1, 1)
            a%so_uvu(add_lm, 1, 1) = s%ps%k(l, 1, 1)
          end if
          
          add_lm = add_lm + 1
        end do
      end do

    else
      add_lm = 1
      do l = 0, s%ps%l_max
      do lm = -l, l
         a%duvu(add_lm, 1:3, 1:3) = s%ps%h(l, 1:3, 1:3)
         a%so_uvu(add_lm, 1:3, 1:3) = s%ps%k(l, 1:3, 1:3)
         add_lm = add_lm + 1
      enddo
      enddo

    end if

    a%zuv = cmplx(a%duv); a%zuvu = cmplx(a%duvu); a%zduv = cmplx(a%dduv)

    call pop_sub(); return
  end subroutine build_nl_part

end subroutine generate_external_pot

subroutine generate_classic_pot(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(IN) :: sys

  integer i, ia
  real(r8) :: r, rc

  sub_name = 'generate_classic_pot'; call push_sub()

  h%Vclassic = 0._r8
  do ia = 1, sys%ncatoms
    do i = 1, sys%m%np
      call mesh_r(sys%m, i, r, a=sys%catom(ia)%x)
      select case(h%classic_pot)
      case(1) ! point charge
        if(r < r_small) r = r_small
        h%Vclassic(i) = h%Vclassic(i) - sys%catom(ia)%charge/r
      case(2) ! gaussion smeared charge
        select case(sys%catom(ia)%label(1:1)) ! covalent radii
        case('H')
          rc = 0.4_r8*P_Ang
        case('C') 
          rc = 0.8_r8*P_Ang
        case default
          rc = 0.7_r8*P_Ang
        end select
        if(abs(r - rc) < r_small) r = rc + sign(r_small, r-rc)
        h%Vclassic(i) = h%Vclassic(i) - sys%catom(ia)%charge*(r**4 - rc**4)/(r**5 - rc**5)
      end select
    end do
  end do

  call pop_sub()
end subroutine generate_classic_pot
