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

subroutine specie_local_fourier_init(ns, s, m, fft, nlcc)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)
  type(mesh_type), intent(IN) :: m
  type(fft_type), intent(in) :: fft
  logical, intent(in) :: nlcc

  integer :: i, j, ix, iy, iz, n, ixx(3), db(3), dbc(3), c(3)
  real(r8) :: x(3), g(3), modg, temp(3)
  real(r8), allocatable :: fr(:,:,:), v(:)
  !complex(r8) :: c

  call push_sub('specie_local_fourier_init')

  call fft_getdim_real   (fft, db)  ! get dimensions of real array
  call fft_getdim_complex(fft, dbc) ! get dimensions of complex array
  c(:) = db(:)/2 + 1                ! get center of double box

  allocate(fr(db(1), db(2), db(3)))

  do i = 1, ns
    allocate(s(i)%local_fw(dbc(1), dbc(2), dbc(3)))

    if (conf%periodic_dim==0) then
      do ix = 1, db(1)
        ixx(1) = ix - c(1)
        do iy = 1, db(2)
          ixx(2) = iy - c(2)
          do iz = 1, db(3)
            ixx(3) = iz - c(3)
            
            x(:) = m%h(:)*ixx(:)
            fr(ix, iy, iz) = specie_get_local(s(i), x)
          end do
        end do
      end do
      call rfft_forward(fft, fr, s(i)%local_fw)
    else
      allocate(v(2:s(i)%ps%g%nrval))
      temp(:) = M_TWO*M_PI/(db(:)*m%h(:))
      do ix = 1, dbc(1)
        ixx(1) = pad_feq(ix, db(1), .true.)
        do iy = 1, dbc(2)
          ixx(2) = pad_feq(iy, db(2), .true.)
          do iz = 1, dbc(3)
            ixx(3) = pad_feq(iz, db(3), .true.)
            modg = sqrt(sum((temp(:)*ixx(:))**2))
            if(modg.ne.M_ZERO) then
              do j = 2, s(i)%ps%g%nrval
                v(j) = (sin(modg*s(i)%ps%g%rofi(j))/(modg*s(i)%ps%g%rofi(j)))*     &
                     s(i)%ps%g%rofi(j)**2*(splint(s(i)%ps%vlocal,s(i)%ps%g%rofi(j)))
              enddo
              s(i)%local_fw(ix, iy, iz) = M_FOUR*M_PI*    &
                   (sum(s(i)%ps%g%drdi(2:s(i)%ps%g%nrval)*v(2:s(i)%ps%g%nrval))-s(i)%ps%z_val/modg)
            else
              s(i)%local_fw(ix, iy, iz) = M_ZERO
            end if
          end do
        end do
      end do
      deallocate(v)
    end if
    
    ! now we built the non-local core corrections in momentum space
    if(nlcc) then
      allocate(s(i)%rhocore_fw(dbc(1), dbc(2), dbc(3)))
      
      fr = 0.0_r8
      do ix = 1, db(1)
        ixx(1) = ix - c(1)
        do iy = 1, db(2)
          ixx(2) = iy - c(2)
          do iz = 1, db(3)
            ixx(3) = iz - c(3)
            
            x(:) = m%h(:)*ixx(:)
            fr(ix, iy, iz) = specie_get_nlcc(s(i), x)
          end do
        end do
      end do

      call rfft_forward(fft, fr, s(i)%rhocore_fw)
    end if

  end do
  
  deallocate(fr)
  call pop_sub()
end subroutine specie_local_fourier_init

subroutine generate_external_pot(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys

  integer :: ia, db(3), dbc(3)
  type(specie_type), pointer :: s
  type(atom_type),   pointer :: a
  real(r8), allocatable :: fr(:,:,:)
  complex(r8), allocatable :: fw(:,:,:)
  complex(r8), allocatable :: fwc(:,:,:) ! for the nl core corrections

  call push_sub('generate_external_pot')

  ! first we assume that we need to recalculate the ion_ion energy
  sys%eii = ion_ion_energy(sys%natoms, sys%atom)

  if(h%vpsl_space == RECIPROCAL_SPACE) then
    call fft_getdim_real   (h%fft, db)  ! get dimensions of real array
    call fft_getdim_complex(h%fft, dbc) ! get dimensions of complex array
    
    allocate(fw(dbc(1), dbc(2), dbc(3)))
    fw = M_z0

    if(sys%nlcc) then
      allocate(fwc(dbc(1), dbc(2), dbc(3)))
      fwc = M_z0
    end if
  end if

  h%Vpsl = M_ZERO
  if(sys%nlcc) sys%st%rho_core = M_ZERO

  do ia = 1, sys%natoms
    a => sys%atom(ia) ! shortcuts
    s => a%spec
    
    call build_local_part(sys%m)
    
    if(.not.s%local) then
      call build_kb_sphere(sys%m)
      call build_nl_part(sys%m)
    end if

  end do

  if(h%vpsl_space == RECIPROCAL_SPACE) then
    allocate(fr(db(1), db(2), db(3)))

    ! first the potential
    call rfft_backward(h%fft, fw, fr)
    call dcube_to_mesh(sys%m, fr, h%Vpsl, db)
      
    ! and the non-local core corrections
    if(sys%nlcc) then
      call rfft_backward(h%fft, fwc, fr)
      call dcube_to_mesh(sys%m, fr, sys%st%rho_core, db)
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

    call push_sub('build_local_part')
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
      call phase_factor(m, db, a%x, s%local_fw, fw)
      if(s%nlcc) then
        call phase_factor(m, db, a%x, s%rhocore_fw, fwc)
      end if
    end if

    call pop_sub()
  end subroutine build_local_part

  subroutine build_kb_sphere(m)
    type(mesh_type), intent(in) :: m
    
    integer :: i, j, k
    real(r8) :: r, x(3)

    call push_sub('build_kb_sphere')

    ! This is for the ions movement; probably it is not too elegant, I will rethink it later.
    if(associated(a%jxyz)) deallocate(&
         a%jxyz, a%duv,  a%dduv,  a%duvu,     &
         a%zuv, a%zduv, a%zuvu,       &
         a%so_uv, a%so_duv, a%so_uvu, &
         a%so_luv)
    if(conf%periodic_dim/=0 .and. associated(a%phases)) deallocate(a%phases)
    if (any(s%ps%rc_max + m%h(1)>=m%lsize(1:conf%periodic_dim))) then
      message(1)='KB sphere is larger than the box size'
      write(message(2),'(a,f12.6,a)')'  rc_max+h = ',s%ps%rc_max + m%h(1),' [b]'
      write(message(3),'(a,3f12.4,a)')'  lsize    = ',m%lsize,' [b]'
      message(4)='Please change pseudopotential'
      call write_fatal(4)
    end if
    j = 0
    do k = 1, h%np
      do i=1,3**conf%periodic_dim
        call mesh_r(m, k, r, a=a%x+m%shift(i,:))
        if(r > s%ps%rc_max + m%h(1)) cycle
        j = j + 1
        exit
      end do
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
      do i=1,3**conf%periodic_dim
        call mesh_r(m, k, r, a=a%x+m%shift(i,:))
        ! we enlarge slightly the mesh (good for the interpolation scheme)
        if(r > s%ps%rc_max + m%h(1)) cycle
        j = j + 1
        a%Jxyz(j) = k
        exit
      end do
    end do

    call pop_sub()
  end subroutine build_kb_sphere

  subroutine build_nl_part(m)
    type(mesh_type), intent(IN) :: m

    integer :: i, j, k, l, lm, n, add_lm, p, ix(3), center(3)
    real(r8) :: r, x(3), x_in(3), ylm
    real(r8) :: so_uv, so_duv(3)

    call push_sub('build_nl_part')

    ! This loop is done always, for spin-orbit is, up to now, only done in real space.
    ! If we want the nl part also in real space, it is read.
    j_loop: do j = 1, a%Mps
      call mesh_xyz(m, a%Jxyz(j), x_in)
      k_loop: do k=1,3**conf%periodic_dim
        x(:) = x_in(:) - m%shift(k,:)          
        r=sqrt(sum((x-a%x)*(x-a%x)))
        if (r > s%ps%rc_max + m%h(1)) cycle
        x = x - a%x
        add_lm = 1
        l_loop: do l = 0, s%ps%L_max
          lm_loop: do lm = -l, l
            i_loop : do i = 1, s%ps%kbc
              if(l .ne. s%ps%L_loc .and. h%vnl_space == REAL_SPACE) then
                call specie_get_nl_part(s, x, l, lm, i, a%duV(j, add_lm, i), a%dduV(:, j, add_lm, i))
              end if
              if(l>0 .and. s%ps%so_l_max>=0) then
                call specie_get_nl_part(s, x, l, lm, i, so_uv, so_duv(:), so=.true.)
                a%so_uv(j, add_lm, i) = so_uv
                a%so_duv(1:3, j, add_lm, i) = so_duv(1:3)
              endif
            end do i_loop
            add_lm = add_lm + 1
          end do lm_loop
        end do l_loop
        exit
      end do k_loop
    end do j_loop

    if(h%reltype == SPIN_ORBIT) then
      do j = 1, a%mps
        call mesh_xyz(m, a%Jxyz(j), x_in)
        do k=1,3**conf%periodic_dim
          x(:) = x_in(:) - m%shift(k,:)          
          r=sqrt(sum((x-a%x)*(x-a%x)))
          if (r > s%ps%rc_max + m%h(1)) cycle
          x = x - a%x
          do add_lm = 1, (a%spec%ps%l_max+1)**2
            a%so_luv(j, add_lm, 1:a%spec%ps%kbc, 1) = &
                  x(2)*a%so_duv(3, j, add_lm, 1:a%spec%ps%kbc) - &
                  x(3)*a%so_duv(2, j, add_lm, 1:a%spec%ps%kbc)
            a%so_luv(j, add_lm, 1:a%spec%ps%kbc, 2) = &
                  x(3)*a%so_duv(1, j, add_lm, 1:a%spec%ps%kbc) - &
                  x(1)*a%so_duv(3, j, add_lm, 1:a%spec%ps%kbc)
            a%so_luv(j, add_lm, 1:a%spec%ps%kbc , 3) = &
                  x(1)*a%so_duv(2, j, add_lm, 1:a%spec%ps%kbc ) - &
                  x(2)*a%so_duv(1, j, add_lm, 1:a%spec%ps%kbc )
          enddo
          exit
        enddo
      enddo
      a%so_luv = -M_zI*a%so_luv
    endif

    ! and here we calculate the uVu
    if(s%ps%flavour(1:2) == 'tm') then
      a%duVu = M_ZERO
      add_lm = 1
      do l = 0, s%ps%L_max
        do lm = -l , l
          if(l .ne. s%ps%L_loc) then
            do j = 1, a%Mps
              call mesh_r(m, a%Jxyz(j), r, x=x_in, a=a%x)
              do k=1,3**conf%periodic_dim
                x(:) = x_in(:) - m%shift(k,:)          
                r=sqrt(sum(x*x))
                if (r > s%ps%rc_max + m%h(1)) cycle
                ylm = oct_ylm(x(1), x(2), x(3), l, lm)
                a%duvu(add_lm, 1, 1) = a%duvu(add_lm, 1, 1) + a%duv(j, add_lm, 1)* &
                                       splint(s%ps%ur(l+1, 1), r)*ylm*m%vol_pp
                exit
              end do
            end do
            a%duvu(add_lm, 1, 1) = M_ONE/(a%duvu(add_lm, 1, 1)*s%ps%dknrm(l))
            if(abs((a%duVu(add_lm, 1, 1) - s%ps%h(l,1,1))/s%ps%h(l,1,1)) > 0.05_r8 .and. s%ps%ispin==1) then
              write(message(1), '(a,i4)') "Low precision in the calculation of the uVu for lm = ", &
                   add_lm
              write(message(2), '(f14.6,a,f14.6)') s%ps%h(l,1,1), ' .ne. ', a%duVu(add_lm, 1, 1)
              message(3) = "Please consider decreasing the spacing, or changing pseudopotential"
              call write_warning(3)
            end if
            ! uVu is in fact calculated in the logarithmic grid, previous stuff is a check.
            a%duvu(add_lm, 1, 1) = s%ps%h(l, 1, 1)
          end if
          a%so_uvu(add_lm, 1, 1) = s%ps%k(l, 1, 1)
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

    call pop_sub()
  end subroutine build_nl_part

end subroutine generate_external_pot

subroutine generate_classic_pot(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(IN) :: sys

  integer i, ia
  real(r8) :: r, rc

  call push_sub('generate_classic_pot')

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
