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

#include "global.h"

module external_pot
  use system
  use lasers

  implicit none

  type epot_type
    integer :: vpsl_space           ! How should the local potential be calculated
    integer :: vnl_space            ! How should the nl    potential be calculated
    integer :: nextra               ! extra points for the interpolation method(s)
    
    integer :: classic_pot          ! How to include the classic charges
    real(r8), pointer :: Vclassic(:) => NULL()  ! We use it to store the potential of the classic charges
 
    ! lasers stuff
    integer :: no_lasers ! number of laser pulses used
    type(laser_type), pointer :: lasers(:) => NULL()    
   
#ifdef HAVE_FFT
    ! For the local pseudopotential in Fourier space...
    type(dcf), pointer :: local_cf(:)   => NULL()
    type(dcf), pointer :: rhocore_cf(:) => NULL() ! and for the core density
#endif
  end type epot_type

contains

  subroutine epot_init(ep, sys)
    type(epot_type),   intent(out) :: ep
    type(system_type), intent(in)  :: sys

    ! should we calculate the local pseudopotentials in Fourier space?
    ep%vpsl_space = REAL_SPACE
    ep%vnl_space  = REAL_SPACE
#ifdef HAVE_FFT
    call oct_parse_int('LocalPotentialSpace', RECIPROCAL_SPACE, ep%vpsl_space)
#endif
    
    select case(ep%vpsl_space)
    case(RECIPROCAL_SPACE)
      message(1) = 'Info: Local Potential in Reciprocal Space.'
    case(REAL_SPACE)
      if (conf%periodic_dim==0) then
        message(1) = 'Info: Local Potential in Real Space.'
      else
        message(1) = 'for periodic systems you must set LocalPotentialSpace = 1'
        call write_fatal(1)
      end if
    case default
      write(message(1), '(a,i5,a)') "Input: '", ep%vpsl_space, &
           "' is not a valid LocalPotentialSpace"
      message(2) = '(LocalPotentialSpace = 0 | 1)'
      call write_fatal(2)
    end select
    call write_info(1)

#ifdef HAVE_FFT
    if(ep%vpsl_space == RECIPROCAL_SPACE) call epot_local_fourier_init(ep, sys%m, sys)
#endif

    ep%classic_pot = 0
    if(sys%ncatoms > 0) then
      call oct_parse_int("ClassicPotential", 0, ep%classic_pot)
      if(ep%classic_pot > 0) then
        allocate(ep%Vclassic(sys%m%np))
        call epot_generate_classic(ep, sys%m, sys)
      end if
    end if

    ! lasers
    call laser_init(sys%m, ep%no_lasers, ep%lasers)
    if(ep%no_lasers>0 ) then
      message(1) = 'Info: Lasers'
      call write_info(1)
      if(conf%verbose > 20 .and. mpiv%node == 0) then
        call laser_write_info(ep%no_lasers, ep%lasers, stdout)
      end if
    end if

  end subroutine epot_init

  subroutine epot_end(ep, sys)
    type(epot_type),   intent(inout) :: ep
    type(system_type), intent(in)    :: sys

#ifdef HAVE_FFT
    integer :: i

    if(ep%vpsl_space == RECIPROCAL_SPACE) then
      do i = 1, sys%nspecies
        call dcf_free(ep%local_cf(i))
        if(sys%nlcc) call dcf_free(ep%rhocore_cf(i))
      end do
      deallocate(ep%local_cf)
      if(sys%nlcc) deallocate(ep%rhocore_cf)
    end if
#endif

    if(ep%classic_pot > 0) then
      ep%classic_pot = 0
      ASSERT(associated(ep%Vclassic)) ! sanity check
      deallocate(ep%Vclassic)         ! and clean up
      nullify(ep%Vclassic)
    end if
    
    call laser_end(ep%no_lasers, ep%lasers)

  end subroutine epot_end

#ifdef HAVE_FFT
  subroutine epot_local_fourier_init(ep, m, sys)
    type(epot_type),   intent(inout) :: ep
    type(mesh_type),   intent(in)    :: m
    type(system_type), intent(in)    :: sys
    
    integer :: i, j, ix, iy, iz, ixx(3), db(3), c(3)
    real(r8) :: x(3), modg, temp(3)
    real(r8), allocatable :: v(:)
    
    type(specie_type), pointer :: s ! shortcuts
    type(dcf), pointer :: cf
    
    call push_sub('specie_local_fourier_init')

    allocate(ep%local_cf(sys%nspecies))

    specie: do i = 1, sys%nspecies
      s  => sys%specie(i)
      cf => ep%local_cf(i)
      
      if(i == 1) then
        call mesh_double_box(m, db)
        call dcf_new(db, cf)    ! initialize the cube
        call dcf_fft_init(cf)   ! and initialize the ffts
        db = cf%n               ! dimensions of the box may have been optimized, so get them
        c(:) = db(:)/2 + 1      ! get center of double box
        
        if(sys%nlcc) then                      ! if we have non-linear core corrections
          call dcf_new(db, ep%rhocore_cf(i))
          call dcf_fft_init(ep%rhocore_cf(i))
        end if
      else
        call dcf_new_from(cf, ep%local_cf(1))   ! we can just copy from the first one
        if(sys%nlcc) call dcf_new_from(ep%rhocore_cf(i), ep%rhocore_cf(1))
      end if
      
      periodic: if (conf%periodic_dim==0) then
        call dcf_alloc_RS(cf)                  ! allocate the cube in real space
        
        do iz = 1, db(3)
          ixx(3) = iz - c(3)
          do iy = 1, db(2)
            ixx(2) = iy - c(2)
            do ix = 1, db(1)
              ixx(1) = ix - c(1)
              
              x(:) = m%h(:)*ixx(:)
              cf%RS(ix, iy, iz) = specie_get_local(s, x)
            end do
          end do
        end do
        
        call dcf_alloc_FS(cf)      ! allocate the tube in Fourier space
        call dcf_RS2FS(cf)         ! Fourier transform
        call dcf_free_RS(cf)       ! we do not need the real space any longer
      else
        call dcf_alloc_FS(cf)      ! allocate the tube in Fourier space
        
        allocate(v(2:s%ps%g%nrval))
        temp(:) = M_TWO*M_PI/(db(:)*m%h(:))
        do ix = 1, cf%nx
          ixx(1) = pad_feq(ix, db(1), .true.)
          do iy = 1, db(2)
            ixx(2) = pad_feq(iy, db(2), .true.)
            do iz = 1, db(3)
              ixx(3) = pad_feq(iz, db(3), .true.)
              
              modg = sqrt(sum((temp(:)*ixx(:))**2))
              if(modg.ne.M_ZERO) then
                do j = 2, s%ps%g%nrval
                  v(j) = (sin(modg*s%ps%g%rofi(j))/(modg*s%ps%g%rofi(j)))*     &
                       s%ps%g%rofi(j)**2*(splint(s%ps%vlocal,s%ps%g%rofi(j)))
                enddo
                cf%FS(ix, iy, iz) = M_FOUR*M_PI*    &
                     (sum(s%ps%g%drdi(2:s%ps%g%nrval)*v(2:s%ps%g%nrval))-s%ps%z_val/modg)
              else
                cf%FS(ix, iy, iz) = M_ZERO
              end if
            end do
          end do
        end do
        deallocate(v)
        
      end if periodic
    
      ! now we built the non-local core corrections in momentum space
      nlcc: if(sys%nlcc) then
        call dcf_alloc_RS(ep%rhocore_cf(i))
        
        do ix = 1, db(1)
          ixx(1) = ix - c(1)
          do iy = 1, db(2)
            ixx(2) = iy - c(2)
            do iz = 1, db(3)
              ixx(3) = iz - c(3)
              
              x(:) = m%h(:)*ixx(:)
              if(sys%nlcc) ep%rhocore_cf(i)%RS(ix, iy, iz) = specie_get_nlcc(s, x)
            end do
          end do
        end do
        call dcf_alloc_FS(ep%rhocore_cf(i))      ! allocate the tube in Fourier space
        call dcf_RS2FS(ep%rhocore_cf(i))         ! Fourier transform
        call dcf_free_RS(ep%rhocore_cf(i))       ! we do not need the real space any longer
      end if nlcc
    
    end do specie
    
    call pop_sub()
  end subroutine epot_local_fourier_init
#endif

  subroutine epot_generate(ep, m, sys, Vpsl, reltype)
    type(epot_type),   intent(inout) :: ep
    type(mesh_type),   intent(in)    :: m
    type(system_type), intent(inout) :: sys
    real(r8),          pointer       :: Vpsl(:)
    integer,           intent(in)    :: reltype

    integer :: ia
    type(specie_type), pointer :: s
    type(atom_type),   pointer :: a
    type(dcf) :: cf_loc, cf_nlcc
    
    call push_sub('generate_external_pot')
    
    ! first we assume that we need to recalculate the ion_ion energy
    sys%eii = ion_ion_energy(sys%natoms, sys%atom)

#ifdef HAVE_FFT
    if(ep%vpsl_space == RECIPROCAL_SPACE) then
      call dcf_new_from(cf_loc, ep%local_cf(1)) ! at least one specie must exist
      call dcf_alloc_FS(cf_loc)
      cf_loc%FS = M_z0

      if(sys%nlcc) then
        call dcf_new_from(cf_nlcc, ep%rhocore_cf(1)) ! at least one specie must exist
        call dcf_alloc_FS(cf_nlcc)
        cf_nlcc%FS = M_z0
      end if
    end if
#endif

    do ia = 1, sys%natoms
      a => sys%atom(ia) ! shortcuts
      s => a%spec
    
      call build_local_part()
      
      if(.not.s%local) then
        call build_kb_sphere()
        call build_nl_part()
      end if
      
    end do

#ifdef HAVE_FFT
    if(ep%vpsl_space == RECIPROCAL_SPACE) then
      ! first the potential
      call dcf_alloc_RS(cf_loc)
      call dcf_FS2RS(cf_loc)
      call dcf2mf(m, cf_loc, Vpsl)
      call dcf_free(cf_loc)
      
      ! and the non-local core corrections
      if(sys%nlcc) then
        call dcf_alloc_RS(cf_nlcc)
        call dcf_FS2RS(cf_nlcc)
        call dcf2mf(m, cf_nlcc, sys%st%rho_core)
        call dcf_free(cf_nlcc)
      end if
    end if
#endif

    if (ep%classic_pot > 0) then
      Vpsl(1:m%np) = Vpsl(1:m%np) + ep%Vclassic(1:m%np)
    end if
    
    call pop_sub()
    
  contains
    subroutine build_local_part()
      integer :: i
      real(r8) :: x(3), r
      
      call push_sub('build_local_part')
      if(ep%vpsl_space == REAL_SPACE) then ! real space
        do i = 1, m%np
          call mesh_xyz(m, i, x)
          x = x - a%x
          Vpsl(i) = Vpsl(i) + specie_get_local(s, x)
          
          if(s%nlcc) then
            sys%st%rho_core(i) = sys%st%rho_core(i) + specie_get_nlcc(s, x)
          end if
        end do
        
#ifdef HAVE_FFT
      else ! momentum space
        call cf_phase_factor(m, a%x, ep%local_cf(s%index), cf_loc)
        if(s%nlcc) then
          call cf_phase_factor(m, a%x, ep%rhocore_cf(s%index), cf_nlcc)
        end if
#endif
      end if
      
      call pop_sub()
    end subroutine build_local_part
  
    subroutine build_kb_sphere()
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
      do k = 1, m%np
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
      do k = 1, m%np
        do i = 1, 3**conf%periodic_dim
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

    subroutine build_nl_part()
      integer :: i, j, k, l, lm, n, add_lm, p, ix(3), center(3)
      real(r8) :: r, x(3), x_in(3), ylm
      real(r8) :: so_uv, so_duv(3)
      
      call push_sub('build_nl_part')
      
      ! This loop is done always, for spin-orbit is, up to now, only done in real space.
      ! If we want the nl part also in real space, it is read.
      j_loop: do j = 1, a%Mps
        call mesh_xyz(m, a%Jxyz(j), x_in)
        k_loop: do k = 1, 3**conf%periodic_dim
          x(:) = x_in(:) - m%shift(k,:)          
          r=sqrt(sum((x-a%x)*(x-a%x)))
          if (r > s%ps%rc_max + m%h(1)) cycle
          x = x - a%x
          add_lm = 1
          l_loop: do l = 0, s%ps%L_max
            lm_loop: do lm = -l, l
              i_loop : do i = 1, s%ps%kbc
                if(l .ne. s%ps%L_loc .and. ep%vnl_space == REAL_SPACE) then
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
      
      if(reltype == 1) then ! SPIN_ORBIT
        do j = 1, a%mps
          call mesh_xyz(m, a%Jxyz(j), x_in)
          do k = 1, 3**conf%periodic_dim
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
    
  end subroutine epot_generate

  subroutine epot_generate_classic(ep, m, sys)
    type(epot_type),   intent(inout) :: ep
    type(mesh_type),   intent(in)    :: m
    type(system_type), intent(IN)    :: sys
    
    integer i, ia
    real(r8) :: r, rc
    
    call push_sub('epot_generate_classic')
    
    ep%Vclassic = M_ZERO
    do ia = 1, sys%ncatoms
      do i = 1, m%np
        call mesh_r(m, i, r, a=sys%catom(ia)%x)
        select case(ep%classic_pot)
        case(1) ! point charge
          if(r < r_small) r = r_small
          ep%Vclassic(i) = ep%Vclassic(i) - sys%catom(ia)%charge/r
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
          ep%Vclassic(i) = ep%Vclassic(i) - sys%catom(ia)%charge*(r**4 - rc**4)/(r**5 - rc**5)
        end select
      end do
    end do

    call pop_sub()
  end subroutine epot_generate_classic

  subroutine epot_laser_field(ep, t, field)
    type(epot_type), intent(in) :: ep
    real(r8),         intent(in)  :: t
    real(r8),         intent(out) :: field(conf%dim)

    call laser_field(ep%no_lasers, ep%lasers, t, field)
  end subroutine epot_laser_field

  subroutine epot_laser_vector_field(ep, t, field)
    type(epot_type), intent(in) :: ep
    real(r8),         intent(in)  :: t
    real(r8),         intent(out) :: field(conf%dim)

    call laser_vector_field(ep%no_lasers, ep%lasers, t, field)
  end subroutine epot_laser_vector_field



#include "undef.F90"
#include "real.F90"
#include "epot_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "epot_inc.F90"

end module external_pot
