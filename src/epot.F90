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
  use global
  use lib_oct_parser
  use mesh
#ifdef HAVE_FFT
  use cube_function
#endif
  use logrid
  use ps
  use specie
  use atom
  use geometry
  use states
  use lasers

  implicit none

  type nonlocal_op
    integer :: n, c
    integer, pointer :: jxyz(:)
    FLOAT,   pointer ::    uv(:, :),    duv(:, :, :),    uvu(:, :)
    CMPLX,   pointer :: so_uv(:, :), so_duv(:, :, :), so_uvu(:, :), so_luv(:, :, :)
    CMPLX,   pointer :: phases(:,:)    ! factors exp(ik*x)
  end type nonlocal_op

  type epot_type
    integer :: vpsl_space           ! How should the local potential be calculated
    integer :: vnl_space            ! How should the nl    potential be calculated
    integer :: nextra               ! extra points for the interpolation method(s)
    
    integer :: classic_pot          ! How to include the classic charges
    FLOAT, pointer :: Vclassic(:)! We use it to store the potential of the classic charges
 
    ! lasers stuff
    integer :: no_lasers ! number of laser pulses used
    type(laser_type), pointer :: lasers(:)

    ! static magnetic field
    FLOAT, pointer :: b(:)

#ifdef HAVE_FFT
    ! For the local pseudopotential in Fourier space...
    type(dcf), pointer :: local_cf(:)
    type(dcf), pointer :: rhocore_cf(:) ! and for the core density
#endif

    ! Nonlocal operators
    integer :: nvnl
    type(nonlocal_op), pointer :: vnl(:)
  end type epot_type

contains

  subroutine epot_init(ep, m, geo)
    type(epot_type),     intent(out) :: ep
    type(mesh_type),     intent(IN)  :: m
    type(geometry_type), intent(IN)  :: geo
    integer :: i

    ! should we calculate the local pseudopotentials in Fourier space?
    ep%vpsl_space = REAL_SPACE
    ep%vnl_space  = REAL_SPACE
#ifdef HAVE_FFT
    call loct_parse_int('LocalPotentialSpace', RECIPROCAL_SPACE, ep%vpsl_space)
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
    if(ep%vpsl_space == RECIPROCAL_SPACE) call epot_local_fourier_init(ep, m, geo)
#endif

    ep%classic_pot = 0
    if(geo%ncatoms > 0) then
      call loct_parse_int("ClassicPotential", 0, ep%classic_pot)
      if(ep%classic_pot > 0) then
        allocate(ep%Vclassic(m%np))
        call epot_generate_classic(ep, m, geo)
      end if
    end if

    ! lasers
    call laser_init(m, ep%no_lasers, ep%lasers)
    if(ep%no_lasers>0 ) then
      message(1) = 'Info: Lasers'
      call write_info(1)
      if(conf%verbose > 20 .and. mpiv%node == 0) then
        call laser_write_info(ep%no_lasers, ep%lasers, stdout)
      end if
    end if

    ! static magnetic field
    nullify(ep%b)
    select case(loct_parse_block_n("StaticMagneticField"))
    case (1)
      allocate(ep%b(conf%dim))
      do i = 1, conf%dim
        call loct_parse_block_float("StaticMagneticField", 0, i-1, ep%b(i))
      end do
    end select

    ! Non local operators
    ep%nvnl = geometry_nvnl(geo)
    nullify(ep%vnl)
    if(ep%nvnl>0) then
       allocate(ep%vnl(ep%nvnl))
       do i = 1, ep%nvnl
          nullify(ep%vnl(i)%jxyz, ep%vnl(i)%uv, ep%vnl(i)%uvu, ep%vnl(i)%duv, &
                  ep%vnl(i)%so_uv, ep%vnl(i)%so_duv, ep%vnl(i)%so_luv, ep%vnl(i)%so_uvu)
       end do
    endif

  end subroutine epot_init

  subroutine epot_end(ep, geo)
    type(epot_type),     intent(inout) :: ep
    type(geometry_type), intent(IN)    :: geo

#ifdef HAVE_FFT
    integer :: i

    if(ep%vpsl_space == RECIPROCAL_SPACE) then
      do i = 1, geo%nspecies
        call dcf_free(ep%local_cf(i))
        if(geo%specie(i)%nlcc) call dcf_free(ep%rhocore_cf(i))
      end do
      deallocate(ep%local_cf)
      if(geo%nlcc) deallocate(ep%rhocore_cf)
    end if
#endif

    if(ep%classic_pot > 0) then
      ep%classic_pot = 0
      ASSERT(associated(ep%Vclassic)) ! sanity check
      deallocate(ep%Vclassic)         ! and clean up
      nullify(ep%Vclassic)
    end if
    
    call laser_end(ep%no_lasers, ep%lasers)

    if(associated(ep%b)) then
      deallocate(ep%b)
    end if

    if(ep%nvnl>0) then
        ASSERT(associated(ep%vnl))
        deallocate(ep%vnl); nullify(ep%vnl)
    endif

  end subroutine epot_end

#ifdef HAVE_FFT
  subroutine epot_local_fourier_init(ep, m, geo)
    type(epot_type),     intent(inout) :: ep
    type(mesh_type),     intent(IN)    :: m
    type(geometry_type), intent(IN)    :: geo
    
    integer :: i, j, ix, iy, iz, ixx(3), db(3), c(3)
    FLOAT :: x(3), modg, temp(3)
    FLOAT, allocatable :: v(:)
    
    type(specie_type), pointer :: s ! shortcuts
    type(dcf), pointer :: cf
    
    call push_sub('epot_local_fourier_init')

    allocate(ep%local_cf(geo%nspecies))
    if(geo%nlcc) allocate(ep%rhocore_cf(geo%nspecies))

    specie: do i = 1, geo%nspecies
      s  => geo%specie(i)
      cf => ep%local_cf(i)
      
      if(i == 1) then
        call mesh_double_box(m, db)
        call dcf_new(db, cf)    ! initialize the cube
        call dcf_fft_init(cf)   ! and initialize the ffts
        db = cf%n               ! dimensions of the box may have been optimized, so get them
        c(:) = db(:)/2 + 1      ! get center of double box
      else
        call dcf_new_from(cf, ep%local_cf(1))   ! we can just copy from the first one
      end if
      
      if(geo%nlcc) call dcf_new_from(ep%rhocore_cf(i), ep%local_cf(1))

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
                       s%ps%g%rofi(j)**2*(loct_splint(s%ps%vlocal,s%ps%g%rofi(j)))
                enddo
                ! this line was not compiling in the intel compiler version 8.
!                cf%FS(ix, iy, iz) = M_FOUR*M_PI*    &
!                     (sum(s%ps%g%drdi(2:s%ps%g%nrval)*v(2:s%ps%g%nrval))-s%ps%z_val/modg)
              else
                cf%FS(ix, iy, iz) = M_ZERO
              end if
            end do
          end do
        end do
        deallocate(v)
        
      end if periodic
    
      ! now we built the non-local core corrections in momentum space
      nlcc: if(s%nlcc) then
        call dcf_alloc_RS(ep%rhocore_cf(i))
        
        do ix = 1, db(1)
          ixx(1) = ix - c(1)
          do iy = 1, db(2)
            ixx(2) = iy - c(2)
            do iz = 1, db(3)
              ixx(3) = iz - c(3)
              
              x(:) = m%h(:)*ixx(:)
              ep%rhocore_cf(i)%RS(ix, iy, iz) = specie_get_nlcc(s, x)
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

  subroutine epot_generate(ep, m, st, geo, Vpsl, reltype)
    type(epot_type),     intent(inout) :: ep
    type(mesh_type),     intent(IN)    :: m
    type(states_type),   intent(inout) :: st
    type(geometry_type), intent(inout) :: geo
    FLOAT,               pointer       :: Vpsl(:)
    integer,             intent(in)    :: reltype

    integer :: ia, i, l, lm, add_lm, k
    type(specie_type), pointer :: s
    type(atom_type),   pointer :: a
    type(dcf) :: cf_loc, cf_nlcc
    
    call push_sub('epot_generate')
    
    ! first we assume that we need to recalculate the ion_ion energy
    geo%eii = ion_ion_energy(geo)

#ifdef HAVE_FFT
    if(ep%vpsl_space == RECIPROCAL_SPACE) then
      call dcf_new_from(cf_loc, ep%local_cf(1)) ! at least one specie must exist
      call dcf_alloc_FS(cf_loc)
      cf_loc%FS = M_z0

      if(geo%nlcc) then
        call dcf_new_from(cf_nlcc, ep%local_cf(1)) ! at least one specie must exist
        call dcf_alloc_FS(cf_nlcc)
        cf_nlcc%FS = M_z0
      end if
    end if
#endif

    ! Local.
    vpsl = M_ZERO
    do ia = 1, geo%natoms
      a => geo%atom(ia) ! shortcuts
      s => a%spec
      call build_local_part()
    enddo
      
    ! Nonlocal part.
    i = 1
    do ia = 1, geo%natoms
      a => geo%atom(ia)
      s => a%spec
      if(s%local) cycle
      add_lm = 1
      do l = 0, s%ps%l_max
         do lm = -l, l
            call nonlocal_op_kill(ep%vnl(i))
            ! This if is a performance hack, necessary for when the ions move.
            ! For each atom, the sphere is the same, so we just calculate it once.
            if(add_lm == 1) then
              k = i
              call build_kb_sphere(i)
            else
              ep%vnl(i)%n = ep%vnl(k)%n
              ep%vnl(i)%c = ep%vnl(k)%c
              ep%vnl(i)%jxyz => ep%vnl(k)%jxyz
            endif
            call allocate_nl_part(i)
            call build_nl_part(i, l, lm, add_lm)
            add_lm = add_lm + 1
            i = i + 1
         enddo
      enddo
    end do

#ifdef HAVE_FFT
    if(ep%vpsl_space == RECIPROCAL_SPACE) then
      ! first the potential
      call dcf_alloc_RS(cf_loc)
      call dcf_FS2RS(cf_loc)
      call dcf2mf(m, cf_loc, Vpsl)
      call dcf_free(cf_loc)
      
      ! and the non-local core corrections
      if(geo%nlcc) then
        call dcf_alloc_RS(cf_nlcc)
        call dcf_FS2RS(cf_nlcc)
        call dcf2mf(m, cf_nlcc, st%rho_core)
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
      FLOAT :: x(3)
      
      call push_sub('build_local_part')
      if(ep%vpsl_space == REAL_SPACE) then ! real space
        do i = 1, m%np
          call mesh_xyz(m, i, x)
          x = x - a%x
          Vpsl(i) = Vpsl(i) + specie_get_local(s, x)
          if(s%nlcc) then
            st%rho_core(i) = st%rho_core(i) + specie_get_nlcc(s, x)
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
  
    ! Note: this should only build the spheres, as it says, not allocating the variables.
    ! I will change it later.
    subroutine build_kb_sphere(ivnl)
      integer, intent(in) :: ivnl
      integer :: i, j, k
      FLOAT :: r
      
      call push_sub('build_kb_sphere')

      if (any(s%ps%rc_max + m%h(1)>=m%lsize(1:conf%periodic_dim))) then
        message(1)='KB sphere is larger than the box size'
        write(message(2),'(a,f12.6,a)')'  rc_max+h = ',s%ps%rc_max + m%h(1),' [b]'
        write(message(3),'(a,3f12.4,a)')'  lsize    = ',m%lsize,' [b]'
        message(4)='Please change pseudopotential'
        call write_fatal(4)
      end if
      j = 0
      do k = 1, m%np
        do i = 1,3**conf%periodic_dim
          call mesh_r(m, k, r, a=a%x+m%shift(i,:))
          if(r > s%ps%rc_max + m%h(1)) cycle
          j = j + 1
          exit
        end do
      end do
      ep%vnl(ivnl)%n = j
      ep%vnl(ivnl)%c = s%ps%kbc

      allocate(ep%vnl(ivnl)%jxyz(j))

      j = 0
      do k = 1, m%np
        do i = 1, 3**conf%periodic_dim
          call mesh_r(m, k, r, a=a%x+m%shift(i,:))
          ! we enlarge slightly the mesh (good for the interpolation scheme)
          if(r > s%ps%rc_max + m%h(1)) cycle
          j = j + 1
          ep%vnl(ivnl)%jxyz(j) = k
          exit
        end do
      end do
      
      call pop_sub()
    end subroutine build_kb_sphere

    subroutine allocate_nl_part(ivnl)
      integer, intent(in) :: ivnl
      integer :: j, c
      call push_sub('allocate_nl_part')

      j = ep%vnl(ivnl)%n
      c = ep%vnl(ivnl)%c

      ! First, the "normal" non-local projectors
      allocate(ep%vnl(ivnl)%uv(j, c), &
               ep%vnl(ivnl)%duv(3, j, c), &
               ep%vnl(ivnl)%uvu(c, c))
      ep%vnl(ivnl)%uv(:,:)     = M_ZERO
      ep%vnl(ivnl)%duv(:,:, :) = M_ZERO
      ep%vnl(ivnl)%uvu(:, :)   = M_ZERO

      ! Then, the spin-orbit projectors (Eventually I will change this).
      allocate(ep%vnl(ivnl)%so_uv(j, c), &
               ep%vnl(ivnl)%so_duv(3, j, c), &
               ep%vnl(ivnl)%so_uvu(s%ps%kbc, c), &
               ep%vnl(ivnl)%so_luv(j, c, 3))
      ep%vnl(ivnl)%so_uv(:, :)     = M_z0
      ep%vnl(ivnl)%so_duv(:, :, :) = M_z0
      ep%vnl(ivnl)%so_uvu(: , :)   = M_z0
      ep%vnl(ivnl)%so_luv(:, :, :) = M_z0

      if(conf%periodic_dim/=0) then
        allocate(ep%vnl(ivnl)%phases(ep%vnl(ivnl)%n, st%d%nik))
      endif

      call pop_sub()
    end subroutine allocate_nl_part

    subroutine build_nl_part(ivnl, l, lm, add_lm)
      integer, intent(in) :: ivnl, l, lm, add_lm

      integer :: i, j, k
      FLOAT :: r, x(3), x_in(3), ylm
      FLOAT :: so_uv, so_duv(3)
      
      call push_sub('build_nl_part')
      
      j_loop: do j = 1, ep%vnl(ivnl)%n
         call mesh_xyz(m, ep%vnl(ivnl)%jxyz(j), x_in)
         k_loop: do k = 1, 3**conf%periodic_dim
          x(:) = x_in(:) - m%shift(k,:)          
          r=sqrt(sum((x-a%x)*(x-a%x)))
          if (r > s%ps%rc_max + m%h(1)) cycle
          x = x - a%x
          i_loop : do i = 1, ep%vnl(ivnl)%c
                if(l .ne. s%ps%L_loc .and. ep%vnl_space == REAL_SPACE) then
                  call specie_get_nl_part(s, x, l, lm, i, ep%vnl(ivnl)%uv(j, i), ep%vnl(ivnl)%duv(:, j, i))
                end if
                if(l>0 .and. s%ps%so_l_max>=0) then
                  call specie_get_nl_part(s, x, l, lm, i, so_uv, so_duv(:), so=.true.)
                  ep%vnl(ivnl)%so_uv(j, i) = so_uv
                  ep%vnl(ivnl)%so_duv(1:3, j, i) = so_duv(1:3)
                endif
          end do i_loop
         end do k_loop
      end do j_loop
      
      if(reltype == 1) then ! SPIN_ORBIT
        do j = 1, ep%vnl(ivnl)%n
          call mesh_xyz(m, ep%vnl(ivnl)%jxyz(j), x_in)
          do k = 1, 3**conf%periodic_dim
            x(:) = x_in(:) - m%shift(k,:)          
            r=sqrt(sum((x-a%x)*(x-a%x)))
            if (r > s%ps%rc_max + m%h(1)) cycle
            x = x - a%x
              ep%vnl(ivnl)%so_luv(j, 1:ep%vnl(ivnl)%c, 1) = &
                   x(2)*ep%vnl(ivnl)%so_duv(3, j, 1:a%spec%ps%kbc) - &
                   x(3)*ep%vnl(ivnl)%so_duv(2, j, 1:a%spec%ps%kbc)
              ep%vnl(ivnl)%so_luv(j, 1:ep%vnl(ivnl)%c, 2) = &
                   x(3)*ep%vnl(ivnl)%so_duv(1, j, 1:a%spec%ps%kbc) - &
                   x(1)*ep%vnl(ivnl)%so_duv(3, j, 1:a%spec%ps%kbc)
              ep%vnl(ivnl)%so_luv(j, 1:ep%vnl(ivnl)%c, 3) = &
                   x(1)*ep%vnl(ivnl)%so_duv(2, j, 1:a%spec%ps%kbc ) - &
                   x(2)*ep%vnl(ivnl)%so_duv(1, j, 1:a%spec%ps%kbc )
            exit
          enddo
        enddo
        ep%vnl(ivnl)%so_luv = -M_zI*ep%vnl(ivnl)%so_luv
      endif
      
      ! and here we calculate the uVu
      if(s%ps%flavour(1:2) == 'tm') then
        ep%vnl(ivnl)%uvu(:, :) = M_ZERO
            if(l .ne. s%ps%L_loc) then
              do j = 1, ep%vnl(ivnl)%n
                call mesh_r(m, ep%vnl(ivnl)%jxyz(j), r, x=x_in, a=a%x)
                do k=1,3**conf%periodic_dim
                  x(:) = x_in(:) - m%shift(k,:)          
                  r=sqrt(sum(x*x))
                  if (r > s%ps%rc_max + m%h(1)) cycle
                  ylm = loct_ylm(x(1), x(2), x(3), l, lm)
                  ep%vnl(ivnl)%uvu(1, 1) = ep%vnl(ivnl)%uvu(1, 1) + ep%vnl(ivnl)%uv(j, 1) * &
                                           loct_splint(s%ps%ur(l+1, 1), r)*ylm*m%vol_pp
                  exit
                end do
              end do
              ep%vnl(ivnl)%uvu(1, 1) = M_ONE/(ep%vnl(ivnl)%uvu(1, 1)*s%ps%dknrm(l))
!!$ This check is temporarily gone. I will put it back.
!!$              if(abs((a%duVu(add_lm, 1, 1) - s%ps%h(l,1,1))/s%ps%h(l,1,1)) > CNST(0.05) .and. s%ps%ispin==1) then
!!$                write(message(1), '(a,i4)') "Low precision in the calculation of the uVu for lm = ", &
!!$                     add_lm
!!$                write(message(2), '(f14.6,a,f14.6)') s%ps%h(l,1,1), ' .ne. ', a%duVu(add_lm, 1, 1)
!!$                message(3) = "Please consider decreasing the spacing, or changing pseudopotential"
!!$                call write_warning(3)
!!$              end if
!!$              ! uVu is in fact calculated in the logarithmic grid, previous stuff is a check.
              ep%vnl(ivnl)%uvu(1, 1) = s%ps%h(l, 1, 1)
            end if
        ep%vnl(ivnl)%so_uvu(1, 1) = s%ps%k(l, 1, 1)
      else
        ep%vnl(ivnl)%uvu(1:3, 1:3) = s%ps%h(l, 1:3, 1:3)
        ep%vnl(ivnl)%so_uvu(1:3, 1:3) = s%ps%k(l, 1:3, 1:3)
      end if

      if(conf%periodic_dim/=0) then
        do j = 1, ep%vnl(ivnl)%n
           call mesh_xyz(m, ep%vnl(ivnl)%jxyz(j), x)
           do k=1, st%d%nik
              ep%vnl(ivnl)%phases(j, k) = exp(M_zI*sum(st%d%kpoints(:, k)*x(:)))
           end do
        end do
      endif

      call pop_sub(); return
    end subroutine build_nl_part
    
  end subroutine epot_generate

  subroutine epot_generate_classic(ep, m, geo)
    type(epot_type),     intent(inout) :: ep
    type(mesh_type),     intent(IN)    :: m
    type(geometry_type), intent(IN)    :: geo
    
    integer i, ia
    FLOAT :: r, rc
    
    call push_sub('epot_generate_classic')
    
    ep%Vclassic = M_ZERO
    do ia = 1, geo%ncatoms
      do i = 1, m%np
        call mesh_r(m, i, r, a=geo%catom(ia)%x)
        select case(ep%classic_pot)
        case(1) ! point charge
          if(r < r_small) r = r_small
          ep%Vclassic(i) = ep%Vclassic(i) - geo%catom(ia)%charge/r
        case(2) ! gaussion smeared charge
          select case(geo%catom(ia)%label(1:1)) ! covalent radii
          case('H')
            rc = CNST(0.4)*P_Ang
          case('C') 
            rc = CNST(0.8)*P_Ang
          case default
            rc = CNST(0.7)*P_Ang
          end select
          if(abs(r - rc) < r_small) r = rc + sign(r_small, r-rc)
          ep%Vclassic(i) = ep%Vclassic(i) - geo%catom(ia)%charge*(r**4 - rc**4)/(r**5 - rc**5)
        end select
      end do
    end do

    call pop_sub()
  end subroutine epot_generate_classic

  subroutine epot_laser_field(ep, t, field)
    type(epot_type), intent(IN) :: ep
    FLOAT,           intent(in)  :: t
    FLOAT,           intent(out) :: field(conf%dim)

    call laser_field(ep%no_lasers, ep%lasers, t, field)
  end subroutine epot_laser_field

  subroutine epot_laser_vector_field(ep, t, field)
    type(epot_type), intent(IN) :: ep
    FLOAT,           intent(in)  :: t
    FLOAT,           intent(out) :: field(conf%dim)

    call laser_vector_field(ep%no_lasers, ep%lasers, t, field)
  end subroutine epot_laser_vector_field

  subroutine nonlocal_op_kill(nlop)
    type(nonlocal_op) :: nlop
    if(associated(nlop%jxyz)) then
      deallocate(nlop%jxyz, nlop%uv, nlop%uvu, nlop%duv, &
                 nlop%so_uv, nlop%so_duv, nlop%so_luv, nlop%so_uvu)
      if(conf%periodic_dim/=0) deallocate(nlop%phases)
    endif
  end subroutine nonlocal_op_kill

#include "undef.F90"
#include "real.F90"
#include "epot_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "epot_inc.F90"

end module external_pot
