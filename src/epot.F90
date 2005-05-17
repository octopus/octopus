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
!!
!! $Id$

#include "global.h"

module external_pot
  use global
  use units
  use lib_oct_parser
  use lib_oct_gsl_spline
  use fft
  use math
  use mesh
  use functions
#ifdef HAVE_FFT
  use cube_function
#endif
  use logrid
  use ps
  use specie
  use geometry
  use states
  use lasers

  implicit none

  private
  public :: nonlocal_op,            &
            epot_type,              &
            epot_init,              &
            epot_end,               &
            epot_generate,          &
            epot_laser_scalar_pot,  &
            epot_laser_field,       &
            epot_laser_vector_pot,  &
            depot_forces, zepot_forces

  type nonlocal_op
    integer :: n, c
    integer, pointer :: jxyz(:)
    FLOAT,   pointer ::    uv(:, :),    duv(:, :, :),    uvu(:, :)
    CMPLX,   pointer :: so_uv(:, :), so_duv(:, :, :), so_uvu(:, :), so_luv(:, :, :)
    CMPLX,   pointer :: phases(:,:)    ! factors exp(ik*x)
  end type nonlocal_op

  type epot_type
    !Classic charges:
    integer :: classic_pot        ! How to include the classic charges
    FLOAT, pointer :: vclassic(:) ! We use it to store the potential of the classic charges
 
    !Ions
    FLOAT, pointer :: vpsl(:)            ! the local part of the pseudopotentials
#ifdef HAVE_FFT
    type(dcf), pointer :: local_cf(:)    ! for the local pseudopotential in Fourier space
    type(dcf), pointer :: rhocore_cf(:)  ! and for the core density
#endif
    integer :: nvnl                      ! number of nonlocal operators
    type(nonlocal_op), pointer :: vnl(:) ! nonlocal operators

    !External e-m fields
    integer :: no_lasers                   ! number of laser pulses used
    type(laser_type), pointer :: lasers(:) ! lasers stuff
    FLOAT, pointer :: b(:)                 ! static magnetic field
    FLOAT, pointer :: a(:,:)               ! static vector potential

  end type epot_type

contains

  subroutine epot_init(ep, m, ispin, geo)
    type(epot_type),     intent(out) :: ep
    type(mesh_type),     intent(in)  :: m
    integer,             intent(in)  :: ispin
    type(geometry_type), intent(in)  :: geo

    integer :: i
    integer(POINTER_SIZE) :: blk
    FLOAT, allocatable :: x(:)

    call push_sub('epot_init')

    !Local part of the pseudopotentials
    allocate(ep%vpsl(m%np))
    ep%vpsl = M_ZERO

    ! should we calculate the local pseudopotentials in Fourier space?
    ! This depends on wether we have periodic dimensions or not
    if(conf%periodic_dim>0) call epot_local_fourier_init(ep, m, geo)

    ep%classic_pot = 0
    if(geo%ncatoms > 0) then
      call loct_parse_int(check_inp('ClassicPotential'), 0, ep%classic_pot)
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

    ! static magnetic fields
    nullify(ep%b, ep%a)
    if(loct_parse_block(check_inp('StaticMagneticField'), blk)==0) then
      select case(loct_parse_block_n(blk))
      case (1)
        select case (conf%dim)
        case (1)
        case (2)
          if (loct_parse_block_cols(blk, 0) /= 1) then
            message(1) = "When performing 2D calculations, the external magnetic field can only have one component,"
            message(2) = "corresponding to the direction perpendicular to the plane."
            call write_fatal(2)
          else
            allocate(ep%b(1))
            call loct_parse_block_float(blk, 0, 0, ep%b(1))
          end if
        case (3)
          allocate(ep%b(3))
          do i = 1, 3
            call loct_parse_block_float(blk, 0, i-1, ep%b(i))
          end do
        end select
      end select
      call loct_parse_block_end(blk)

      !Compute the vector potential
      allocate(ep%a(m%np, conf%dim), x(conf%dim))
      do i = 1, m%np
        call mesh_xyz(m, i, x)
        select case (conf%dim)
        case (2)
          ep%a(i, :) = -M_HALF*(/x(2)*ep%b(1), x(1)*ep%b(1)/)
        case (3)
          ep%a(i, :) = -M_HALF*(/x(2)*ep%b(3) - x(3)*ep%b(2), x(3)*ep%b(1) - x(1)*ep%b(3), x(1)*ep%b(2) - x(2)*ep%b(1)/)
        end select
      end do
      deallocate(x)

    end if

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

    call pop_sub()
  end subroutine epot_init

  subroutine epot_end(ep, geo)
    type(epot_type),     intent(inout) :: ep
    type(geometry_type), intent(IN)    :: geo

#ifdef HAVE_FFT
    integer :: i

    if(conf%periodic_dim>0) then
      do i = 1, geo%nspecies
        call dcf_free(ep%local_cf(i))
        if(geo%specie(i)%nlcc) call dcf_free(ep%rhocore_cf(i))
      end do
      deallocate(ep%local_cf)
      if(geo%nlcc) deallocate(ep%rhocore_cf)
    endif

#endif

    if(associated(ep%vpsl)) then
      deallocate(ep%vpsl)
      nullify(ep%vpsl)
    end if

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
    if(associated(ep%a)) then
      deallocate(ep%a)
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
    
    integer :: vlocal_cutoff
    integer :: i, j, ix, iy, iz, ixx(3), db(3), c(3)
    FLOAT :: x(3)
    FLOAT :: g(3), gpar, gperp, gx, gz, modg
    FLOAT :: a_erf, r_0, temp(3), tmp, norm
    
    type(specie_type), pointer :: s ! shortcuts
    type(dcf), pointer :: cf
    
    call push_sub('epot_local_fourier_init')

    call loct_parse_int(check_inp('VlocalCutoff'), conf%periodic_dim , vlocal_cutoff)
    if (vlocal_cutoff /= conf%periodic_dim) then
      write(message(1), '(a,i1,a)')'The System is periodic in ',conf%periodic_dim ,' dimension(s),'
      write(message(2), '(a,i1,a)')'but VlocalCutoff is set for ',vlocal_cutoff,' dimensions.'
      call write_warning(2)
    end if
    select case(vlocal_cutoff)
      case(0)
        message(1) = 'Info: Vlocal Cutoff = sphere'
        call write_info(1)
      case(1)
        message(1) = 'Info: Vlocal Cutoff = cylinder'
        call write_info(1)
      case(2)
        message(1) = 'Info: Vlocal Cutoff = slab'
        call write_info(1)
      case(3)
        message(1) = 'Info: Vlocal Cutoff = no cutoff'
        call write_info(1)
      case default
        write(message(1), '(a,i2,a)') 'Input: ', vlocal_cutoff, &
           ' is not a valid VlocalCutoff'
        message(2) = 'VlocalCutoff = 0(spherical) | 1(cylindrical) | 2(planar)'
        message(3) = '               3(no cutoff)'
        call write_fatal(3)
    end select

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
        if (vlocal_cutoff == 3) then
          r_0 = M_ZERO
        else
          call loct_parse_float(check_inp('VlocalCutoffRadius'),&
               maxval(db(:)*m%h(:)/M_TWO)/units_inp%length%factor , r_0)
          r_0 = r_0*units_inp%length%factor
          write(message(1),'(3a,f12.6)')'Info: Vlocal Cutoff Radius [',  &
                            trim(units_out%length%abbrev), '] = ',       &
                            r_0/units_out%length%factor
          call write_info(1)
        end if
      else
        call dcf_new_from(cf, ep%local_cf(1))   ! we can just copy from the first one
      end if
      
      if(geo%nlcc) call dcf_new_from(ep%rhocore_cf(i), ep%local_cf(1))

      call dcf_alloc_FS(cf)      ! allocate the tube in Fourier space

        a_erf = M_TWO
        norm = M_FOUR*M_PI/m%vol_pp(1)
        temp(:) = M_TWO*M_PI/(db(:)*m%h(:))
        cf%FS = M_Z0
        do ix = 1, cf%nx
          ixx(1) = pad_feq(ix, db(1), .true.)
          do iy = 1, db(2)
            ixx(2) = pad_feq(iy, db(2), .true.)
            do iz = 1, db(3)
              ixx(3) = pad_feq(iz, db(3), .true.)

              x = temp(:)*ixx(:)
              modg = sqrt(sum((temp(:)*ixx(:))**2))

              tmp = specie_get_local_fourier(s, x)
              if(modg /= M_ZERO) then
                tmp = tmp - s%z_val*exp(-(modg/(2*a_erf))**2)/modg**2
                select case(vlocal_cutoff)
                case(0) 
                  cf%FS(ix, iy, iz) = tmp*cutoff0(modg*r_0)
                case(1)
                  gx = abs(temp(1)*ixx(1))
                  gperp = sqrt((temp(2)*ixx(2))**2+(temp(3)*ixx(3))**2)
                  cf%FS(ix, iy, iz) = tmp*cutoff1(gx*r_0, gperp*r_0)
                case(2)
                  gz = abs(temp(3)*ixx(3))
                  gpar = sqrt((temp(1)*ixx(1))**2+(temp(2)*ixx(2))**2)
                  cf%FS(ix, iy, iz) = tmp*cutoff2(gpar*r_0, gz*r_0)
                case(3)
                  cf%FS(ix, iy, iz) = tmp
                end select
              else
                select case(vlocal_cutoff)
                case(0)  ; cf%FS(ix, iy, iz) = -r_0**2/M_TWO
                case(1,2); cf%FS(ix, iy, iz) = M_ZERO
                case(3)  ; cf%FS(ix, iy, iz) = tmp
                end select
              end if
              
 ! multiply by normalization factor and a phase shift to get the center of the box
              cf%FS(ix, iy, iz) = norm*                         &
                                  exp(M_PI*M_ZI*sum(ixx(:)))*   &
                                  cf%FS(ix, iy, iz)
 
            end do
          end do
        end do
        
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

  subroutine epot_generate(ep, m, st, geo, reltype, fast_generation)
    type(epot_type),     intent(inout) :: ep
    type(mesh_type),     intent(IN)    :: m
    type(states_type),   intent(inout) :: st
    type(geometry_type), intent(inout) :: geo
    integer,             intent(in)    :: reltype
    logical, optional, intent(in)      :: fast_generation

    logical :: fast_generation_
    integer :: ia, i, l, lm, add_lm, k, ierr, p
    type(specie_type), pointer :: s
    type(atom_type),   pointer :: a
    type(dcf) :: cf_loc, cf_nlcc
    
    call push_sub('epot_generate')

    fast_generation_ = .false.
    if (present(fast_generation)) fast_generation_ = fast_generation
    
    ! first we assume that we need to recalculate the ion_ion energy
    geo%eii = ion_ion_energy(geo)

#ifdef HAVE_FFT
    if(conf%periodic_dim>0) then
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
    ep%vpsl = M_ZERO
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
      p = 1
      do l = 0, s%ps%l_max
         if(s%ps%l_loc == l) then
            add_lm = add_lm + 2*l + 1
            cycle
         endif
         do lm = -l, l
            if(.not.fast_generation_) then
               call nonlocal_op_kill(ep%vnl(i))
               ! This if is a performance hack, necessary for when the ions move.
               ! For each atom, the sphere is the same, so we just calculate it once.
               if(p == 1) then
                 k = i
                 p = 2
                 deallocate(ep%vnl(i)%jxyz, stat = ierr)
                 call build_kb_sphere(i)
               else
                 ep%vnl(i)%n = ep%vnl(k)%n
                 ep%vnl(i)%c = ep%vnl(k)%c
                 ep%vnl(i)%jxyz => ep%vnl(k)%jxyz
               endif
               call allocate_nl_part(i)
            endif
            call build_nl_part(i, l, lm, add_lm)
            add_lm = add_lm + 1
            i = i + 1
         enddo
      enddo
    end do

#ifdef HAVE_FFT
    if(conf%periodic_dim>0) then
      ! first the potential
      call dcf_alloc_RS(cf_loc)
      call dcf_FS2RS(cf_loc)
      call dcf2mf(m, cf_loc, ep%vpsl)
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
      ep%vpsl(1:m%np) = ep%vpsl(1:m%np) + ep%vclassic(1:m%np)
    end if
    
    call pop_sub()
    
  contains
    subroutine build_local_part()
      integer :: i
      FLOAT :: x(3)
      
      call push_sub('build_local_part')
      if(conf%periodic_dim==0) then
        do i = 1, m%np
          call mesh_xyz(m, i, x)
          x = x - a%x
          ep%vpsl(i) = ep%vpsl(i) + specie_get_local(s, x)
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
                if(l .ne. s%ps%L_loc) then
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
        ep%vnl(ivnl)%so_luv(:, :, :) = -M_zI*ep%vnl(ivnl)%so_luv(:, :, :)
      endif
      
      ! and here we calculate the uVu
      if(s%ps%flavour == PS_TM2) then
        ep%vnl(ivnl)%uvu(:, :) = M_ZERO
            if(l .ne. s%ps%L_loc) then
              do j = 1, ep%vnl(ivnl)%n
                call mesh_r(m, ep%vnl(ivnl)%jxyz(j), r, x=x_in, a=a%x)

                do k = 1, 3**conf%periodic_dim
                  x(:) = x_in(:) - m%shift(k,:)          
                  r    = sqrt(sum(x*x))
                  if (r > s%ps%rc_max + m%h(1)) cycle

                  ylm = loct_ylm(x(1), x(2), x(3), l, lm)
                  ep%vnl(ivnl)%uvu(1, 1) = ep%vnl(ivnl)%uvu(1, 1) + ep%vnl(ivnl)%uv(j, 1) * &
                     loct_splint(s%ps%ur(l+1, 1), r) * ylm * &
                     m%vol_pp(ep%vnl(ivnl)%jxyz(j))
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

  subroutine epot_laser_scalar_pot(ep, x, t, v)
    type(epot_type), intent(in)  :: ep
    FLOAT,           intent(in)  :: t
    FLOAT,           intent(in)  :: x(conf%dim)
    FLOAT,           intent(out) :: v

    FLOAT, allocatable :: e(:)

    allocate(e(conf%dim))
    call laser_field(ep%no_lasers, ep%lasers, t, e)
    v = sum(e*x)
    deallocate(e)

  end subroutine epot_laser_scalar_pot

  subroutine epot_laser_vector_pot(ep, x, t, a)
    type(epot_type), intent(in)  :: ep
    FLOAT,           intent(in)  :: t
    FLOAT,           intent(in)  :: x(conf%dim)
    FLOAT,           intent(out) :: a(conf%dim)

    call laser_vector_field(ep%no_lasers, ep%lasers, t, a)

  end subroutine epot_laser_vector_pot

  subroutine epot_laser_field(ep, t, e)
    type(epot_type), intent(in)  :: ep
    FLOAT,           intent(in)  :: t
    FLOAT,           intent(out) :: e(conf%dim)

    call laser_field(ep%no_lasers, ep%lasers, t, e)

  end subroutine epot_laser_field

  ! This subroutine deallocates (if associated) every component of a nonlocal operator
  ! type variabel, except for jxyz. This is not very elegant, but it helps with performance.
  subroutine nonlocal_op_kill(nlop)
    type(nonlocal_op) :: nlop
    if(associated(nlop%uv)) then
      deallocate(nlop%uv, nlop%uvu, nlop%duv, &
                 nlop%so_uv, nlop%so_duv, nlop%so_luv, nlop%so_uvu)
      if(conf%periodic_dim/=0 .and. associated(nlop%phases)) then
        deallocate(nlop%phases)
      end if
    end if
  end subroutine nonlocal_op_kill

#include "undef.F90"
#include "real.F90"
#include "epot_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "epot_inc.F90"

end module external_pot


