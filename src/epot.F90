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
  use messages
  use units
  use lib_oct
  use lib_oct_parser
  use lib_oct_gsl_spline
  use fft
  use math
  use mesh
  use grid
  use simul_box
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
  public :: epot_type,              &
            epot_init,              &
            epot_end,               &
            epot_generate,          &
            epot_laser_scalar_pot,  &
            epot_laser_field,       &
            epot_laser_vector_pot,  &
            depot_forces,           &
            zepot_forces,           &
            dproject, zproject,     &
            projector

  ! /* The projector data type is intended to hold the non-local part of the
  ! pseudopotentials. The definition of the action of a projector (which is
  ! done through the X(project) subroutine) is:
  ! \hat{P} |phi> = \sum_{ij} uvu(i, j) | a_i >< b_j |phi>
  ! i and j run from one to c. The phases are applied whenever the system is
  ! periodic in some dimension.
  ! For the pseudopotentials themselves, a = b. But to calculate the forces,
  ! one needs to put in b the gradient of a. And in the case of the spin-orbit
  ! coupling term, I have put in b the angular momentum of a. */
  type projector
    integer :: n, c
    integer, pointer :: jxyz(:)
    FLOAT, pointer   :: uvu(:, :)
    FLOAT, pointer   :: a(:, :), b(:, :)
    CMPLX, pointer   :: phases(:, :)
    integer :: index
  end type projector

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
    type(projector), pointer :: p(:)
    type(projector), pointer :: dp(:, :)
    type(projector), pointer :: lso(:, :)

    !External e-m fields
    integer :: no_lasers                   ! number of laser pulses used
    type(laser_type), pointer :: lasers(:) ! lasers stuff
    FLOAT, pointer :: e(:)                 ! static electric field
    FLOAT, pointer :: v(:)                 ! static scalar potential
    FLOAT, pointer :: b(:)                 ! static magnetic field
    FLOAT, pointer :: a(:,:)               ! static vector potential

  end type epot_type

contains

  subroutine epot_init(ep, gr)
    type(epot_type),     intent(out) :: ep
    type(grid_type),     intent(in)  :: gr

    integer :: i, j
    integer(POINTER_SIZE) :: blk
    FLOAT, allocatable :: x(:)

    call push_sub('epot.epot_init')

    !Local part of the pseudopotentials
    allocate(ep%vpsl(NP))
    ep%vpsl = M_ZERO

    ! should we calculate the local pseudopotentials in Fourier space?
    ! This depends on wether we have periodic dimensions or not
    if(simul_box_is_periodic(gr%sb).and.(.not.gr%geo%only_user_def)) then
      call epot_local_fourier_init(ep, gr%m, gr%sb, gr%geo)
    end if

    ep%classic_pot = 0
    if(gr%geo%ncatoms > 0) then
      call loct_parse_int(check_inp('ClassicPotential'), 0, ep%classic_pot)
      if(ep%classic_pot > 0) then
        message(1) = 'Info: generating classic external potential'
        call write_info(1)

        allocate(ep%Vclassic(NP))
        call epot_generate_classic(ep, gr%m, gr%geo)
      end if
    end if

    ! lasers
    call laser_init(ep%no_lasers, ep%lasers, gr%m)
    if(ep%no_lasers>0 ) then
      message(1) = 'Info: Lasers'
      call write_info(1)
      if(mpiv%node == 0) then
        call laser_write_info(ep%no_lasers, ep%lasers, stdout)
      end if
    end if

    ! static electric field
    nullify(ep%e, ep%v)
    if(loct_parse_block(check_inp('StaticElectricField'), blk)==0) then
      select case(loct_parse_block_n(blk))
      case (1)
        do i = 1, NDIM
          call loct_parse_block_float(blk, 0, i-1, ep%e(i))
        end do
      end select
      call loct_parse_block_end(blk)

      !Compute the scalar potential
      allocate(ep%v(NP), x(NDIM))
      do i = 1, NP
        ep%v(i) = sum(gr%m%x(i,:)*ep%e(:))
      end do
      deallocate(x)
    end if

    ! static magnetic field
    nullify(ep%b, ep%a)
    if(loct_parse_block(check_inp('StaticMagneticField'), blk)==0) then
      select case(loct_parse_block_n(blk))
      case (1)
        select case (NDIM)
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
      allocate(ep%a(NP, NDIM), x(NDIM))
      do i = 1, NP
        x = gr%m%x(i, :)
        select case (NDIM)
        case (2)
          ep%a(i, :) = -M_HALF*(/x(2)*ep%b(1), x(1)*ep%b(1)/)
        case (3)
          ep%a(i, :) = -M_HALF*(/x(2)*ep%b(3) - x(3)*ep%b(2), x(3)*ep%b(1) - x(1)*ep%b(3), x(1)*ep%b(2) - x(2)*ep%b(1)/)
        end select
      end do
      deallocate(x)

    end if

    ! The projectors
    ep%nvnl = geometry_nvnl(gr%geo)
    nullify(ep%p)
    if(ep%nvnl > 0) then
       allocate(ep%p(ep%nvnl))
       do i = 1, ep%nvnl
          nullify(ep%p(i)%jxyz, ep%p(i)%a, ep%p(i)%b, ep%p(i)%uvu, ep%p(i)%phases)
       enddo
       allocate(ep%dp(NDIM, ep%nvnl))
       do i = 1, ep%nvnl
          do j = 1, NDIM
             nullify(ep%dp(j, i)%jxyz, ep%dp(j, i)%a, ep%dp(j, i)%b, ep%dp(j, i)%uvu, ep%dp(j, i)%phases)
          enddo
       enddo
       allocate(ep%lso(NDIM, ep%nvnl))
       do i = 1, ep%nvnl
          do j = 1, NDIM
             nullify(ep%lso(j, i)%jxyz, ep%lso(j, i)%a, ep%lso(j, i)%b, ep%lso(j, i)%uvu, ep%lso(j, i)%phases)
          enddo
       enddo
    endif

    call pop_sub()
  end subroutine epot_init

  subroutine epot_end(ep, sb, geo)
    type(epot_type),      intent(inout) :: ep
    type(simul_box_type), intent(in)    :: sb
    type(geometry_type),  intent(in)    :: geo

#ifdef HAVE_FFT
    integer :: i

    if(simul_box_is_periodic(sb).and.(.not.geo%only_user_def)) then
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
    if(associated(ep%e)) then
      deallocate(ep%e)
    end if
    if(associated(ep%v)) then
      deallocate(ep%v)
    end if
    if(associated(ep%b)) then
      deallocate(ep%b)
    end if
    if(associated(ep%a)) then
      deallocate(ep%a)
    end if

    if(ep%nvnl>0) then
        ASSERT(associated(ep%p))
        deallocate(ep%p); nullify(ep%p)

        ASSERT(associated(ep%dp))
        deallocate(ep%dp); nullify(ep%dp)

        ! Here the spin-orbit should be finalized, but only if they have been built...
    endif

  end subroutine epot_end

#ifdef HAVE_FFT
  subroutine epot_local_fourier_init(ep, m, sb, geo)
    type(epot_type),     intent(inout) :: ep
    type(mesh_type),     intent(in)    :: m
    type(simul_box_type), intent(in)   :: sb
    type(geometry_type), intent(in)    :: geo

    integer :: vlocal_cutoff
    integer :: i, ix, iy, iz, ixx(3), db(3), c(3)
    FLOAT :: x(3)
    FLOAT :: gpar, gperp, gx, gz, modg
    FLOAT :: r_0, temp(3), tmp, norm

    type(specie_type), pointer :: s ! shortcuts
    type(dcf), pointer :: cf

    call push_sub('epot.epot_local_fourier_init')

    call loct_parse_int(check_inp('VlocalCutoff'), sb%periodic_dim , vlocal_cutoff)
    if (vlocal_cutoff /= sb%periodic_dim) then
      write(message(1), '(a,i1,a)')'The System is periodic in ', sb%periodic_dim, ' dimension(s),'
      write(message(2), '(a,i1,a)')'but VlocalCutoff is set for ', vlocal_cutoff, ' dimensions.'
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
        call mesh_double_box(sb, m, db)
        call dcf_new(db, cf)    ! initialize the cube
        call dcf_fft_init(cf, sb)   ! and initialize the ffts
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

              tmp = specie_get_local_fourier(sb%dim, s, x)
              if(modg /= M_ZERO) then
                tmp = tmp - s%z_val*exp(-(modg/(2*s%ps%a_erf))**2)/modg**2
                select case(vlocal_cutoff)
                case(0)
                  cf%FS(ix, iy, iz) = tmp*cutoff0(modg,r_0)
                case(1)
                  gx = abs(temp(1)*ixx(1))
                  gperp = sqrt((temp(2)*ixx(2))**2+(temp(3)*ixx(3))**2)
                  cf%FS(ix, iy, iz) = tmp*cutoff1(gx, gperp,r_0)
                case(2)
                  gz = abs(temp(3)*ixx(3))
                  gpar = sqrt((temp(1)*ixx(1))**2+(temp(2)*ixx(2))**2)
                  cf%FS(ix, iy, iz) = tmp*cutoff2(gpar, gz,r_0)
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

  subroutine epot_generate(ep, m, sb, geo, st, reltype, fast_generation)
    type(epot_type),      intent(inout) :: ep
    type(mesh_type),      intent(in)    :: m
    type(simul_box_type), intent(in)    :: sb
    type(geometry_type),  intent(inout) :: geo
    type(states_type),    intent(inout) :: st
    integer,              intent(in)    :: reltype
    logical, optional,    intent(in)    :: fast_generation

    logical :: fast_generation_
    integer :: ia, i, l, lm, add_lm, k, ierr, p
    type(specie_type), pointer :: s
    type(atom_type),   pointer :: a
    type(dcf) :: cf_loc, cf_nlcc

    integer :: j

    call push_sub('epot.epot_generate')

    fast_generation_ = .false.
    if (present(fast_generation)) fast_generation_ = fast_generation

    ! first we assume that we need to recalculate the ion_ion energy
    geo%eii = ion_ion_energy(geo)

#ifdef HAVE_FFT
    if(simul_box_is_periodic(sb).and.(.not.geo%only_user_def)) then
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
               !The projectors maybe should be killed, in the same way that the nonlocal_op where.
               !call nonlocal_op_kill(ep%vnl(i), sb)
               ! This if is a performance hack, necessary for when the ions move.
               ! For each atom, the sphere is the same, so we just calculate it once.
               if(p == 1) then
                 k = i
                 p = 2
                 deallocate(ep%p(i)%jxyz, stat = ierr)
                 do j = 1, sb%dim
                    deallocate(ep%dp(j, i)%jxyz, stat = ierr)
                 enddo
                 do j = 1, sb%dim
                    deallocate(ep%lso(j, i)%jxyz, stat = ierr)
                 enddo
                 call build_kb_sphere(i)
               else
                 ep%p(i)%n = ep%p(k)%n
                 ep%p(i)%c = ep%p(k)%c
                 ep%p(i)%jxyz => ep%p(k)%jxyz
                 do j = 1, sb%dim
                    ep%dp(j, i)%n = ep%dp(j, k)%n
                    ep%dp(j, i)%c = ep%dp(j, k)%c
                    ep%dp(j, i)%jxyz => ep%dp(j, k)%jxyz
                 enddo
                 do j = 1, sb%dim
                    ep%lso(j, i)%n = ep%lso(j, k)%n
                    ep%lso(j, i)%c = ep%lso(j, k)%c
                    ep%lso(j, i)%jxyz => ep%lso(j, k)%jxyz
                 enddo
               endif
               call allocate_nl_part(i)
            endif
            call build_nl_part(i, l, lm)
            add_lm = add_lm + 1

            ep%p(i)%index = ia
            ep%dp(1:sb%dim, i)%index = ia
            ep%lso(1:sb%dim, i)%index = ia

            i = i + 1
         enddo
      enddo
    end do

#ifdef HAVE_FFT
    if(simul_box_is_periodic(sb).and.(.not.geo%only_user_def)) then
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

      call push_sub('epot.build_local_part')

      if((.not.simul_box_is_periodic(sb)).or.geo%only_user_def) then
        do i = 1, m%np
          x(:) = m%x(i, :) - a%x(:)
          ep%vpsl(i) = ep%vpsl(i) + specie_get_local(s, x)
          if(s%nlcc) then
            st%rho_core(i) = st%rho_core(i) + specie_get_nlcc(s, x)
          end if
        end do
#ifdef HAVE_FFT
      else ! momentum space
        call cf_phase_factor(sb, m, a%x, ep%local_cf(s%index), cf_loc)
        if(s%nlcc) then
          call cf_phase_factor(sb, m, a%x, ep%rhocore_cf(s%index), cf_nlcc)
        end if
#endif
      end if

      call pop_sub()
    end subroutine build_local_part

    subroutine build_kb_sphere(ivnl)
      integer, intent(in) :: ivnl
      integer :: i, j, k, d
      FLOAT :: r

      call push_sub('epot.build_kb_sphere')

      if (any(s%ps%rc_max + m%h(1) >= sb%lsize(1:sb%periodic_dim))) then
        message(1)='KB sphere is larger than the box size'
        write(message(2),'(a,f12.6,a)')  '  rc_max+h = ', s%ps%rc_max + m%h(1), ' [b]'
        write(message(3),'(a,3f12.4,a)') '  lsize    = ', sb%lsize, ' [b]'
        message(4)='Please change pseudopotential'
        call write_fatal(4)
      end if
      j = 0
      do k = 1, m%np
        do i = 1, 3**sb%periodic_dim
          call mesh_r(m, k, r, a=a%x + sb%shift(i,:))
          if(r > s%ps%rc_max + m%h(1)) cycle
          j = j + 1
          exit
        end do
      end do

      ep%p(ivnl)%n = j
      ep%p(ivnl)%c = s%ps%kbc

      do d = 1, sb%dim
         ep%dp(d, ivnl)%n = j
         ep%dp(d, ivnl)%c = s%ps%kbc
      enddo

      do d = 1, sb%dim
         ep%lso(d, ivnl)%n = j
         ep%lso(d, ivnl)%c = s%ps%kbc
      enddo

      allocate(ep%p(ivnl)%jxyz(j))
      do d = 1, sb%dim
         allocate(ep%dp(d, ivnl)%jxyz(j))
      enddo
      do d = 1, sb%dim
         allocate(ep%lso(d, ivnl)%jxyz(j))
      enddo

      j = 0
      do k = 1, m%np
        do i = 1, 3**sb%periodic_dim
          call mesh_r(m, k, r, a=a%x + sb%shift(i,:))
          ! we enlarge slightly the mesh (good for the interpolation scheme)
          if(r > s%ps%rc_max + m%h(1)) cycle
          j = j + 1
          ep%p(ivnl)%jxyz(j) = k
          do d = 1, sb%dim
             ep%dp(d, ivnl)%jxyz(j) = k
          enddo
          do d = 1, sb%dim
             ep%lso(d, ivnl)%jxyz(j) = k
          enddo
          exit
        end do
      end do

      call pop_sub()
    end subroutine build_kb_sphere

    subroutine allocate_nl_part(ivnl)
      integer, intent(in) :: ivnl
      integer :: j, c, d
      call push_sub('epot.allocate_nl_part')

      j = ep%p(ivnl)%n
      c = ep%p(ivnl)%c

      allocate(ep%p(ivnl)%a(j, c), &
               ep%p(ivnl)%uvu(c, c))
      ep%p(ivnl)%a(:,:)     = M_ZERO
      ep%p(ivnl)%uvu(:,:) = M_ZERO
      ep%p(ivnl)%b => ep%p(ivnl)%a

      do d = 1, sb%dim
         allocate(ep%dp(d, ivnl)%a(j, c), ep%dp(d, ivnl)%b(j, c), &
                  ep%dp(d, ivnl)%uvu(c, c))
         ep%dp(d, ivnl)%a(:,:)     = M_ZERO
         ep%dp(d, ivnl)%b(:,:)     = M_ZERO
         ep%dp(d, ivnl)%uvu(:,:) = M_ZERO
      enddo

      do d = 1, sb%dim
         allocate(ep%lso(d, ivnl)%a(j, c), ep%lso(d, ivnl)%b(j, c), &
                  ep%lso(d, ivnl)%uvu(c, c))
         ep%lso(d, ivnl)%a(:,:)     = M_ZERO
         ep%lso(d, ivnl)%b(:,:)     = M_ZERO
         ep%lso(d, ivnl)%uvu(:,:) = M_ZERO
      enddo

      if(sb%periodic_dim/=0) then
        allocate(ep%p(ivnl)%phases(ep%p(ivnl)%n, st%d%nik))
        do d = 1, sb%dim
           allocate(ep%dp(d,ivnl)%phases(ep%dp(d, ivnl)%n, st%d%nik))
        enddo
        do d = 1, sb%dim
           allocate(ep%lso(d,ivnl)%phases(ep%lso(d, ivnl)%n, st%d%nik))
        enddo
      endif

      call pop_sub()
    end subroutine allocate_nl_part

    subroutine build_nl_part(ivnl, l, lm)
      integer, intent(in) :: ivnl, l, lm

      integer :: i, j, k, d
      FLOAT :: r, x(3), x_in(3)
      FLOAT :: so_uv, so_duv(3)
      FLOAT :: v, dv(3)

      integer :: c, n
      CMPLX, allocatable :: grad_so(:, :, :)

      call push_sub('epot.build_nl_part')

      c = ep%lso(1, ivnl)%c
      n = ep%lso(1, ivnl)%n
      allocate(grad_so(n, 3, c))

      j_loop: do j = 1, n
        x_in(:) = m%x(ep%p(ivnl)%jxyz(j), :)
         k_loop: do k = 1, 3**sb%periodic_dim
          x(:) = x_in(:) - sb%shift(k,:)
          r=sqrt(sum((x-a%x)*(x-a%x)))
          if (r > s%ps%rc_max + m%h(1)) cycle
          x = x - a%x
          i_loop : do i = 1, c
                if(l .ne. s%ps%L_loc) then
                  call specie_get_nl_part(s, x, l, lm, i, v, dv(1:3))
                  ep%p(ivnl)%a(j, i) = v
                  do d = 1, sb%dim
                     ep%dp(d, ivnl)%a(j, i) = dv(d)
                     ep%dp(d, ivnl)%b(j, i) = v
                  enddo
                end if
                if(l>0 .and. s%ps%so_l_max>=0) then
                  call specie_get_nl_part(s, x, l, lm, i, so_uv, so_duv(:), so=.true.)
                  grad_so(j, 1:3, i) = so_duv(1:3)
                  do d = 1, sb%dim
                     ep%lso(d, ivnl)%a(j, i) = so_uv
                  enddo
                endif
          end do i_loop
         end do k_loop
      end do j_loop

      if(reltype == 1) then ! SPIN_ORBIT
         c = ep%lso(1, ivnl)%c
         do j = 1, ep%lso(1, ivnl)%n
            x_in(1:3) = m%x(ep%lso(1, ivnl)%jxyz(j), 1:3)
            x(1:3) = x_in(1:3) !- sb%shift(k, 1:3) !???
            x = x - a%x
            r = sqrt(sum(x*x))
            ep%lso(1, ivnl)%b(j, 1:c) = (x(2)*grad_so(j, 3, 1:c) - x(3)*grad_so(j, 2, 1:c))
            ep%lso(2, ivnl)%b(j, 1:c) = (x(3)*grad_so(j, 1, 1:c) - x(1)*grad_so(j, 3, 1:c))
            ep%lso(3, ivnl)%b(j, 1:c) = (x(1)*grad_so(j, 2, 1:c) - x(2)*grad_so(j, 1, 1:c))
         enddo
      endif

      ! and here we calculate the uVu
      if(s%ps%flavour == PS_TM2) then
            if(l .ne. s%ps%L_loc) then
              ep%p(ivnl)%uvu(1, 1) = s%ps%h(l, 1, 1)
              do d = 1, sb%dim
                 ep%dp(d, ivnl)%uvu(1, 1) = s%ps%h(l, 1, 1)
              enddo
            end if
      else
        ep%p(ivnl)%uvu(1:3, 1:3) = s%ps%h(l, 1:3, 1:3)
        do d = 1, sb%dim
           ep%dp(d, ivnl)%uvu(1:3, 1:3) = s%ps%h(l, 1:3, 1:3)
        enddo
        do d = 1, sb%dim
           ep%lso(d, ivnl)%uvu(1:3, 1:3) = s%ps%k(l, 1:3, 1:3)
        enddo
      end if

      if(sb%periodic_dim/=0) then
        do j = 1, ep%p(ivnl)%n
           x(:) = m%x(ep%p(ivnl)%jxyz(j), :)
           do k=1, st%d%nik
              ep%p(ivnl)%phases(j, k) = exp(M_zI*sum(st%d%kpoints(:, k)*x(:)))
              do d = 1, sb%dim
                 ep%dp(d, ivnl)%phases(j, k) = exp(M_zI*sum(st%d%kpoints(:, k)*x(:)))
              enddo
              do d = 1, sb%dim
                 ep%lso(d, ivnl)%phases(j, k) = exp(M_zI*sum(st%d%kpoints(:, k)*x(:)))
              enddo
           end do
        end do
      endif

      call pop_sub(); return
    end subroutine build_nl_part

  end subroutine epot_generate

  subroutine epot_generate_classic(ep, m, geo)
    type(epot_type),     intent(inout) :: ep
    type(mesh_type),     intent(in)    :: m
    type(geometry_type), intent(in)    :: geo

    integer i, ia
    FLOAT :: r, rc

    call push_sub('epot.epot_generate_classic')

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


  function epot_laser_scalar_pot(np, gr, ep, t) result(v)
    integer, intent(in) :: np
    type(grid_type), intent(in) :: gr
    type(epot_type), intent(in)  :: ep
    FLOAT, intent(in)  :: t
    FLOAT :: v(np)

    call laser_potential(gr%sb, ep%no_lasers, ep%lasers, t, gr%m, v)

  end function epot_laser_scalar_pot


  subroutine epot_laser_vector_pot(sb, ep, t, a)
    type(simul_box_type), intent(in) :: sb
    type(epot_type), intent(in)  :: ep
    FLOAT,           intent(in)  :: t
    FLOAT,           intent(out) :: a(sb%dim)

    call laser_vector_field(sb, ep%no_lasers, ep%lasers, t, a)

  end subroutine epot_laser_vector_pot


  subroutine epot_laser_field(sb, ep, t, e)
    type(simul_box_type), intent(in) :: sb
    type(epot_type), intent(in)  :: ep
    FLOAT,           intent(in)  :: t
    FLOAT,           intent(out) :: e(sb%dim)

    call laser_field(sb, ep%no_lasers, ep%lasers, t, e)

  end subroutine epot_laser_field

#include "undef.F90"
#include "real.F90"
#include "epot_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "epot_inc.F90"

end module external_pot


