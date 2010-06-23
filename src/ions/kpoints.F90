!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: kpoints.F90 5635 2009-06-25 14:55:35Z nitsche $

#include "global.h"
  
module kpoints_m
  use datasets_m
  use geometry_m
  use global_m
  use loct_m
  use parser_m
  use messages_m
  use profiling_m
  use symmetries_m
  use unit_m
  use unit_system_m
  
  implicit none
  
  private
  
  public ::                  &
    kpoints_grid_t,          &
    kpoints_t,               &
    kpoints_init,            &
    kpoints_end,             &
    kpoints_copy,            &
    kpoints_number,          &
    kpoints_get_weight,      &
    kpoints_get_point,       &
    kpoints_write_info,      &
    kpoints_point_is_gamma

  type kpoints_grid_t
    FLOAT, pointer :: point(:, :)
    FLOAT, pointer :: red_point(:, :)
    FLOAT, pointer :: weight(:)
    integer        :: npoints
  end type kpoints_grid_t

  type kpoints_t
    type(kpoints_grid_t) :: full
    type(kpoints_grid_t) :: reduced

    integer        :: method

    logical        :: use_symmetries
    logical        :: use_time_reversal

    ! For the modified Monkhorst-Pack  scheme
    integer        :: nik_axis(MAX_DIM)    ! number of MP divisions
    FLOAT          :: shifts(MAX_DIM)      ! 
  end type kpoints_t

  integer, parameter ::                &
    KPOINTS_GAMMA       =  1,          &
    KPOINTS_MONKH_PACK  =  2,          &
    KPOINTS_USER        =  3

contains

  subroutine kpoints_grid_init(this, npoints)
    type(kpoints_grid_t), intent(out) :: this
    integer,             intent(in)  :: npoints

    this%npoints = npoints
    SAFE_ALLOCATE(this%red_point(1:3, npoints))
    SAFE_ALLOCATE(this%point(1:3, npoints))
    SAFE_ALLOCATE(this%weight(npoints))
  end subroutine kpoints_grid_init

  ! ---------------------------------------------------------

  subroutine kpoints_grid_end(this)
    type(kpoints_grid_t), intent(out) :: this

    SAFE_DEALLOCATE_P(this%red_point)
    SAFE_DEALLOCATE_P(this%point)
    SAFE_DEALLOCATE_P(this%weight)
  end subroutine kpoints_grid_end

  ! ---------------------------------------------------------

  subroutine kpoints_grid_copy(bb, aa)
    type(kpoints_grid_t), intent(in)  :: bb
    type(kpoints_grid_t), intent(out) :: aa

    call push_sub('kpoints.kpoints_grid_copy')
    
    call kpoints_grid_init(aa, bb%npoints)
    
    aa%weight = bb%weight
    aa%point  = bb%point
    aa%red_point  = bb%red_point

    call pop_sub('kpoints.kpoints_grid_copy')
  end subroutine kpoints_grid_copy

  ! ---------------------------------------------------------

  subroutine kpoints_init(this, symm, dim, periodic_dim, rlattice, klattice, geo)
    type(kpoints_t),    intent(out) :: this
    type(symmetries_t), intent(in)  :: symm
    integer,            intent(in)  :: dim
    integer,            intent(in)  :: periodic_dim
    FLOAT,              intent(in)  :: rlattice(:,:), klattice(:,:)
    type(geometry_t),   intent(in)  :: geo

    integer :: ik

    call push_sub('kpoints.kpoints_init')

    !%Variable KPointsUseSymmetries
    !%Type logical
    !%Default no
    !%Section Mesh::KPoints
    !%Description
    !% This variable defines whether symmetries are taken into account
    !% or not for the choice of <i>k</i>-points. If it is set to no, the <i>k</i>-point
    !% sampling will range over the full Brillouin zone.
    !% Symmetries should not be used whenever a perturbation is applied to the system.
    !%End
    call parse_logical(datasets_check('KPointsUseSymmetries'), .false., this%use_symmetries)

    !%Variable KPointsUseTimeReversal
    !%Type logical
    !%Default yes
    !%Section Mesh::KPoints
    !%Description
    !% <b>WARNING: This variable does not seem to work.</b>
    !% This variable defines whether time-reversal symmetry is taken into account
    !% or not for the choice of <i>k</i>-points. If it is set to no, the <i>k</i>-point
    !% sampling will not be reduced according to time-reversal symmetry. The default is
    !% yes. If <tt>KPointsUseSymmetries = no</tt>, this variable is ignored, and time-reversal
    !% symmetry is not used.
    !%End
    call parse_logical(datasets_check('KPointsUseTimeReversal'), .true., this%use_time_reversal)

    if(periodic_dim == 0) then
      this%method = KPOINTS_GAMMA
      call read_MP(.true.)
    else
      if(read_user_kpoints()) then
        this%method = KPOINTS_USER
      else
        this%method = KPOINTS_MONKH_PACK
        call read_MP(.false.)
      end if
    end if

    write(message(1),'(a)') ' '
    write(message(2),'(1x,i3,a)') this%reduced%npoints, ' k-points generated from parameters :'
    write(message(3),'(1x,a)') '---------------------------------------------------'
    write(message(4),'(4x,a,3i5,6x,a,3f6.2)') 'n =', this%nik_axis(1:3), 's = ', this%shifts(1:3)
    write(message(5),'(a)') ' '
    write(message(6),'(a)') ' index |    weight    |              coordinates              |'
    call write_info(6)

    do ik = 1, this%reduced%npoints
      write(message(1), &
        '(i6,a,f12.6,a,3f12.6, a)') ik, " | ", this%reduced%weight(ik), " | ", this%reduced%red_point(1:3, ik), "  |"
      call write_info(1)
    end do
  
    write(message(1),'(a)') ''
    call write_info(1)

    call pop_sub('kpoints.kpoints_init')

  contains
    ! ---------------------------------------------------------
    subroutine read_mp(gamma_only)
      logical, intent(in) :: gamma_only

      logical       :: gamma_only_
      integer       :: ii
      type(block_t) :: blk

      call push_sub('kpoints.kpoints_init.read_MP')

      call messages_obsolete_variable('KPointsMonkhorstPack', 'KPointsGrid')

      !%Variable KPointsGrid
      !%Type block
      !%Default Gamma-point only
      !%Section Mesh::KPoints
      !%Description
      !% When this block is given (and the <tt>KPoints</tt> block is not present),
      !% <i>k</i>-points are distributed in a uniform grid.
      !%
      !% The first row of the block is a triplet of integers defining
      !% the number of <i>k</i>-points to be used along each direction
      !% in reciprocal space. The numbers refer to the whole Brillouin
      !% zone, and the actual number of <i>k</i>-points is usually
      !% reduced exploiting the symmetries of the system.  By default
      !% the grid will always include the Gamma point. An optional
      !% second row can specify a shift in the <i>k</i>-points.
      !%
      !% For example, the following input samples the BZ with 100 points in the 
      !% <i>xy</i>-plane of reciprocal space:
      !%
      !% <tt>%KPointsGrid
      !% <br>&nbsp;&nbsp;10 | 10 | 1
      !% <br>%</tt>
      !%
      !%End

      gamma_only_ = gamma_only
      if(.not.gamma_only_) gamma_only_ = (parse_block(datasets_check('KPointsGrid'), blk) .ne. 0)

      this%nik_axis(:) = 1
      this%shifts(:)   = M_ZERO

      if(.not.gamma_only_) then
        do ii = 1, periodic_dim
          call parse_block_integer(blk, 0, ii - 1, this%nik_axis(ii))
        end do

        if (any(this%nik_axis < 1)) then
          message(1) = 'Input: KPointsGrid is not valid'
          call write_fatal(1)
        end if

        if(parse_block_n(blk) > 1) then ! we have a shift
          do ii = 1, periodic_dim
            call parse_block_float(blk, 1, ii - 1, this%shifts(ii))
          end do
        end if

        call parse_block_end(blk)
      end if

      call kpoints_grid_init(this%full, product(this%nik_axis))

      call kpoints_grid_generate(periodic_dim, this%nik_axis, this%shifts, this%full%npoints, this%full%red_point)

      this%full%weight = M_ONE/real(this%full%npoints, REAL_PRECISION)

      call kpoints_grid_copy(this%full, this%reduced)

      if(this%use_symmetries) then
        call kpoints_grid_reduce(symm, this%use_time_reversal, &
          this%reduced%npoints, this%reduced%red_point, this%reduced%weight)
      end if

      do ik = 1, this%full%npoints
        call kpoints_to_absolute(klattice, this%full%red_point(:, ik), this%full%point(:, ik))
      end do

      do ik = 1, this%reduced%npoints
        call kpoints_to_absolute(klattice, this%reduced%red_point(:, ik), this%reduced%point(:, ik))
      end do

      call pop_sub('kpoints.kpoints_init.read_MP')
    end subroutine read_mp

    ! ---------------------------------------------------------
    logical function read_user_kpoints()
      type(block_t) :: blk
      logical :: reduced
      integer :: ik, idir

      call push_sub('kpoints.kpoints_init.read_user_kpoints')

      !%Variable KPoints
      !%Type block
      !%Section Mesh::KPoints
      !%Description
      !% This block defines an explicit set of <i>k</i>-points and their weights for
      !% a periodic-system calculation. The first column is the weight
      !% of each <i>k</i>-point and the following are the components of the <i>k</i>-point
      !% vector. You only need to specify the components for the
      !% periodic directions. Note that the <i>k</i>-points should be given in
      !% Cartesian coordinates (not in reduced coordinates), <i>i.e.</i>
      !% what Octopus writes in a line in the ground-state standard output as
      !% <tt>#k =   1, k = (    0.154000,    0.154000,    0.154000)</tt>.
      !%
      !% For example, if you want to include only the Gamma point, you can
      !% use:
      !%
      !% <tt>%KPoints
      !% <br>&nbsp;&nbsp;1.0 | 0 | 0 | 0
      !% <br>%</tt>
      !%
      !%End

      !%Variable KPointsReduced
      !%Type block
      !%Section Mesh::KPoints
      !%Description
      !% Same as the block <tt>KPoints</tt> but this time the input is given in reduced 
      !% coordinates.
      !%End

      reduced = .false.
      if(parse_block(datasets_check('KPoints'), blk) /= 0) then
        if(parse_block(datasets_check('KPointsReduced'), blk) == 0) then
          reduced = .true.
        else
          read_user_kpoints = .false.
          call pop_sub('kpoints.kpoints_init.read_user_kpoints')
          return
        end if
      end if
      read_user_kpoints = .true.

      call kpoints_grid_init(this%full, parse_block_n(blk))

      this%full%red_point = M_ZERO
      this%full%weight = M_ZERO

      if(reduced) then
        do ik = 1, this%full%npoints
          call parse_block_float(blk, ik - 1, 0, this%full%weight(ik))
          do idir = 1, periodic_dim
            call parse_block_float(blk, ik - 1, idir, this%full%red_point(idir, ik))
          end do
          ! generate also the absolute coordinates
          call kpoints_to_absolute(klattice, this%full%red_point(:, ik), this%full%point(:, ik))
        end do
      else
        do ik = 1, this%full%npoints
          call parse_block_float(blk, ik - 1, 0, this%full%weight(ik))
          do idir = 1, periodic_dim
            call parse_block_float(blk, ik - 1, idir, this%full%point(idir, ik), unit_one/units_inp%length)
          end do
          ! generate also the reduced coordinates
          call kpoints_to_reduced(rlattice, this%full%point(:, ik), this%full%red_point(:, ik))
        end do
      end if
      call parse_block_end(blk)


      ! for the moment we do not apply symmetries to user kpoints
      call kpoints_grid_copy(this%full, this%reduced)

      write(message(1), '(a,i4,a)') 'Input: ', this%full%npoints, ' k-points were read from the input file'
      call write_info(1)

      call pop_sub('kpoints.kpoints_init.read_user_kpoints')
    end function read_user_kpoints

  end subroutine kpoints_init


  ! ---------------------------------------------------------
  subroutine kpoints_end(this)
    type(kpoints_t), intent(inout) :: this

    call kpoints_grid_end(this%full)
    call kpoints_grid_end(this%reduced)

  end subroutine kpoints_end


  ! ---------------------------------------------------------
  subroutine kpoints_to_absolute(klattice, kin, kout)
    FLOAT, intent(in)  :: klattice(:,:), kin(:)
    FLOAT, intent(out) :: kout(:)

    integer :: ii
    
    ! short

    kout(:) = M_ZERO
    do ii = 1, MAX_DIM
      kout(:) = kout(:) + kin(ii)*klattice(:, ii)
    end do

  end subroutine kpoints_to_absolute


  ! ---------------------------------------------------------
  subroutine kpoints_to_reduced(rlattice, kin, kout)
    FLOAT, intent(in)  :: rlattice(:,:)
    FLOAT, intent(in)  :: kin(:)
    FLOAT, intent(out) :: kout(:)

    integer :: ii

    ! short

    kout(:) = M_ZERO
    do ii = 1, MAX_DIM
      kout(:) = kout(:) + kin(ii)*rlattice(ii, :)
    end do
    kout(:) = kout(:) / (M_TWO*M_PI)

  end subroutine kpoints_to_reduced


  ! ---------------------------------------------------------
  subroutine kpoints_copy(kin, kout)
    type(kpoints_t), intent(in)  :: kin
    type(kpoints_t), intent(out) :: kout

    call push_sub('kpoints.kpoints_copy')

    kout%method = kin%method

    call kpoints_grid_copy(kin%full, kout%full)
    call kpoints_grid_copy(kin%reduced, kout%reduced)

    kout%use_symmetries = kin%use_symmetries
    kout%use_time_reversal = kin%use_time_reversal

    kout%nik_axis = kin%nik_axis
    kout%shifts = kin%shifts

    call pop_sub('kpoints.kpoints_copy')
  end subroutine kpoints_copy

  ! ----------------------------------------------------------

  integer pure function kpoints_number(this) result(number)
    type(kpoints_t), intent(in) :: this

    number = this%reduced%npoints
    
  end function kpoints_number

  ! ----------------------------------------------------------

  pure function kpoints_get_point(this, ik) result(point)
    type(kpoints_t), intent(in) :: this
    integer,         intent(in) :: ik
    FLOAT                       :: point(1:3)

    point(1:3) = this%reduced%point(1:3, ik)

  end function kpoints_get_point

  ! ----------------------------------------------------------

  FLOAT pure function kpoints_get_weight(this, ik) result(weight)
    type(kpoints_t), intent(in) :: this
    integer,         intent(in) :: ik

    weight = this%reduced%weight(ik)

  end function kpoints_get_weight

  ! ----------------------------------------------------------
  ! Generates the k-points grid
  ! Sets up a uniform array of k-points. Use a modification of the normal Monkhorst-Pack scheme, 
  ! which is equivalent to the normal MP scheme in the case of even number of kpoints (i.e. naxis (i) even)  
  ! used with a shift of (1/2, 1/2, 1/2)
  ! For the original MP scheme, see (PRB 13, 518 (1976))
  ! and (PRB 16, 1748 (1977))
  ! naxis(i) are the number of points in the three
  ! directions dermined by the lattice vectors.
  ! shift(i) and sz shift the grid of integration points from the origin.

  subroutine kpoints_grid_generate(periodic_dim, naxis, shift, nkpoints, kpoints)  
    integer,           intent(in)  :: periodic_dim
    integer,           intent(in)  :: naxis(1:MAX_DIM)
    FLOAT,             intent(in)  :: shift(1:MAX_DIM)
    integer,           intent(out) :: nkpoints
    FLOAT,             intent(out) :: kpoints(:, :)
  
    FLOAT :: dx(1:MAX_DIM)
    integer :: ii, jj, kk, idir, ix(1:MAX_DIM)

    call push_sub('kpoints.kpoints_grid_generate')
    
    dx(1:MAX_DIM) = M_ONE/real(2*naxis(1:MAX_DIM), REAL_PRECISION)

    kpoints = M_ZERO
    nkpoints = 0
    ix=1
    do ii = 1, naxis(1)
       do jj = 1, naxis(2)
          do kk = 1, naxis(3)
             ix(1:3) = (/ii, jj, kk/)
             nkpoints = nkpoints + 1
             do idir = 1, periodic_dim
                if ( (mod(naxis(idir),2) .ne. 0 ) ) then    
                   kpoints(idir, nkpoints) = &
                        & (real(2*ix(idir) - naxis(idir) - 1, REAL_PRECISION) + 2*shift(idir))*dx(idir)
                else 
                   kpoints(idir, nkpoints) = &
                        & (real(2*ix(idir) - naxis(idir) , REAL_PRECISION) + 2*shift(idir))*dx(idir)
                end if
		! bring back point to first Brillouin zone
                kpoints(idir, nkpoints) = mod(kpoints(idir, nkpoints) + M_HALF, M_ONE) - M_HALF
             end do
          end do
       end do
    end do
    
    call pop_sub('kpoints.kpoints_grid_generate')

  end subroutine kpoints_grid_generate
  
  ! --------------------------------------------------------------------------------------------

  subroutine kpoints_grid_reduce(symm, time_reversal, nkpoints, kpoints, weights)
    type(symmetries_t), intent(in)    :: symm
    logical,            intent(in)    :: time_reversal
    integer,            intent(inout) :: nkpoints
    FLOAT,              intent(inout) :: kpoints(:, :)
    FLOAT,              intent(out)   :: weights(:)

    integer :: nreduced
    FLOAT, allocatable :: reduced(:, :)
    
    FLOAT :: dw
    integer ik, iop, ik2
    FLOAT :: tran(3), tran_inv(3)
    integer, allocatable :: kmap(:)

    call push_sub('kpoints.kpoints_grid_reduce')

    ! reduce to irreducible zone

    ! kmap is used to mark reducible k-points and also to
    ! map reducible to irreducible k-points

    SAFE_ALLOCATE(kmap(1:nkpoints))
    SAFE_ALLOCATE(reduced(1:3, 1:nkpoints))

    forall(ik = 1:nkpoints) kmap(ik) = ik
    
    dw = M_ONE/real(nkpoints, REAL_PRECISION)

    nreduced = 0

    do ik = 1, nkpoints
      if (kmap(ik) /= ik) cycle
      
      ! new irreducible point
      ! mareduced with negative kmap
      
      nreduced = nreduced + 1
      reduced(1:3, nreduced) = kpoints(1:3, ik)
      
      kmap(ik) = -nreduced
      
      weights(nreduced) = dw

      if (ik == nkpoints) cycle
      
      ! operate with the symmetry operations
      
      do iop = 1, symmetries_number(symm)
        call symmetries_apply_kpoint(symm, iop, reduced(:, nreduced), tran)
        tran_inv(1:3) = -tran(1:3)
           
        ! remove (mark) k-points related to irreducible reduced by symmetry
        do ik2 = ik + 1, nkpoints
          if (kmap(ik2) /= ik2) cycle
          
          ! both the transformed rk ...
          if (all( abs(tran(1:3) - kpoints(1:3, ik2)) <= CNST(1.0e-5))) then
            weights(nreduced) = weights(nreduced) + dw
            kmap(ik2) = nreduced
            cycle
          end if

          ! and its inverse
          if (time_reversal .and. all(abs(tran_inv(1:3) - kpoints(1:3, ik2)) <= CNST(1.0e-5)) ) then
            weights(nreduced) = weights(nreduced) + dw
            kmap(ik2) = nreduced
          end if
          
        end do
      end do
    end do
    
    nkpoints = nreduced
    do ik = 1, nreduced
      kpoints(1:3, ik) = reduced(1:3, ik)
    end do

    call pop_sub('kpoints.kpoints_grid_reduce')
  end subroutine kpoints_grid_reduce

  ! ---------------------------------------------------------

  subroutine kpoints_write_info(this, iunit)
    type(kpoints_t),    intent(in) :: this
    integer,            intent(in) :: iunit
    
    integer :: ik, idir
    
    call push_sub('kpoints.kpoints_write_info')
    
    if(this%method == KPOINTS_MONKH_PACK) then
      write(message(1),'(a,9(i3,1x))') 'Number of k-points in each direction = ', this%nik_axis(1:3)
      call write_info(1, iunit)
    else
      ! a Monkhorst-Pack grid was not used
      write(message(1),'(a,9(i3,1x))') 'Number of k-points = ', kpoints_number(this)
      call write_info(1, iunit)
    endif
    
    write(message(1), '(3x,a,7x,a,9x,a,9x,a,8x,a)') 'ik', 'K_x', 'K_y', 'K_z', 'Weight'
    message(2) = '   --------------------------------------------------'
    call write_info(2, iunit, verbose_limit=80)
    
    do ik = 1, kpoints_number(this)
      write(message(1),'(i4,1x,4f12.4)') ik, &
        (units_from_atomic(unit_one/units_out%length, this%reduced%red_point(idir, ik)), idir=1, 3), kpoints_get_weight(this, ik)
      call write_info(1, iunit, verbose_limit=80)
    end do
    
    call pop_sub('kpoints.kpoints_write_info')
  end subroutine kpoints_write_info
  

  ! ---------------------------------------------------------
  logical pure function kpoints_point_is_gamma(this, ik) result(is_gamma)
    type(kpoints_t),    intent(in) :: this
    integer,            intent(in) :: ik
    
    is_gamma = (maxval(abs(kpoints_get_point(this, ik))) < M_EPSILON)
    
  end function kpoints_point_is_gamma

end module kpoints_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
