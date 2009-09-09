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
  use loct_parser_m
  use messages_m
  use profiling_m
  use units_m

  implicit none

  private
  public ::       &
    kpoints_t,    &
    kpoints_init, &
    kpoints_end,  &
    kpoints_copy


  type kpoints_t
    integer        :: method
    
    integer        :: nik_full             ! number of k-points in full Brillouin zone
    FLOAT, pointer :: points_full(:,:)     ! kpoints in absolute coordinates
    FLOAT, pointer :: points_full_red(:,:) ! kpoints in reduced coordinates
    FLOAT, pointer :: weights_full(:)      ! weights for the k-point integrations

    integer        :: nik                  ! number of k-points in the irreducible zone
    FLOAT, pointer :: points(:,:)          ! kpoints in absolute coordinates
    FLOAT, pointer :: points_red(:,:)      ! kpoints in relative coordinates
    FLOAT, pointer :: weights(:)           ! weights for the k-point integrations

    logical        :: use_symmetries
    logical        :: use_inversion

    ! For the Monkhorst-Pack scheme
    integer        :: nik_axis(MAX_DIM)    ! number of MP divisions
    FLOAT          :: shifts(MAX_DIM)      ! 

  end type kpoints_t

  integer, parameter ::        &
    KPOINTS_GAMMA       =  1,  &
    KPOINTS_MONKH_PACK  =  2,  &
    KPOINTS_USER        =  3

contains


! ---------------------------------------------------------
subroutine kpoints_init(this, dim, periodic_dim, rlattice, klattice, geo)
  type(kpoints_t),    intent(out) :: this
  integer,            intent(in)  :: dim, periodic_dim
  FLOAT,              intent(in)  :: rlattice(:,:), klattice(:,:)
  type(geometry_t),   intent(in)  :: geo

  call push_sub('kpoints.kpoints_init')

  if(periodic_dim==0) then
    this%method = KPOINTS_GAMMA
    call read_MP(.true.)
    call generate_MP()
  else
    if(read_user_kpoints()) then
      this%method = KPOINTS_USER
    else
      this%method = KPOINTS_MONKH_PACK
      call read_MP(.false.)
      call generate_MP()
    end if
  end if

  call pop_sub()

contains
  ! ---------------------------------------------------------
  subroutine read_MP(gamma_only)
    logical, intent(in) :: gamma_only

    logical       :: gamma_only_
    integer       :: ii
    type(block_t) :: blk

    call push_sub('kpoints.kpoints_init.read_MP')

    !%Variable KPointsMonkhorstPack
    !%Type block
    !%Default Gamma-point only
    !%Section Mesh::KPoints
    !%Description
    !% When this block is given (and the <tt>KPoints</tt> block is not present),
    !% <i>k</i>-points are generated in a Monkhorst-Pack grid.
    !%
    !% The first row of the block is a triplet of integers defining the
    !% number of <i>k</i>-points to be used along each direction in 
    !% reciprocal space. The numbers refer to the whole Brillouin zone,
    !% and the actual number of <i>k</i>-points is usually reduced exploiting
    !% the symmetries of the system.  An optional second row can specify
    !% a shift in the <i>k</i>-points.
    !%
    !% For example, the following input samples the BZ with 100 points in the 
    !% <i>xy</i>-plane of reciprocal space:
    !%
    !% <tt>%KPointsMonkhorstPack
    !% <br>&nbsp;&nbsp;10 | 10 | 1
    !% <br>%</tt>
    !%
    !%End
    
    gamma_only_ = gamma_only
    if(.not.gamma_only_) gamma_only_ = (loct_parse_block(datasets_check('KPointsMonkhorstPack'), blk) .ne. 0)

    this%nik_axis(:) = 1
    this%shifts(:)   = M_ZERO

    if(.not.gamma_only_) then
      do ii = 1, periodic_dim
        call loct_parse_block_int(blk, 0, ii-1, this%nik_axis(ii))
      end do

      if (any(this%nik_axis < 1)) then
        message(1) = 'Input: KPointsMonkhorstPack is not valid'
        call write_fatal(1)
      end if

      if(loct_parse_block_n(blk) > 1) then ! we have a shift
        do ii = 1, periodic_dim
          call loct_parse_block_float(blk, 1, ii-1, this%shifts(ii))
        end do
      end if

      call loct_parse_block_end(blk)
    end if
      
    this%nik_full = product(this%nik_axis(:))

    SAFE_ALLOCATE(this%points_full    (1:MAX_DIM, 1:this%nik_full))
    SAFE_ALLOCATE(this%points_full_red(1:MAX_DIM, 1:this%nik_full))
    SAFE_ALLOCATE(this%weights_full(1:this%nik_full))
    
    this%points_full(:,:)     = M_ZERO
    this%points_full_red(:,:) = M_ZERO
    this%weights_full(:)      = M_ONE / TOFLOAT(this%nik_full)

    call pop_sub()
  end subroutine read_MP

  
  ! ---------------------------------------------------------
  subroutine generate_MP()
    !generate k-points using the MP scheme
    FLOAT :: dx(1:MAX_DIM)
    integer :: ii, jj, kk, ix(1:MAX_DIM), idir, ikp, odd_shifts(1:MAX_DIM)
    
    call push_sub('kpoints.kpoints_init.generate_MP')

    dx(:) = M_ONE/TOFLOAT(2*this%nik_axis(:))
    odd_shifts(:) = M_ZERO
    do idir = 1, periodic_dim
      if(mod(this%nik_axis(idir), 2) == 1) odd_shifts(idir) = 1
    end do

    ikp = 0
    ix(:) = M_ZERO
    do ii = 1, this%nik_axis(1)
      do jj = 1, this%nik_axis(2)
        do kk = 1, this%nik_axis(3)
          ikp = ikp + 1
          ix(1:3) = (/ii, jj, kk/)

          forall(idir = 1:periodic_dim)
            this%points_full_red(idir, ikp) = TOFLOAT(2*ix(idir) - this%nik_axis(idir) - odd_shifts(idir))
            this%points_full_red(idir, ikp) = (this%points_full(idir, ikp) + this%shifts(idir))*dx(idir)
          end forall
        end do
      end do
    end do
    
    call pop_sub()

  end subroutine generate_MP


  ! ---------------------------------------------------------
  logical function read_user_kpoints()
    type(block_t) :: blk
    logical :: reduced
    integer :: ik, idim

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

    this%nik_full = 0
    reduced = .false.
    if(loct_parse_block(datasets_check('KPoints'), blk).ne.0) then
      if(loct_parse_block(datasets_check('KPointsReduced'), blk) == 0) then
        reduced = .true.
      else
        read_user_kpoints = .false.
        call pop_sub(); return
      end if
    end if
    read_user_kpoints = .true.
    
    this%nik_full = loct_parse_block_n(blk)

    SAFE_ALLOCATE(this%points_full(1:MAX_DIM, 1:this%nik_full))
    SAFE_ALLOCATE(this%points_full_red(1:MAX_DIM, 1:this%nik_full))
    SAFE_ALLOCATE(this%weights_full(1:this%nik_full))

    this%points_full  = M_ZERO
    this%weights_full = M_ZERO

    if(reduced) then
      do ik = 1, this%nik_full
        call loct_parse_block_float(blk, ik - 1, 0, this%weights_full(ik))
        do idim = 1, periodic_dim
          call loct_parse_block_float(blk, ik - 1, idim, this%points_full_red(idim, ik))
        end do
      end do

      ! generate also the absolute coordinates
      do ik = 1, this%nik_full
        call kpoints_to_absolute(klattice, this%points_full_red(:,ik), this%points_full(:,ik))
      end do
    else
      do ik = 1, this%nik_full
        call loct_parse_block_float(blk, ik - 1, 0, this%weights_full(ik))
        do idim = 1, periodic_dim
          call loct_parse_block_float(blk, ik - 1, idim, this%points_full(idim, ik))
        end do
      end do

      do ik = 1, this%nik_full
        ! k-points have 1/length units
        this%points_full(:, ik) = units_to_atomic(unit_one/units_inp%length, this%points_full(:, ik))

        ! generate also the reduced coordinates
        call kpoints_to_reduced(rlattice, this%points_full(:,ik), this%points_full_red(:,ik))
      end do
    end if
    call loct_parse_block_end(blk)

    write(message(1), '(a,i4,a)') 'Input: ', this%nik_full, ' k-points were read from the input file'
    call write_info(1)

    call pop_sub()
  end function read_user_kpoints

end subroutine kpoints_init


! ---------------------------------------------------------
subroutine kpoints_end(this)
  type(kpoints_t), intent(inout) :: this

  SAFE_DEALLOCATE_P(this%points_full)
  SAFE_DEALLOCATE_P(this%points_full_red)
  SAFE_DEALLOCATE_P(this%weights_full)
  SAFE_DEALLOCATE_P(this%points)
  SAFE_DEALLOCATE_P(this%points_red)
  SAFE_DEALLOCATE_P(this%weights)
  
end subroutine kpoints_end


! ---------------------------------------------------------
subroutine kpoints_to_absolute(klattice, kin, kout)
  FLOAT, intent(in)  :: klattice(:,:), kin(:)
  FLOAT, intent(out) :: kout(:)

  integer :: ii

  call push_sub('kpoints.kpoints_to_absolute')

  kout(:) = M_ZERO
  do ii = 1, MAX_DIM
    kout(:) = kout(:) + kin(ii)*klattice(:, ii)
  end do

  call pop_sub()

end subroutine kpoints_to_absolute


! ---------------------------------------------------------
subroutine kpoints_to_reduced(rlattice, kin, kout)
  FLOAT, intent(in)  :: rlattice(:,:)
  FLOAT, intent(in)  :: kin(:)
  FLOAT, intent(out) :: kout(:)

  integer :: ii

  call push_sub('kpoints.kpoints_to_reduced')

  kout(:) = M_ZERO
  do ii = 1, MAX_DIM
    kout(:) = kout(:) + kin(ii)*rlattice(ii, :)
  end do
  kout(:) = kout(:) / (M_TWO*M_PI)

  call pop_sub()
end subroutine kpoints_to_reduced


! ---------------------------------------------------------
subroutine kpoints_copy(kout, kin)
  type(kpoints_t), intent(out) :: kout
  type(kpoints_t), intent(in)  :: kin

  call push_sub('kpoints.kpoints_copy')

  kout%method = kin%method
  kout%nik_full = kin%nik_full
  call loct_pointer_copy(kout%points_full, kin%points_full)
  call loct_pointer_copy(kout%points_full_red, kin%points_full_red)
  call loct_pointer_copy(kout%weights_full, kin%weights_full)

  kout%nik = kin%nik
  call loct_pointer_copy(kout%points, kin%points)
  call loct_pointer_copy(kout%points_red, kin%points_red)
  call loct_pointer_copy(kout%weights, kin%weights)

  kout%use_symmetries = kin%use_symmetries
  kout%use_inversion  = kin%use_inversion

  kout%nik_axis = kin%nik_axis
  kout%shifts = kin%shifts

  call pop_sub()
end subroutine kpoints_copy

end module kpoints_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
