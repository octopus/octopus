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
  use math_m
  use messages_m
  use parser_m
  use profiling_m
  use symmetries_m
  use unit_m
  use unit_system_m
  use utils_m
  
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
    kpoints_set_point,       &
    kpoints_set_transport_mode,&
    kpoints_write_info,      &
    kpoints_point_is_gamma

  type kpoints_grid_t
    FLOAT, pointer :: point(:, :)
    FLOAT, pointer :: red_point(:, :)
    FLOAT, pointer :: weight(:)
    integer        :: npoints
    integer        :: dim
  end type kpoints_grid_t

  type kpoints_t
    type(kpoints_grid_t) :: full
    type(kpoints_grid_t) :: reduced

    integer        :: method

    logical        :: use_symmetries
    logical        :: use_time_reversal

    ! For the modified Monkhorst-Pack scheme
    integer        :: nik_axis(MAX_DIM)    !< number of MP divisions
    FLOAT          :: shifts(MAX_DIM)      ! 
    integer, pointer :: symmetry_ops(:, :)
    integer, pointer :: num_symmetry_ops(:)
  end type kpoints_t

  integer, parameter ::                &
    KPOINTS_GAMMA       =  1,          &
    KPOINTS_MONKH_PACK  =  2,          &
    KPOINTS_USER        =  3

contains

  subroutine kpoints_grid_init(dim, this, npoints)
    integer,              intent(in)  :: dim
    type(kpoints_grid_t), intent(out) :: this
    integer,              intent(in)  :: npoints

    PUSH_SUB(kpoints_grid_init)

    this%dim = dim
    this%npoints = npoints
    SAFE_ALLOCATE(this%red_point(1:dim, 1:npoints))
    SAFE_ALLOCATE(this%point(1:dim, 1:npoints))
    SAFE_ALLOCATE(this%weight(npoints))

    POP_SUB(kpoints_grid_init)
  end subroutine kpoints_grid_init

  ! ---------------------------------------------------------
  subroutine kpoints_grid_end(this)
    type(kpoints_grid_t), intent(out) :: this

    PUSH_SUB(kpoints_grid_end)

    SAFE_DEALLOCATE_P(this%red_point)
    SAFE_DEALLOCATE_P(this%point)
    SAFE_DEALLOCATE_P(this%weight)

    POP_SUB(kpoints_grid_end)
  end subroutine kpoints_grid_end


  ! ---------------------------------------------------------
  subroutine kpoints_grid_copy(bb, aa)
    type(kpoints_grid_t), intent(in)  :: bb
    type(kpoints_grid_t), intent(out) :: aa

    PUSH_SUB(kpoints_grid_copy)
    
    call kpoints_grid_init(bb%dim, aa, bb%npoints)
    
    aa%weight = bb%weight
    aa%point  = bb%point
    aa%red_point = bb%red_point

    POP_SUB(kpoints_grid_copy)
  end subroutine kpoints_grid_copy


  ! ---------------------------------------------------------
  subroutine kpoints_init(this, symm, dim, rlattice, klattice, only_gamma)
    type(kpoints_t),    intent(out) :: this
    type(symmetries_t), intent(in)  :: symm
    integer,            intent(in)  :: dim
    FLOAT,              intent(in)  :: rlattice(:,:), klattice(:,:)
    logical,            intent(in)  :: only_gamma

    integer :: ik, idir
    character(len=100) :: str_tmp

    PUSH_SUB(kpoints_init)

    ASSERT(dim <= MAX_DIM)

    nullify(this%symmetry_ops)
    nullify(this%num_symmetry_ops)

    !%Variable KPointsUseSymmetries
    !%Type logical
    !%Default no
    !%Section Mesh::KPoints
    !%Description
    !% This variable defines whether symmetries are taken into account
    !% or not for the choice of <i>k</i>-points. If it is set to no, the <i>k</i>-point
    !% sampling will range over the full Brillouin zone.
    !%
    !% The default is no.
    !%
    !% When a perturbation is applied to the system, the full
    !% symmetries of the system cannot be used. In this case you must
    !% not use symmetries or use the <tt>SymmetryBreakDir</tt> to tell
    !% Octopus the direction of the perturbation (for the moment this
    !% has to be done by hand by the user, in the future it will be
    !% automatic).
    !%
    !%End
    call parse_logical(datasets_check('KPointsUseSymmetries'), .false., this%use_symmetries)

    !%Variable KPointsUseTimeReversal
    !%Type logical
    !%Default yes
    !%Section Mesh::KPoints
    !%Description
    !% If symmetries are used to reduce the number of <i>k</i>-points,
    !% this variable defines whether time-reversal symmetry is taken
    !% into account or not. If it is set to no, the <i>k</i>-point
    !% sampling will not be reduced according to time-reversal
    !% symmetry.
    !%
    !% The default is yes, unless symmetries are broken in one
    !% direction by the SymmetryBreakDir block.
    !% 
    !% Warning: For time propagation runs with an external field,
    !% time-reversal symmetry should not be used.
    !%
    !%End
    call parse_logical(datasets_check('KPointsUseTimeReversal'), .not. symmetries_have_break_dir(symm), this%use_time_reversal)

    if(only_gamma) then
      this%method = KPOINTS_GAMMA
      call read_MP(gamma_only = .true.)
    else 
      if(read_user_kpoints()) then
        this%method = KPOINTS_USER
      else
        this%method = KPOINTS_MONKH_PACK
        call read_MP(gamma_only = .false.)

        write(message(1),'(a)') ' '
        write(message(2),'(1x,i3,a)') this%reduced%npoints, ' k-points generated from parameters :'
        write(message(3),'(1x,a)') '---------------------------------------------------'

        write(message(4),'(4x,a)') 'n ='    
        do idir = 1, dim
          write(str_tmp,'(i5)') this%nik_axis(idir)
          message(4) = trim(message(4)) // trim(str_tmp)
        enddo
        write(str_tmp,'(6x,a)') 's ='
        message(4) = trim(message(4)) // trim(str_tmp)
        do idir = 1, dim
          write(str_tmp,'(f6.2)') this%shifts(idir)
          message(4) = trim(message(4)) // trim(str_tmp)
        enddo
        call messages_info(4)

      end if

      write(message(1),'(a)') ' '
      write(message(2),'(a)') ' index |    weight    |             coordinates              |'
      call messages_info(2)

      do ik = 1, this%reduced%npoints
        write(str_tmp,'(i6,a,f12.6,a)') ik, " | ", this%reduced%weight(ik), " |"
        message(1) =  str_tmp
        do idir = 1, dim
          write(str_tmp,'(f12.6)') this%reduced%red_point(idir, ik)
          message(1) = trim(message(1)) // trim(str_tmp)
        enddo
        write(str_tmp,'(a)') "  |"
        message(1) = trim(message(1)) // trim(str_tmp)
        call messages_info(1)
      end do

      write(message(1),'(a)') ' '
      call messages_info(1)

    end if

    POP_SUB(kpoints_init)

  contains

    ! ---------------------------------------------------------
    subroutine read_MP(gamma_only)
      logical, intent(in) :: gamma_only

      logical       :: gamma_only_
      integer       :: ii, ncols
      type(block_t) :: blk
      integer, allocatable :: symm_ops(:, :), num_symm_ops(:)
      

      PUSH_SUB(kpoints_init.read_MP)

      call messages_obsolete_variable('KPointsMonkhorstPack', 'KPointsGrid')

      !%Variable KPointsGrid
      !%Type block
      !%Default Gamma-point only
      !%Section Mesh::KPoints
      !%Description
      !% When this block is given (and the <tt>KPoints</tt> block is not present),
      !% <i>k</i>-points are distributed in a uniform grid.
      !%
      !% The first row of the block is a set of integers defining
      !% the number of <i>k</i>-points to be used along each direction
      !% in reciprocal space. The numbers refer to the whole Brillouin
      !% zone, and the actual number of <i>k</i>-points is usually
      !% reduced exploiting the symmetries of the system.  By default
      !% the grid will always include the Gamma point. An optional
      !% second row can specify a shift in the <i>k</i>-points.
      !% The number of columns should be equal to <tt>Dimensions</tt>,
      !% but the grid and shift numbers should be 1 and zero in finite directions.
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
      if(.not. gamma_only_) &
        gamma_only_ = (parse_block(datasets_check('KPointsGrid'), blk) .ne. 0)

      this%nik_axis(1:MAX_DIM) = 1
      this%shifts(1:MAX_DIM) = M_ZERO

      if(.not. gamma_only_) then
        ncols = parse_block_cols(blk, 0)
        if(ncols /= dim) then
          write(message(1),'(a,i3,a,i3)') 'KPointsGrid first row has ', ncols, ' columns but must have ', dim
          call messages_fatal(1)
        endif
        do ii = 1, dim
          call parse_block_integer(blk, 0, ii - 1, this%nik_axis(ii))
        end do

        if (any(this%nik_axis(1:dim) < 1)) then
          message(1) = 'Input: KPointsGrid is not valid.'
          call messages_fatal(1)
        end if

        if(parse_block_n(blk) > 1) then ! we have a shift
          ncols = parse_block_cols(blk, 1)
          if(ncols /= dim) then
            write(message(1),'(a,i3,a,i3)') 'KPointsGrid second row has ', ncols, ' columns but must have ', dim
            call messages_fatal(1)
          endif
          do ii = 1, dim
            call parse_block_float(blk, 1, ii - 1, this%shifts(ii))
          end do
        end if

        call parse_block_end(blk)
      end if

      call kpoints_grid_init(dim, this%full, product(this%nik_axis(1:dim)))

      call kpoints_grid_generate(dim, this%nik_axis(1:dim), this%shifts(1:dim), this%full%red_point)

      this%full%weight = M_ONE / this%full%npoints

      call kpoints_grid_copy(this%full, this%reduced)


      if(this%use_symmetries) then

        SAFE_ALLOCATE(num_symm_ops(1:this%full%npoints))
        SAFE_ALLOCATE(symm_ops(1:this%full%npoints, 1:symmetries_number(symm)))
        
        call kpoints_grid_reduce(symm, this%use_time_reversal, &
          this%reduced%npoints, dim, this%reduced%red_point, this%reduced%weight, symm_ops, num_symm_ops)
        
        ASSERT(maxval(num_symm_ops) >= 0)
        ASSERT(maxval(num_symm_ops) <= symmetries_number(symm))
        
        SAFE_ALLOCATE(this%num_symmetry_ops(1:this%reduced%npoints))
        SAFE_ALLOCATE(this%symmetry_ops(1:this%reduced%npoints, 1:maxval(num_symm_ops)))
        
        this%num_symmetry_ops(1:this%reduced%npoints) = num_symm_ops(1:this%reduced%npoints)
        this%symmetry_ops(1:this%reduced%npoints, 1:maxval(num_symm_ops)) = &
          symm_ops(1:this%reduced%npoints, 1:maxval(num_symm_ops))
        
        SAFE_DEALLOCATE_A(num_symm_ops)
        SAFE_DEALLOCATE_A(symm_ops)
      end if
      
      do ik = 1, this%full%npoints
        call kpoints_to_absolute(klattice, this%full%red_point(:, ik), this%full%point(:, ik), dim)
      end do

      do ik = 1, this%reduced%npoints
        call kpoints_to_absolute(klattice, this%reduced%red_point(:, ik), this%reduced%point(:, ik), dim)
      end do

      POP_SUB(kpoints_init.read_MP)
    end subroutine read_MP


    ! ---------------------------------------------------------
    logical function read_user_kpoints()
      type(block_t) :: blk
      logical :: reduced
      integer :: ik, idir

      PUSH_SUB(kpoints_init.read_user_kpoints)

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
      !% what <tt>Octopus</tt> writes in a line in the ground-state standard output as
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
          POP_SUB(kpoints_init.read_user_kpoints)
          return
        end if
      end if
      read_user_kpoints = .true.

      call kpoints_grid_init(dim, this%full, parse_block_n(blk))

      this%full%red_point = M_ZERO
      this%full%point = M_ZERO
      this%full%weight = M_ZERO

      if(reduced) then
        do ik = 1, this%full%npoints
          call parse_block_float(blk, ik - 1, 0, this%full%weight(ik))
          do idir = 1, dim
            call parse_block_float(blk, ik - 1, idir, this%full%red_point(idir, ik))
          end do
          ! generate also the absolute coordinates
          call kpoints_to_absolute(klattice, this%full%red_point(:, ik), this%full%point(:, ik), dim)
        end do
      else
        do ik = 1, this%full%npoints
          call parse_block_float(blk, ik - 1, 0, this%full%weight(ik))
          do idir = 1, dim
            call parse_block_float(blk, ik - 1, idir, this%full%point(idir, ik), unit_one/units_inp%length)
          end do
          ! generate also the reduced coordinates
          call kpoints_to_reduced(rlattice, this%full%point(:, ik), this%full%red_point(:, ik), dim)
        end do
      end if
      call parse_block_end(blk)

      ! for the moment we do not apply symmetries to user kpoints
      call kpoints_grid_copy(this%full, this%reduced)

      write(message(1), '(a,i4,a)') 'Input: ', this%full%npoints, ' k-points were read from the input file'
      call messages_info(1)

      POP_SUB(kpoints_init.read_user_kpoints)
    end function read_user_kpoints

  end subroutine kpoints_init

  ! ---------------------------------------------------------
  subroutine kpoints_end(this)
    type(kpoints_t), intent(inout) :: this

    PUSH_SUB(kpoints_end)

    call kpoints_grid_end(this%full)
    call kpoints_grid_end(this%reduced)

    SAFE_DEALLOCATE_P(this%symmetry_ops)
    SAFE_DEALLOCATE_P(this%num_symmetry_ops)

    POP_SUB(kpoints_end)
  end subroutine kpoints_end


  ! ---------------------------------------------------------
  !> sets the kpoints to zero (only to be used in transport mode)
  subroutine kpoints_set_transport_mode(this)
    type(kpoints_t), intent(inout) :: this

    PUSH_SUB(kpoints_set_transport_mode)

    this%full%point = M_ZERO
    this%full%red_point = M_ZERO
    this%reduced%point = M_ZERO
    this%reduced%red_point = M_ZERO

    POP_SUB(kpoints_set_transport_mode)
  end subroutine kpoints_set_transport_mode


  ! ---------------------------------------------------------
  subroutine kpoints_to_absolute(klattice, kin, kout, dim)
    FLOAT,   intent(in)  :: klattice(:,:)
    FLOAT,   intent(in)  :: kin(:)
    FLOAT,   intent(out) :: kout(:)
    integer, intent(in)  :: dim

    integer :: ii
    
    PUSH_SUB(kpoints_to_absolute)

    kout(1:dim) = M_ZERO
    do ii = 1, dim
      kout(1:dim) = kout(1:dim) + kin(ii) * klattice(1:dim, ii)
    end do

    POP_SUB(kpoints_to_absolute)
  end subroutine kpoints_to_absolute


  ! ---------------------------------------------------------
  subroutine kpoints_to_reduced(rlattice, kin, kout, dim)
    FLOAT,   intent(in)  :: rlattice(:,:)
    FLOAT,   intent(in)  :: kin(:)
    FLOAT,   intent(out) :: kout(:)
    integer, intent(in)  :: dim

    integer :: ii

    PUSH_SUB(kpoints_to_reduced)

    kout(1:dim) = M_ZERO
    do ii = 1, dim
      kout(1:dim) = kout(1:dim) + kin(ii) * rlattice(ii, 1:dim)
    end do
    kout(1:dim) = kout(1:dim) / (M_TWO * M_PI)

    POP_SUB(kpoints_to_reduced)
  end subroutine kpoints_to_reduced


  ! ---------------------------------------------------------
  subroutine kpoints_copy(kin, kout)
    type(kpoints_t), intent(in)  :: kin
    type(kpoints_t), intent(out) :: kout

    PUSH_SUB(kpoints_copy)

    kout%method = kin%method

    call kpoints_grid_copy(kin%full, kout%full)
    call kpoints_grid_copy(kin%reduced, kout%reduced)

    kout%use_symmetries = kin%use_symmetries
    kout%use_time_reversal = kin%use_time_reversal

    kout%nik_axis(1:kin%full%dim) = kin%nik_axis(1:kin%full%dim)
    kout%shifts  (1:kin%full%dim) = kin%shifts  (1:kin%full%dim)

    POP_SUB(kpoints_copy)
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
    FLOAT                       :: point(1:this%full%dim)

    point(1:this%full%dim) = this%reduced%point(1:this%full%dim, ik)

  end function kpoints_get_point

  ! ----------------------------------------------------------
  !> sets the ik-th kpoint to a given value (only to be used in transport mode)
  subroutine kpoints_set_point(this, ik, point)
    type(kpoints_t), intent(inout) :: this
    integer,         intent(in)    :: ik
    FLOAT  ,         intent(in)    :: point(1:this%full%dim)

    PUSH_SUB(kpoints_set_point)

    this%reduced%point(1:this%full%dim, ik) = point(1:this%full%dim)
    this%full%point(1:this%full%dim, ik)    = point(1:this%full%dim)

    POP_SUB(kpoints_set_point)
  end subroutine kpoints_set_point

  ! ----------------------------------------------------------
  FLOAT pure function kpoints_get_weight(this, ik) result(weight)
    type(kpoints_t), intent(in) :: this
    integer,         intent(in) :: ik

    weight = this%reduced%weight(ik)

  end function kpoints_get_weight


  ! ----------------------------------------------------------
  !> Generates the k-points grid.
  !! Sets up a uniform array of k-points. Use a modification of the normal Monkhorst-Pack scheme, 
  !! which is equivalent to the normal MP scheme in the case of even number of kpoints (i.e. naxis (i) even)  
  !! used with a shift of (1/2, 1/2, 1/2).
  !! For the original MP scheme, see (PRB 13, 518 (1976)) and (PRB 16, 1748 (1977))
  !! naxis(i) are the number of points in the three directions determined by the lattice vectors.
  !! shift(i) and sz shift the grid of integration points from the origin.
  subroutine kpoints_grid_generate(dim, naxis, shift, kpoints)  
    integer,           intent(in)  :: dim
    integer,           intent(in)  :: naxis(:)
    FLOAT,             intent(in)  :: shift(:)
    FLOAT,             intent(out) :: kpoints(:, :)
  
    FLOAT :: dx(1:MAX_DIM)
    integer :: ii, jj, divisor, ik, idir, ix(1:MAX_DIM), npoints
    FLOAT, allocatable :: nrm(:)

    PUSH_SUB(kpoints_grid_generate)
    
    dx(1:dim) = M_ONE/(M_TWO*naxis(1:dim))

    npoints = product(naxis(1:dim))

    SAFE_ALLOCATE(nrm(1:npoints))
    
    do ii = 0, npoints - 1

      ik = npoints - ii
      jj = ii
      divisor = npoints

      nrm(ik) = M_ZERO
      do idir = 1, dim
        divisor = divisor / naxis(idir)
        ix(idir) = jj / divisor + 1
        jj = mod(jj, divisor)

        kpoints(idir, ik) = (M_TWO*ix(idir) - M_ONE*naxis(idir) + M_TWO*shift(idir))*dx(idir)

        if(mod(naxis(idir), 2) /= 0) then    
          kpoints(idir, ik) = kpoints(idir, ik) - dx(idir)
        end if
  
        ! get the norm for sorting
        nrm(ik) = nrm(ik) + (kpoints(idir, ik)/dx(idir))**2

        ! tweak the norm a bit so the points have a specific order
        if(kpoints(idir, ik) >  M_EPSILON) nrm(ik) = nrm(ik) - CNST(0.01)*(dim - idir + 1)
        if(kpoints(idir, ik) < -M_EPSILON) nrm(ik) = nrm(ik) + CNST(0.01)*(dim - idir + 1)

        !bring back point to first Brillouin zone
        if ( kpoints(idir,ik) /= M_HALF )  kpoints(idir, ik) = mod(kpoints(idir, ik) + M_HALF, M_ONE) - M_HALF
      end do

    end do

    call sort(nrm, kpoints)

    POP_SUB(kpoints_grid_generate)
  end subroutine kpoints_grid_generate

  
  ! --------------------------------------------------------------------------------------------
  subroutine kpoints_grid_reduce(symm, time_reversal, nkpoints, dim, kpoints, weights, symm_ops, num_symm_ops)
    type(symmetries_t), intent(in)    :: symm
    logical,            intent(in)    :: time_reversal
    integer,            intent(inout) :: nkpoints
    integer,            intent(in)    :: dim
    FLOAT,              intent(inout) :: kpoints(:, :)
    FLOAT,              intent(out)   :: weights(:)
    integer,            intent(out)   :: symm_ops(:, :)
    integer,            intent(out)   :: num_symm_ops(:)

    integer :: nreduced
    FLOAT, allocatable :: reduced(:, :)
    
    FLOAT :: dw
    integer ik, iop, ik2
    FLOAT :: tran(MAX_DIM), tran_inv(MAX_DIM)
    integer, allocatable :: kmap(:)

    PUSH_SUB(kpoints_grid_reduce)

    ! reduce to irreducible zone

    ! kmap is used to mark reducible k-points and also to
    ! map reducible to irreducible k-points

    SAFE_ALLOCATE(kmap(1:nkpoints))
    SAFE_ALLOCATE(reduced(1:dim, 1:nkpoints))

    forall(ik = 1:nkpoints) kmap(ik) = ik
    
    dw = M_ONE / nkpoints

    nreduced = 0

    num_symm_ops = 0

    do ik = 1, nkpoints
      if (kmap(ik) /= ik) cycle
      
      ! new irreducible point
      ! mark reduced with negative kmap
      
      nreduced = nreduced + 1
      reduced(1:dim, nreduced) = kpoints(1:dim, ik)
      
      kmap(ik) = -nreduced
      
      weights(nreduced) = dw

      if (ik == nkpoints) cycle
      
      ! operate with the symmetry operations
      
      do iop = 1, symmetries_number(symm)
        call symmetries_apply_kpoint(symm, iop, reduced(1:dim, nreduced), tran)
        tran_inv(1:dim) = -tran(1:dim)
           
        ! remove (mark) k-points related to irreducible reduced by symmetry
        do ik2 = ik + 1, nkpoints
          if (kmap(ik2) /= ik2) cycle
          
          ! both the transformed rk ...
          if (all( abs(tran(1:dim) - kpoints(1:dim, ik2)) <= CNST(1.0e-5))) then
            weights(nreduced) = weights(nreduced) + dw
            kmap(ik2) = nreduced
            INCR(num_symm_ops(nreduced), 1)
            symm_ops(nreduced, num_symm_ops(nreduced)) = iop
            cycle
          end if

          ! and its inverse
          if (time_reversal .and. all(abs(tran_inv(1:dim) - kpoints(1:dim, ik2)) <= CNST(1.0e-5)) ) then
            weights(nreduced) = weights(nreduced) + dw
            kmap(ik2) = nreduced
            INCR(num_symm_ops(nreduced), 1)
            symm_ops(nreduced, num_symm_ops(nreduced)) = iop
          end if
          
        end do
      end do
    end do
    
    nkpoints = nreduced
    do ik = 1, nreduced
      kpoints(1:dim, ik) = reduced(1:dim, ik)
    end do

    POP_SUB(kpoints_grid_reduce)
  end subroutine kpoints_grid_reduce


  ! ---------------------------------------------------------
  subroutine kpoints_write_info(this, iunit)
    type(kpoints_t),    intent(in) :: this
    integer,            intent(in) :: iunit
    
    integer :: ik, idir
    character(len=100) :: str_tmp
    character :: index
    
    PUSH_SUB(kpoints_write_info)
    
    if(this%method == KPOINTS_MONKH_PACK) then
      write(message(1),'(a)') 'Number of k-points in each direction = '
      do idir = 1, this%full%dim
        write(str_tmp,'(i3,1x)') this%nik_axis(idir)
        message(1) = trim(message(1)) // trim(str_tmp)
      enddo
      call messages_info(1, iunit)
    else
      ! a Monkhorst-Pack grid was not used
      write(message(1),'(a,i3)') 'Number of k-points = ', kpoints_number(this)
      call messages_info(1, iunit)
    endif

    write(message(1), '(6x,a)') 'ik'
    do idir = 1, this%full%dim
      index = index2axis(idir)
      write(str_tmp, '(9x,2a)') 'k_', index
      message(1) = trim(message(1)) // trim(str_tmp)
    enddo
    write(str_tmp, '(6x,a)') 'Weight'
    message(1) = trim(message(1)) // trim(str_tmp)
    message(2) = '---------------------------------------------------------'
    call messages_info(2, iunit, verbose_limit = .true.)
    
    do ik = 1, kpoints_number(this)
      write(message(1),'(i8,1x)') ik
      do idir = 1, this%full%dim
        write(str_tmp,'(f12.4)') units_from_atomic(unit_one/units_out%length, this%reduced%red_point(idir, ik))
        message(1) = trim(message(1)) // trim(str_tmp)
      enddo
      write(str_tmp,'(f12.4)') kpoints_get_weight(this, ik)
      message(1) = trim(message(1)) // trim(str_tmp)
      call messages_info(1, iunit, verbose_limit = .true.)
    end do
    
    POP_SUB(kpoints_write_info)
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
