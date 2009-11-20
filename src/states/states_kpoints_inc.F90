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
!! $Id$

! ---------------------------------------------------------
subroutine states_choose_kpoints(d, sb, geo)
  type(states_dim_t), intent(inout) :: d
  type(simul_box_t),  intent(in)    :: sb
  type(geometry_t),   intent(in)    :: geo

  integer   :: i, nkmax, ik, idim
  type(block_t) :: blk

  ! local variables for the crystal_init call
  logical :: use_symmetries, use_time_reversal
  integer :: is, nk
  integer, allocatable :: natom(:)      ! natom(i) is the number of atoms of species i
  FLOAT :: kshifts(MAX_DIM)
  FLOAT, allocatable :: coorat(:,:,:)   ! coorat(i,j,k) is the k-th component (lattice coordinates)
  ! of the position of the j-th atom of type i.
  FLOAT, allocatable :: kp(:,:),kw(:)

  call push_sub('states_kpoints_inc.states_choose_kpoints')

  ! if not periodic just return the Gamma point
  if (sb%periodic_dim == 0) then
    select case(d%ispin)
    case(1,3)
      d%nik = 1
    case(2)
      d%nik = 2
    case default
      message(1) = 'Input: invalid SpinComponents'
      call write_fatal(1)
    end select

    SAFE_ALLOCATE(d%kpoints (1:MAX_DIM, 1:d%nik))
    SAFE_ALLOCATE(d%kweights(1:d%nik))
    d%kpoints  = M_ZERO
    d%kweights = M_ONE

    call pop_sub(); return
  end if

  ! if Monkhorst-Pack used, this variable will be reset
  d%nik_axis(:) = 0

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

  if(parse_block(datasets_check('KPoints'), blk) == 0) then

    if (d%ispin == 2) then
      message(1) = 'Not implemented yet.'
      call write_fatal(1)
    end if

    !we have a block with the k-points given explicitly
    d%nik = parse_block_n(blk)
    write(message(1), '(a,i4,a)') 'Input: ', d%nik, ' k-points will be read from the input file'
    call write_info(1)

    SAFE_ALLOCATE( d%kpoints(1:MAX_DIM, 1:d%nik))
    SAFE_ALLOCATE(d%kweights(1:d%nik))

    d%kpoints  = M_ZERO
    d%kweights = M_ZERO

    do ik = 1, d%nik
      call parse_block_float(blk, ik - 1, 0, d%kweights(ik))
      do idim = 1, sb%periodic_dim
        call parse_block_float(blk, ik - 1, idim, d%kpoints(idim, ik))
      end do
    end do

    d%kpoints = units_to_atomic(unit_one / units_inp%length, d%kpoints) !k-points have 1/length units

    call print_kpoints_debug
    call pop_sub(); return
  end if

  ! default when nothing specified: Gamma point only
  if(parse_block(datasets_check('KPointsGrid'), blk) .ne. 0) then

    if (d%ispin == 2) then
      message(1) = 'Not implemented yet.'
      call write_fatal(1)
    end if

    d%nik = 1

    SAFE_ALLOCATE(d%kpoints(1:MAX_DIM, 1:d%nik))
    SAFE_ALLOCATE(d%kweights(1:d%nik))

    d%kpoints     = M_ZERO
    d%kweights    = M_ZERO
    d%kweights(1) = M_ONE

    call print_kpoints_debug
    call pop_sub(); return
  end if

! now deal with Monkhorst-Pack
  d%nik_axis(:) = 1
  do i = 1, sb%periodic_dim
    call parse_block_integer(blk, 0, i-1, d%nik_axis(i))
  end do
  if (any(d%nik_axis < 1)) then
    message(1) = 'Input: KPointsMonkhorstPack is not valid'
    call write_fatal(1)
  end if
  nkmax = PRODUCT(d%nik_axis)

  kshifts = M_ZERO

  if(parse_block_n(blk) > 1) then ! we have a shift
    do i = 1, sb%periodic_dim
      call parse_block_float(blk, 1, i-1, kshifts(i))
    end do
  end if
  call parse_block_end(blk)

  !%Variable KPointsUseSymmetries
  !%Type logical
  !%Default no
  !%Section Mesh::KPoints
  !%Description
  !% This variable defines whether symmetries are taken into account
  !% or not for the choice of <i>k</i>-points. If it is set to no, the <i>k</i>-point
  !% sampling will range over the full Brillouin zone. The default is yes.
  !%End
  call parse_logical(datasets_check('KPointsUseSymmetries'), .false., use_symmetries)

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
  call parse_logical(datasets_check('KPointsUseTimeReversal'), .true., use_time_reversal)
  

  SAFE_ALLOCATE(natom(1:geo%nspecies))
  SAFE_ALLOCATE(coorat(1:geo%nspecies, 1:geo%natoms, 1:3))

  natom  = 0
  coorat = M_ZERO
  do i = 1, geo%natoms
    is = species_index(geo%atom(i)%spec)
    natom(is) = natom(is) + 1
    coorat(is, natom(is), 1:sb%dim) = geo%atom(i)%x(1:sb%dim)
  end do

  do i = 1, sb%dim
    coorat(:,:,i) = coorat(:,:,i) / sb%rlattice(i, i) + M_HALF
  end do

  SAFE_ALLOCATE(kp(1:3, 1:nkmax))
  SAFE_ALLOCATE(kw(1:nkmax))

  ! choose k-points according to Monkhorst-Pack scheme
  call crystal_init(geo, sb, d%nik_axis, kshifts, use_symmetries, use_time_reversal, nkmax, kp, kw)
  nk = nkmax

  ! double d%nik and copy points for spin-polarized calc
  select case(d%ispin)
  case(UNPOLARIZED, SPINORS)
    d%nik = nk
    SAFE_ALLOCATE( d%kpoints(1:3, 1:d%nik))
    SAFE_ALLOCATE(d%kweights(1:d%nik))
    do i = 1, 3
      d%kpoints(i, 1:d%nik) = kp(i, 1:d%nik)*sb%klattice(i,i)
    end do
    d%kweights(1:d%nik) = kw(1:d%nik)
  case(SPIN_POLARIZED)
    d%nik = 2 * nk
    SAFE_ALLOCATE( d%kpoints(1:3, 1:d%nik))
    SAFE_ALLOCATE(d%kweights(1:d%nik))
    do i = 1,3
      d%kpoints(i,::2)  = kp(i, 1:nk)*sb%klattice(i,i)
      d%kpoints(i,2::2) = kp(i, 1:nk)*sb%klattice(i,i)
    end do
    d%kweights(::2)  = kw(1:nk)
    d%kweights(2::2) = kw(1:nk)
  end select

  SAFE_DEALLOCATE_A(natom)
  SAFE_DEALLOCATE_A(coorat)
  SAFE_DEALLOCATE_A(kp)
  SAFE_DEALLOCATE_A(kw)
  
  call print_kpoints_debug
  call pop_sub()

contains
  subroutine print_kpoints_debug
    integer :: iunit

    call push_sub('states_kpoints_inc.states_choose_kpoints.print_kpoints_debug')

    if(in_debug_mode) then
      
      iunit = io_open('debug/kpoints', action = 'write')
      call kpoints_write_info(d, sb%dim, iunit)      
      call io_close(iunit)

    end if

    call pop_sub()
  end subroutine print_kpoints_debug

end subroutine states_choose_kpoints


! ---------------------------------------------------------
subroutine kpoints_write_info(d, sbdim, iunit)
  type(states_dim_t), intent(in) :: d
  integer,            intent(in) :: sbdim
  integer,            intent(in) :: iunit

  integer :: ik, idim

  call push_sub('states_kpoints_inc.kpoints_write_info')

  if(d%nik_axis(1) .ne. 0) then
    write(message(1),'(a,9(i3,1x))') 'Number of k-points in each direction = ', d%nik_axis(1:sbdim)
    call write_info(1, iunit)
  else
    ! a Monkhorst-Pack grid was not used
    write(message(1),'(a,9(i3,1x))') 'Number of k-points = ', kpoint_index(d, d%nik)
    call write_info(1, iunit)
  endif

  write(message(1), '(3x,a,7x,a,9x,a,9x,a,8x,a)') 'ik', 'K_x', 'K_y', 'K_z', 'Weight'
  message(2) = '   --------------------------------------------------'
  call write_info(2, iunit, verbose_limit=80)

  do ik = 1, d%nik
     write(message(1),'(i4,1x,4f12.4)') &
       ik, (units_from_atomic(unit_one/units_out%length, d%kpoints(idim,ik)), idim=1,sbdim), d%kweights(ik)
     call write_info(1, iunit, verbose_limit=80)
  end do

  call pop_sub()
end subroutine kpoints_write_info

logical pure function kpoint_is_gamma(this, ik)
  type(states_dim_t), intent(in) :: this
  integer,            intent(in) :: ik
  
  kpoint_is_gamma = (maxval(abs(this%kpoints(:, ik))) < M_EPSILON)

end function kpoint_is_gamma

integer function kpoint_index(this, ik) result(index)
  type(states_dim_t), intent(in) :: this
  integer,            intent(in) :: ik

  call push_sub('states_kpoints_inc.kpoint_index')

  if (this%nspin == 2) then
    if (states_dim_get_spin_index(this, ik) == 1) then
      index = (ik + 1) / 2
    else
      index = ik / 2
    endif
  else
     index = ik
  endif

  call pop_sub()
end function kpoint_index

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
