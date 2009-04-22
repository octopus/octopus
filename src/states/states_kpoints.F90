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

  call push_sub('states_kpoints.states_choose_kpoints')

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

    call pop_sub()
    return
  end if

  !%Variable KPoints
  !%Type block
  !%Section Mesh::KPoints
  !%Description
  !% This block defines an explicit set of k-points and their weights for
  !% a periodic-system calculation. The first column is the weight
  !% of each k-point and the following are the components of the k-point
  !% vector. You only need to specify the components for the
  !% periodic directions. Note that the k-points should be given in
  !% reciprocal-space coordinates (not in reduced coordinates), i.e.
  !% what Octopus writes in a line in the standard output as
  !% <tt>#k =   1, k = (    0.154000,    0.154000,    0.154000)</tt>.
  !%
  !% For example, if you want to include only the gamma point, you can
  !% use:
  !%
  !% <tt>%KPoints
  !% <br>&nbsp;&nbsp;1.0 | 0 | 0 | 0
  !% <br>%</tt>
  !%
  !%End

  if(loct_parse_block(datasets_check('KPoints'), blk) == 0) then

    if (d%ispin == 2) then
      message(1) = 'Not implemented yet.'
      call write_fatal(1)
    end if

    !we have a block with the k-points given explicitly
    d%nik = loct_parse_block_n(blk)
    write(message(1), '(a,i4,a)') 'Input: ', d%nik, ' k-points will be read from the input file'
    call write_info(1)

    SAFE_ALLOCATE( d%kpoints(1:MAX_DIM, 1:d%nik))
    SAFE_ALLOCATE(d%kweights(1:d%nik))

    d%kpoints = M_ZERO

    do ik = 1, d%nik
      call loct_parse_block_float(blk, ik - 1, 0, d%kweights(ik))
      do idim = 1, sb%periodic_dim
        call loct_parse_block_float(blk, ik - 1, idim, d%kpoints(idim, ik))
      end do
    end do

    d%kpoints = d%kpoints / units_inp%length%factor !K points have 1/length units

    call print_kpoints_debug
    call pop_sub()
    return
  end if

  !%Variable KPointsMonkhorstPack
  !%Type block
  !%Default 1,1,1
  !%Section Mesh::KPoints
  !%Description
  !% When this block is given (and the KPoints block is not present),
  !% k-points are arranged in a Monkhorst-Pack grid.
  !%
  !% The first row of the block is a triplet of integers defining the
  !% number of k-points to be used along each direction in the
  !% reciprocal space. The numbers refer to the whole Brillouin zone,
  !% and the actual number of k-points is usually reduced exploiting
  !% the symmetries of the system.  An optional second row can specify
  !% a shift in the k-points.
  !%
  !% For example, the following input samples the BZ with 100 points in the 
  !% xy plane of the reciprocal space:
  !%
  !% <tt>%KPointsMP
  !% <br>&nbsp;&nbsp;10 | 10 | 1
  !% <br>%</tt>
  !%
  !%End

  ! default when nothing specified: gamma point only
  if(loct_parse_block(datasets_check('KPointsMonkhorstPack'), blk) .ne. 0) then

    if (d%ispin == 2) then
      message(1) = 'Not implemented yet.'
      call write_fatal(1)
    end if

    d%nik = 1

    SAFE_ALLOCATE( d%kpoints(1:MAX_DIM, 1:d%nik))
    SAFE_ALLOCATE(d%kweights(1:d%nik))

    d%kweights(1) = M_ONE
    d%kpoints(1:MAX_DIM, 1) = M_ZERO

    call print_kpoints_debug
    call pop_sub()
    return
  end if

! now deal with Monkhorst-Pack
  d%nik_axis = 1
  do i = 1, sb%periodic_dim
    call loct_parse_block_int(blk, 0, i-1, d%nik_axis(i))
  end do
  if (any(d%nik_axis < 1)) then
    message(1) = 'Input: KPointsMonkhorstPack is not valid'
    call write_fatal(1)
  end if
  nkmax = PRODUCT(d%nik_axis)

  kshifts = M_ZERO

  if(loct_parse_block_n(blk) > 1) then ! we have a shift
    do i = 1, sb%periodic_dim
      call loct_parse_block_float(blk, 1, i-1, kshifts(i))
    end do
  end if
  call loct_parse_block_end(blk)

  !%Variable KPointsUseSymmetries
  !%Type logical
  !%Default yes
  !%Section Mesh::KPoints
  !%Description
  !% This variable defines whether symmetries are taken into account
  !% or not for the choice of k-points. If it is set to no, the k-point
  !% sampling will range over the full Brillouin zone. The default is
  !% yes.
  !%End
  call loct_parse_logical(datasets_check('KPointsUseSymmetries'), .true., use_symmetries)

  !%Variable KPointsUseTimeReversal
  !%Type logical
  !%Default yes
  !%Section Mesh::KPoints
  !%Description
  !% WARNING: This variable does not seem to work.
  !% This variable defines whether time-reversal symmetry is taken into account
  !% or not for the choice of k-points. If it is set to no, the k-point
  !% sampling will not be reduced according to time-reversal symmetry. The default is
  !% yes. If KPointsUseSymmetries = no, this variable is ignored, and time-reversal
  !% symmetry is not used.
  !%End
  call loct_parse_logical(datasets_check('KPointsUseTimeReversal'), .true., use_time_reversal)
  

  SAFE_ALLOCATE(natom(1:geo%nspecies))
  SAFE_ALLOCATE(coorat(1:geo%nspecies, 1:geo%natoms, 1:3))

  natom  = 0
  coorat = M_ZERO
  do i = 1, geo%natoms
    is = geo%atom(i)%spec%index
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
    integer :: iunit, ik

    if(in_debug_mode) then
      
      iunit = io_open('debug/kpoints', action = 'write')
      
      do ik = 1, d%nik
        write(iunit, '(4e30.20)') d%kweights(ik), d%kpoints(:, ik)
      end do
      
      call io_close(iunit)

    end if

  end subroutine print_kpoints_debug

end subroutine states_choose_kpoints


! ---------------------------------------------------------
subroutine kpoints_write_info(d, sbdim, iunit)
  type(states_dim_t), intent(in) :: d
  integer,            intent(in) :: sbdim
  integer,            intent(in) :: iunit

  integer :: ik

  call push_sub('states_kpoints.kpoints_write_info')

  write(message(1),'(a,9(i3,1x))') 'Number of k points in each direction = ', d%nik_axis(1:sbdim)
  call write_info(1, iunit)

  write(message(1), '(3x,a,7x,a,9x,a,9x,a,8x,a)') 'ik', 'K_x', 'K_y', 'K_z', 'Weight'
  message(2) = '       --------------------------------------------------'
  call write_info(2, iunit, verbose_limit=80)

  do ik = 1, d%nik
     write(message(1),'(i4,1x,4f12.4)') &
       ik, d%kpoints(1:sbdim,ik)*units_out%length%factor, d%kweights(ik)
     call write_info(1, iunit, verbose_limit=80)
  end do

  call pop_sub()
end subroutine kpoints_write_info

! ---------------------------------------------------------
logical function in_wigner_seitz_cell(k_point, klattice) result(in_cell)
  FLOAT, intent(in) :: k_point(:)
  FLOAT, intent(in) :: klattice(:,:)

  integer, allocatable :: bragg_normal(:,:)
  integer :: ib, ws_cell_type
  FLOAT   :: half_dist(3), normal(3), asize

  ! Wigner Seitz cells of cubic lattices
  integer, parameter :: &
    CUBE         = 1,   &
    OCTAHEDRON   = 2,   &
    DODECAHEDRON = 3

  call push_sub('states_kpoints.in_wigner_seitz_cell')

  ! Warning: currently the routine allows only for cubic systems
  ! Tetragonal, orthorombic and monoclinic systems are excluded
  asize = klattice(1, 1)

  half_dist(:) = (/               &
    asize/M_TWO,                  &  ! cube
    asize/M_FOUR * sqrt(M_THREE), &  ! octahedron
    asize/M_FOUR * sqrt(  M_TWO)  &  ! dodecahedron
    /)
  
  !%Variable KPointsCellType
  !%Type integer
  !%Default cube
  !%Section Mesh::KPoints
  !%Description
  !% Determines which shape of Wigner-Seitz cell Octopus should use
  !% for a full k-point sampling (currently this is not automatically
  !% determined from atomic positions and has to be specified
  !% manually). By default a cube is chosen.
  !%Option cube 1
  !% The has is simple cubic shape.
  !%Option octahedron 2
  !% The cell is octahedral.
  !%Option dodecahedron 3
  !% The cell is dodecahedral.
  !%End
  call loct_parse_int(datasets_check('KPointsCellType'), CUBE, ws_cell_type)

  ! the number of Bragg planes corresponds to the coordination number
  ! of the respective Bravais lattice
  select case(ws_cell_type)
  case(CUBE)

    SAFE_ALLOCATE(bragg_normal(1:3, 1:6))

    bragg_normal(1:3, 1) = (/  2,   0,   0/) 
    bragg_normal(1:3, 2) = (/ -2,   0,   0/) 
    bragg_normal(1:3, 3) = (/  0,   2,   0/) 
    bragg_normal(1:3, 4) = (/  0,  -2,   0/) 
    bragg_normal(1:3, 5) = (/  0,   0,   2/) 
    bragg_normal(1:3, 6) = (/  0,   0,  -2/) 

  case(OCTAHEDRON)

    SAFE_ALLOCATE(bragg_normal(1:3, 1:8))

    bragg_normal(1:3, 1) = (/ -1,  -1,  -1 /) 
    bragg_normal(1:3, 2) = (/ -1,   1,  -1 /) 
    bragg_normal(1:3, 3) = (/  1,  -1,  -1 /) 
    bragg_normal(1:3, 4) = (/  1,   1,  -1 /) 
    bragg_normal(1:3, 5) = (/ -1,  -1,   1 /) 
    bragg_normal(1:3, 6) = (/ -1,   1,   1 /) 
    bragg_normal(1:3, 7) = (/  1,  -1,   1 /) 
    bragg_normal(1:3, 8) = (/  1,   1,   1 /) 
    
  case(DODECAHEDRON)

    SAFE_ALLOCATE(bragg_normal(1:3, 1:12))

    bragg_normal(1:3,  1) = (/  0,  -1,  -1 /) 
    bragg_normal(1:3,  2) = (/  0,  -1,   1 /) 
    bragg_normal(1:3,  3) = (/ -1,  -1,   0 /) 
    bragg_normal(1:3,  4) = (/  1,  -1,   0 /) 
    bragg_normal(1:3,  5) = (/  0,   1,  -1 /) 
    bragg_normal(1:3,  6) = (/  0,   1,   1 /) 
    bragg_normal(1:3,  7) = (/ -1,   1,   0 /) 
    bragg_normal(1:3,  8) = (/  1,   1,   0 /) 
    bragg_normal(1:3,  9) = (/ -1,   0,  -1 /) 
    bragg_normal(1:3, 10) = (/ -1,   0,   1 /) 
    bragg_normal(1:3, 11) = (/  1,   0,  -1 /) 
    bragg_normal(1:3, 12) = (/  1,   0,   1 /) 

  end select


  in_cell = .true.

  ! loop over all Bragg planes  
  do ib = 1, size(bragg_normal, 2) 

    ! compute normalized normal vector of Bragg plane 
    normal = bragg_normal(:, ib) / &
      sqrt(dot_product(real(bragg_normal(:, ib)), real(bragg_normal(:, ib)))) 

    ! if the projection of the current k-point onto the normal of the 
    ! Bragg plane is larger than half_dist the k-point lies outside the WS cell 
    if (dot_product(k_point, normal) .gt. half_dist(ws_cell_type)) then
      in_cell = .false. 
    end if

  end do
  
  SAFE_DEALLOCATE_A(bragg_normal)

  call pop_sub()
end function in_wigner_seitz_cell

logical pure function kpoint_is_gamma(this, ik)
  type(states_dim_t), intent(in) :: this
  integer,            intent(in) :: ik
  
  kpoint_is_gamma = (maxval(abs(this%kpoints(:, ik))) < M_EPSILON)

end function kpoint_is_gamma

integer function kpoint_index(this, ik) result(index)
  type(states_dim_t), intent(in) :: this
  integer,            intent(in) :: ik

  if (this%nspin == 2) then
    if (states_dim_get_spin_index(this, ik) == 1) then
      index = (ik + 1) / 2
    else
      index = ik / 2
    endif
  else
     index = ik
  endif

end function kpoint_index

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
