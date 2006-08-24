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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! ---------------------------------------------------------
subroutine states_choose_kpoints(d, sb, geo)
  type(states_dim_t), intent(inout) :: d
  type(simul_box_t),  intent(in)    :: sb
  type(geometry_t),   intent(in)    :: geo

  integer  :: coi, i, nkmax
  integer(POINTER_SIZE) :: blk
  FLOAT :: total_weight, kmax

  ! local variables for the crystal_init call
  logical :: full_ws_cell
  integer :: is, nk
  integer, allocatable :: natom(:)      ! natom(i) is the number of atoms of specie i
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

    ALLOCATE(d%kpoints (MAX_DIM, d%nik), MAX_DIM*d%nik)
    ALLOCATE(d%kweights   (d%nik),   d%nik)
    d%kpoints  = M_ZERO
    d%kweights = M_ONE

    call pop_sub()
    return
  end if

  !%Variable NumberKPoints
  !%Type block
  !%Default 1,1,1
  !%Section States
  !%Description
  !% A triplet of integers defining the number of kpoints to be used
  !% along each direction in the reciprocal space.
  !% The numbers refer to the whole BZ, and the actual number of kpoints
  !% is usually reduced exploiting the symmetries of the system
  !% For example, the following input samples the BZ with 100 points in the 
  !% xy plane of the reciprocal space
  !%
  !% <tt>%NumberKPoints
  !% <br>&nbsp;&nbsp;10 | 10 | 1
  !% <br>%</tt>
  !%
  !%End

  if(loct_parse_block(check_inp('NumberKPoints'), blk) .ne. 0) then
    message(1) = 'Block "NumberKPoints" not found in input file.'
    call write_fatal(1)
  end if

  d%nik_axis = 1
  do i = 1, sb%periodic_dim
    call loct_parse_block_int(blk, 0, i-1, d%nik_axis(i))
  end do
  call loct_parse_block_end(blk)
  if (any(d%nik_axis < 1)) then
    message(1) = 'Input: NumberKPoints is not valid'
    message(2) = '(NumberKPoints >= 1)'
    call write_fatal(2)
  end if
  nkmax = PRODUCT(d%nik_axis)

  !%Variable ShiftKpoints
  !%Type block
  !%Default 0.0,0.0,0.0
  !%Section States
  !%Description
  !% A triplet of real numbers to shift the Monkhorst-Pack k-points grid from its default position
  !%End

  if(loct_parse_block(check_inp('ShiftKPoints'), blk) .ne. 0) then
    kshifts = M_ZERO
  else
    do i = 1, sb%periodic_dim
      call loct_parse_block_float(blk, 0, i-1, kshifts(i))
    end do
    call loct_parse_block_end(blk)
  end if

  !%Variable FullWignerSeitzCell
  !%Type logical
  !%Default no
  !%Section States
  !%Description
  !% If true, no symmetry is taken into account for the choice of k-points.
  !% Instead the k-point sampling will range over the full Wigner-Seitz cell.
  !%End
  call loct_parse_logical(check_inp('FullWignerSeitzCell'), .false., full_ws_cell)

  !%Variable CenterOfInversion
  !%Type integer
  !%Default no
  !%Section States
  !%Description
  !% Only used in 1D periodic calculation to enforce the correspondig symmetry in the Brillouin Zone
  !%Option no 0
  !% The system has no center of inversion: use the whole BZ
  !%Option yes 1
  !% The system has a center of inversion: use half BZ
  !%End
  if(.not. full_ws_cell) then

    if(sb%periodic_dim == 1) then
      call loct_parse_int(check_inp('CenterOfInversion'), 0, coi)
      d%nik = d%nik_axis(1)/(1 + coi) + 1
      ALLOCATE(d%kpoints(3, d%nik), 3*d%nik)
      ALLOCATE(d%kweights  (d%nik),   d%nik)
      d%kpoints     = M_ZERO
      d%kweights    = M_ONE
      kmax          = (M_ONE - coi*M_HALF)*sb%klat(1,1)
      total_weight  = M_ZERO

      do i = 1, d%nik
        d%kpoints(1, i) = (i-1)*kmax/real(d%nik-1)
        if (i /= 1 .and. i /= d%nik) then
          d%kweights(i) = d%kweights(i) + coi
        end if
        total_weight = total_weight + d%kweights(i)
      end do
      d%kweights = d%kweights/total_weight
      if (d%ispin == 2) then
        message(1) = 'Not implemented yet.'
        call write_fatal(1)
      end if

      call pop_sub()
      return
    end if

  end if

  ALLOCATE(natom(geo%nspecies), geo%nspecies)
  ALLOCATE(coorat(geo%nspecies, geo%natoms, 3), geo%nspecies*geo%natoms*3)

  natom  = 0
  coorat = M_ZERO
  do i = 1, geo%natoms
    is = geo%atom(i)%spec%index
    natom(is) = natom(is) + 1
    coorat(is, natom(is), :) = geo%atom(i)%x(:)
  end do

  do i = 1, sb%dim
    coorat(:,:,i) = coorat(:,:,i) / sb%rlat(i, i) + M_HALF
  end do

  if (full_ws_cell) then
    ! find out how many points we have inside the Wigner-Seitz cell
    nkmax = points_inside_ws_cell(sb%periodic_dim, d%nik_axis, sb%klat) 

    write(message(1), '(a,i5)' ) 'Info: Total number of k-points inside Wigner-Seitz cell:', nkmax
    call write_info(1)
    
    ALLOCATE(kp(3, nkmax), 3*nkmax)
    ALLOCATE(kw(nkmax), nkmax)

    ! get kpoints
    call get_points_in_ws_cell(sb%periodic_dim, d%nik_axis, sb%klat, kp) 

    ! equal weights for all points
    nk = nkmax
    kw = M_ONE/nk
  else
    ALLOCATE(kp(3, nkmax), 3*nkmax)
    ALLOCATE(kw(nkmax), nkmax)

    ! choose k-points according to Monkhorst-Pack scheme
    ! all points are equal, but some are more equal than others :)
    call crystal_init(sb%rlat, geo%nspecies, natom, geo%natoms, coorat, d%nik_axis, &
      kshifts, nk, nkmax, kp, kw)
  end if

  ! double d%nik and copy points for spin polarized calc
  select case(d%ispin)
  case(1)
    d%nik = nk
    ALLOCATE(d%kpoints(3, d%nik), 3*d%nik)
    ALLOCATE(d%kweights  (d%nik),   d%nik)
    do i = 1, 3
      d%kpoints(i,:) = kp(i,:)*sb%klat(i,i)
    end do
    d%kweights = kw
  case(2)
    d%nik = 2 * nk
    ALLOCATE(d%kpoints(3, d%nik), 3*d%nik)
    ALLOCATE(d%kweights  (d%nik),   d%nik)
    do i = 1,3
      d%kpoints(i,::2)  = kp(i,:)*sb%klat(i,i)
      d%kpoints(i,2::2) = kp(i,:)*sb%klat(i,i)
    end do
    d%kweights(::2)  = kw(:)
    d%kweights(2::2) = kw(:)
  end select

  deallocate(natom, coorat)
  deallocate(kp, kw)

  call pop_sub()
end subroutine states_choose_kpoints


! ---------------------------------------------------------
subroutine kpoints_write_info(d, iunit)
  type(states_dim_t), intent(in) :: d
  integer, intent(in) :: iunit

  integer :: ik

  call push_sub('states_kpoints.kpoints_write_info')

  write(message(1),'(a,3(i3,1x))') 'Number of k points in each direction = ', d%nik_axis(:)
  call write_info(1, iunit)

  write(message(1), '(3x,a,7x,a,9x,a,9x,a,8x,a)') 'ik', 'K_x', 'K_y', 'K_z', 'Weight'
  message(2) = '       --------------------------------------------------'
  call write_info(2, iunit, verbose_limit=80)

  do ik = 1, d%nik
     write(message(1),'(i4,1x,4f12.4)') ik,d%kpoints(:,ik)*units_out%length%factor, d%kweights(ik)
     call write_info(1, iunit, verbose_limit=80)
  end do

  call pop_sub()
end subroutine kpoints_write_info


! ---------------------------------------------------------
subroutine get_points_in_ws_cell(dim, nik_axis, klat, kp) 
  integer, intent(in) :: dim
  integer, intent(in) :: nik_axis(:)
  FLOAT,   intent(in) :: klat(:,:)
  FLOAT,  intent(out) :: kp(:,:)

  integer :: ix, iy, iz, inik, iunit, ik(3), id
  FLOAT   :: k_point(3)

  call push_sub('states_kpoints.get_points_in_ws_cell')

  iunit = io_open('static/wigner_seitz_cell.dat', action = 'write')
  write(message(1), '(a)') '# index  k_x, k_y, k_z (unscaled)  k_x, k_y, k_z (scaled)'
  call write_info(1, iunit)

  inik = 1

  do ix = 0, nik_axis(1) - 1
    do iy = 0, nik_axis(2) - 1
      do iz = 0, nik_axis(3) - 1

        k_point = M_ZERO
        ik = (/ ix, iy, iz /)

        do id = 1, dim
          k_point(id) = -klat(id, id) + ik(id) * ( M_TWO*klat(id, id) / (nik_axis(id) - 1) ) 
        end do

        if (in_wigner_seitz_cell(k_point, klat)) then
          ! return a scaled point (to be compatible with crystal_init)
          kp(:, inik) = M_ZERO
          do id = 1, dim
            kp(id, inik) = k_point(id)/klat(id, id)
          end do
          write(message(1),'(i4, 6f18.12)') inik, k_point(:), kp(:, inik)
          call write_info(1, iunit)
          inik = inik + 1          
        end if
        
      end do
    end do
  end do

  call io_close(iunit)

  call pop_sub()
end subroutine get_points_in_ws_cell


! ---------------------------------------------------------
integer function points_inside_ws_cell(dim, nik_axis, klat) result(no_of_points)
  integer, intent(in) :: dim
  integer, intent(in) :: nik_axis(:)
  FLOAT,   intent(in) :: klat(:,:)

  integer :: ix, iy, iz, inik, ik(3), id
  FLOAT   :: k_point(3)

  call push_sub('states_kpoints.points_inside_ws_cell')

  ! find out how many k-points are contained in the Wigner-Seitz cell

  inik = 0

  do ix = 0, nik_axis(1) - 1 
    do iy = 0, nik_axis(2) - 1
      do iz = 0, nik_axis(3) - 1

        k_point = M_ZERO
        ik = (/ ix, iy, iz /)

        do id = 1, dim
          k_point(id) = -klat(id, id) + ik(id) * ( M_TWO*klat(id, id) / (nik_axis(id) - 1) ) 
        end do

        if (in_wigner_seitz_cell(k_point, klat)) inik = inik + 1

      end do
    end do
  end do

  no_of_points = inik

  call pop_sub()
end function points_inside_ws_cell


! ---------------------------------------------------------
logical function in_wigner_seitz_cell(k_point, klat) result(in_cell)
  FLOAT, intent(in) :: k_point(:)
  FLOAT, intent(in) :: klat(:,:)

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
  asize = klat(1, 1)

  half_dist(:) = (/               &
    asize/M_TWO,                  &  ! cube
    asize/M_FOUR * sqrt(M_THREE), &  ! octahedron
    asize/M_FOUR * sqrt(  M_TWO)  &  ! dodecahedron
    /)
  
  !%Variable WignerSeitzCellType
  !%Type integer
  !%Default cube
  !%Section States
  !%Description
  !% Determines which form of Wigner-Seitz cell octopus should
  !% use for a full k-point sampling (currently this is not
  !% automatically determined from atomic positions and has to
  !% be specified manually)
  !%Option cube 1
  !% The Wigner-Seitz cell is has simple cubic form
  !%Option octahedron 2
  !% The Wigner-Seitz cell is has the form of a octahedron
  !%Option dodecahedron 3
  !% The Wigner-Seitz cell is has the form of a octahedron
  !%End
  call loct_parse_int(check_inp('WignerSeitzCellType'), CUBE, ws_cell_type)

  ! the number of Bragg planes corresponds to the coordination number
  ! of the respective Bravais lattice
  select case(ws_cell_type)
  case(CUBE)

    ALLOCATE(bragg_normal(3, 6), 3*6)

    bragg_normal(1:3, 1) = (/  2,   0,   0/) 
    bragg_normal(1:3, 2) = (/ -2,   0,   0/) 
    bragg_normal(1:3, 3) = (/  0,   2,   0/) 
    bragg_normal(1:3, 4) = (/  0,  -2,   0/) 
    bragg_normal(1:3, 5) = (/  0,   0,   2/) 
    bragg_normal(1:3, 6) = (/  0,   0,  -2/) 

  case(OCTAHEDRON)

    ALLOCATE(bragg_normal(3, 8), 3*8)

    bragg_normal(1:3, 1) = (/ -1,  -1,  -1 /) 
    bragg_normal(1:3, 2) = (/ -1,   1,  -1 /) 
    bragg_normal(1:3, 3) = (/  1,  -1,  -1 /) 
    bragg_normal(1:3, 4) = (/  1,   1,  -1 /) 
    bragg_normal(1:3, 5) = (/ -1,  -1,   1 /) 
    bragg_normal(1:3, 6) = (/ -1,   1,   1 /) 
    bragg_normal(1:3, 7) = (/  1,  -1,   1 /) 
    bragg_normal(1:3, 8) = (/  1,   1,   1 /) 
    
  case(DODECAHEDRON)

    ALLOCATE(bragg_normal(3, 12), 3*12)

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
  
  deallocate(bragg_normal)

  call pop_sub()
end function in_wigner_seitz_cell
