!! Copyright (C) 2021 N. Tancogne-Dejean, M. Oliveira
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module lattice_vectors_oct_m
  use global_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private

  public ::                   &
    lattice_vectors_t,        &
    reciprocal_lattice

  type lattice_vectors_t
    ! Components are public by default
    type(space_t), private :: space
    FLOAT :: rlattice_primitive(MAX_DIM,MAX_DIM)   !< lattice primitive vectors
    FLOAT :: rlattice          (MAX_DIM,MAX_DIM)   !< lattice vectors
    FLOAT :: klattice_primitive(MAX_DIM,MAX_DIM)   !< reciprocal-lattice primitive vectors
    FLOAT :: klattice          (MAX_DIM,MAX_DIM)   !< reciprocal-lattice vectors
    FLOAT :: alpha, beta, gamma                    !< the angles defining the cell
    FLOAT :: rcell_volume                          !< the volume of the cell defined by the lattice vectors in real spac
    logical :: nonorthogonal = .false.
  contains
    procedure :: copy => lattice_vectors_copy
    generic   :: assignment(=) => copy
    procedure :: scale => lattice_vectors_scale
    procedure :: write_info => lattice_vectors_write_info
  end type lattice_vectors_t

  interface lattice_vectors_t
    module procedure lattice_vectors_constructor
  end interface lattice_vectors_t

contains

  !--------------------------------------------------------------
  type(lattice_vectors_t) function lattice_vectors_constructor(namespace, space) result(latt)
    type(namespace_t), intent(in)    :: namespace
    type(space_t),     intent(in)    :: space

    type(block_t) :: blk
    FLOAT :: norm, lparams(3), volume_element, rlatt(MAX_DIM, MAX_DIM)
    integer :: idim, jdim, ncols
    logical :: has_angles
    FLOAT :: angles(1:MAX_DIM)

    PUSH_SUB(lattice_vectors_constructor)

    latt%space = space

    latt%alpha = CNST(90.0)
    latt%beta  = CNST(90.0)
    latt%gamma = CNST(90.0)

    has_angles = .false.
    angles = CNST(90.0)

    if (space%is_periodic()) then

      !%Variable LatticeParameters
      !%Type block
      !%Section Mesh::Simulation Box
      !%Description
      !% The lattice parameters (a, b, c).
      !% This variable is mandatory for periodic systems and is ignored otherwise.
      !% When PeriodicDimensions = 3, a second optional line can be used to
      !% define the angles between the lattice vectors. If the angles are not
      !% provided, then the variable LatticeVectors must be set.
      !% The number of parameters specified in the block must be at least equal
      !% to the number of periodic dimensions, but it is not mandatory to
      !% specify parameters for the non-periodic dimensions (in that case they
      !% are set to 1).
      !%End
      if (parse_block(namespace, 'LatticeParameters', blk) == 0) then
        ncols = parse_block_cols(blk, 0) 
        if (ncols < space%periodic_dim) then
          call messages_input_error(namespace, 'LatticeParameters', 'The number of columns must be at least PeriodicDimensions')
        end if
        do idim = 1, ncols
          call parse_block_float(blk, 0, idim - 1, lparams(idim))
        end do

        ! If some parameters for non-periodic dimensions are not set in the input file, with set them to 1.
        do idim = ncols + 1, space%dim
          lparams(idim) = M_ONE
        end do

        ! Parse angles, if available
        if (parse_block_n(blk) > 1) then
          if (space%dim /= 3) then
            call messages_input_error(namespace, 'LatticeParameters', 'Angles can only be specified when Dimensions = 3')
          end if

          ncols = parse_block_cols(blk, 1)
          if (ncols /= space%dim) then
            call messages_input_error(namespace, 'LatticeParameters', 'You must specify three angles')
          end if
          do idim = 1, space%dim
            call parse_block_float(blk, 1, idim - 1, angles(idim))
          end do
          has_angles = .true.
        end if
        call parse_block_end(blk)

      else       
        call messages_input_error(namespace, 'LatticeParameters', 'Variable is mandatory for periodic systems')
      end if


      if (has_angles) then
        latt%alpha = angles(1)
        latt%beta  = angles(2)
        latt%gamma = angles(3)

        if (parse_is_defined(namespace, 'LatticeVectors')) then
          message(1) = 'LatticeParameters with angles is incompatible with LatticeVectors'
          call messages_print_var_info(stdout, "LatticeParameters")
          call messages_fatal(1, namespace=namespace)
        end if

        call build_metric_from_angles(latt, angles)

      else

        !%Variable LatticeVectors
        !%Type block
        !%Default simple cubic
        !%Section Mesh::Simulation Box
        !%Description
        !% Primitive lattice vectors. Vectors are stored in rows.
        !% Default:
        !% <br><br><tt>%LatticeVectors
        !% <br>&nbsp;&nbsp;1.0 | 0.0 | 0.0
        !% <br>&nbsp;&nbsp;0.0 | 1.0 | 0.0
        !% <br>&nbsp;&nbsp;0.0 | 0.0 | 1.0
        !% <br>%<br></tt>
        !%End
        latt%rlattice_primitive = M_ZERO
        latt%nonorthogonal = .false.
        do idim = 1, space%dim
          latt%rlattice_primitive(idim, idim) = M_ONE
        end do

        if (parse_block(namespace, 'LatticeVectors', blk) == 0) then
          do idim = 1, space%dim
            do jdim = 1, space%dim
              call parse_block_float(blk, idim - 1,  jdim - 1, latt%rlattice_primitive(jdim, idim))
              if (idim /= jdim .and. abs(latt%rlattice_primitive(jdim, idim)) > M_EPSILON) then
                latt%nonorthogonal = .true.
              end if
            enddo
          end do
          call parse_block_end(blk)

        end if
      end if

    else
      ! Non-periodic
      lparams = M_ONE
      latt%rlattice_primitive = M_ZERO
      do idim = 1, space%dim
        latt%rlattice_primitive(idim, idim) = M_ONE
      end do
    end if

    latt%rlattice = M_ZERO
    do idim = 1, space%dim
      norm = sqrt(sum(latt%rlattice_primitive(1:space%dim, idim)**2))
      lparams(idim) = lparams(idim)*norm
      do jdim = 1, space%dim
        latt%rlattice_primitive(jdim, idim) = latt%rlattice_primitive(jdim, idim) / norm
        latt%rlattice(jdim, idim) = latt%rlattice_primitive(jdim, idim) * lparams(idim)
      end do
    end do

    call reciprocal_lattice(latt%rlattice, latt%klattice, latt%rcell_volume, space%dim, namespace)
    latt%klattice = latt%klattice * M_TWO*M_PI

    call reciprocal_lattice(latt%rlattice_primitive, latt%klattice_primitive, volume_element, space%dim, namespace)

    ! rlattice_primitive is the A matrix from Chelikowski PRB 78 075109 (2008)
    ! klattice_primitive is the transpose (!) of the B matrix, with no 2 pi factor included
    ! klattice is the proper reciprocal lattice vectors, with 2 pi factor, and in units of 1/bohr
    ! The F matrix of Chelikowski is matmul(transpose(latt%klattice_primitive), latt%klattice_primitive)
    rlatt = matmul(transpose(latt%rlattice_primitive), latt%rlattice_primitive)
    if (.not. has_angles .and. space%dim == 3) then
      !We compute the angles from the lattice vectors
      latt%alpha = acos(rlatt(2,3)/sqrt(rlatt(2,2)*rlatt(3,3)))/M_PI*CNST(180.0)
      latt%beta  = acos(rlatt(1,3)/sqrt(rlatt(1,1)*rlatt(3,3)))/M_PI*CNST(180.0)
      latt%gamma = acos(rlatt(1,2)/sqrt(rlatt(1,1)*rlatt(2,2)))/M_PI*CNST(180.0)
    end if

    POP_SUB(lattice_vectors_constructor)
  end function lattice_vectors_constructor

  !--------------------------------------------------------------
  subroutine lattice_vectors_copy(this, source)
    class(lattice_vectors_t), intent(out) :: this
    class(lattice_vectors_t), intent(in)  :: source

    PUSH_SUB(lattice_vectors_copy)

    this%space = source%space
    this%rlattice_primitive = source%rlattice_primitive
    this%rlattice = source%rlattice
    this%klattice_primitive = source%klattice_primitive
    this%klattice = source%klattice
    this%alpha = source%alpha
    this%beta = source%beta
    this%gamma = source%gamma
    this%rcell_volume = source%rcell_volume
    this%nonorthogonal = source%nonorthogonal

    POP_SUB(lattice_vectors_copy)
  end subroutine lattice_vectors_copy

    !--------------------------------------------------------------
  subroutine lattice_vectors_scale(this, factor)
    class(lattice_vectors_t), intent(inout) :: this
    FLOAT,                    intent(in)    :: factor(this%space%dim)

    integer :: idir

    PUSH_SUB(lattice_vectors_scale)

    ! Scale the lattice in real space
    do idir = 1, this%space%dim
      this%rlattice(1:this%space%dim, idir) = this%rlattice(1:this%space%dim, idir)*factor(idir)
    end do

    ! Regenerate the lattice in reciprocal space
    call reciprocal_lattice(this%rlattice, this%klattice, this%rcell_volume, this%space%dim)
    this%klattice = this%klattice * M_TWO*M_PI
    
    POP_SUB(lattice_vectors_scale)
  end subroutine lattice_vectors_scale

  !--------------------------------------------------------------
  subroutine lattice_vectors_write_info(this, iunit)
    class(lattice_vectors_t), intent(in) :: this
    integer,                  intent(in) :: iunit

    integer :: idir, idir2

    PUSH_SUB(lattice_vectors_write_info)

    write(message(1),'(a,3a,a)') '  Lattice Vectors [', trim(units_abbrev(units_out%length)), ']'
    do idir = 1, this%space%dim
      write(message(1+idir),'(9f12.6)') (units_from_atomic(units_out%length, this%rlattice(idir2, idir)), &
        idir2 = 1, this%space%dim)
    end do
    call messages_info(1+this%space%dim, iunit)

    write(message(1),'(a,f18.4,3a,i1.1,a)') &
      '  Cell volume = ', units_from_atomic(units_out%length**this%space%dim, this%rcell_volume), &
      ' [', trim(units_abbrev(units_out%length**this%space%dim)), ']'
    call messages_info(1, iunit)

    write(message(1),'(a,3a,a)') '  Reciprocal-Lattice Vectors [', trim(units_abbrev(units_out%length**(-1))), ']'
    do idir = 1, this%space%dim
      write(message(1+idir),'(3f12.6)') (units_from_atomic(unit_one / units_out%length, this%klattice(idir2, idir)), &
        idir2 = 1, this%space%dim)
    end do
    call messages_info(1+this%space%dim, iunit)

    if (this%space%dim == 3) then
      write(message(1),'(a)') '  Cell angles [degree]'
      write(message(2),'(a, f8.3)') '    alpha = ', this%alpha
      write(message(3),'(a, f8.3)') '    beta  = ', this%beta
      write(message(4),'(a, f8.3)') '    gamma = ', this%gamma
      call messages_info(4, iunit)
    end if

    POP_SUB(lattice_vectors_write_info)
  end subroutine lattice_vectors_write_info

  !--------------------------------------------------------------
  subroutine build_metric_from_angles(this, angles)
    type(lattice_vectors_t), intent(inout) :: this
    FLOAT,                   intent(in)    :: angles(3)

    FLOAT :: cosang, a2, aa, cc
    FLOAT, parameter :: tol_angle = CNST(1.0e-6)

    PUSH_SUB(build_metric_from_angles)
   
    !Converting the angles to LatticeVectors
    !See 57_iovars/ingeo.F90 in Abinit for details
    if (abs(angles(1) - angles(2)) < tol_angle .and. abs(angles(2) - angles(3)) < tol_angle .and.  &
      (abs(angles(1) - CNST(90.0)) + abs(angles(2) - CNST(90.0)) + abs(angles(3) - CNST(90.0))) > tol_angle) then

      cosang = cos(M_PI*angles(1)/CNST(180.0))
      a2 = M_TWO/M_THREE*(M_ONE - cosang)
      aa = sqrt(a2)
      cc = sqrt(M_ONE - a2)
      this%rlattice_primitive(1, 1) =  aa
      this%rlattice_primitive(2, 1) =  M_ZERO
      this%rlattice_primitive(3, 1) =  cc
      this%rlattice_primitive(1, 2) = -M_HALF*aa
      this%rlattice_primitive(2, 2) =  M_HALF*sqrt(M_THREE)*aa
      this%rlattice_primitive(3, 2) =  cc
      this%rlattice_primitive(1, 3) = -M_HALF*aa
      this%rlattice_primitive(2, 3) = -M_HALF*sqrt(M_THREE)*aa
      this%rlattice_primitive(3, 3) =  cc
    else
      this%rlattice_primitive(1, 1) = M_ONE
      this%rlattice_primitive(2, 1) = M_ZERO
      this%rlattice_primitive(3, 1) = M_ZERO
      this%rlattice_primitive(1, 2) = cos(M_PI*angles(3)/CNST(180.0))
      this%rlattice_primitive(2, 2) = sin(M_PI*angles(3)/CNST(180.0))
      this%rlattice_primitive(3, 2) = M_ZERO
      this%rlattice_primitive(1, 3) = cos(M_PI*angles(2)/CNST(180.0))
      this%rlattice_primitive(2, 3) = (cos(M_PI*angles(1)/CNST(180.0)) - &
                                      this%rlattice_primitive(1, 2)*this%rlattice_primitive(1, 3))/this%rlattice_primitive(2,2)
      this%rlattice_primitive(3, 3) = sqrt(M_ONE - this%rlattice_primitive(1, 3)**2 - this%rlattice_primitive(2, 3)**2)
    end if

    this%nonorthogonal = any(abs(angles - CNST(90.0)) > M_EPSILON)

    POP_SUB(build_metric_from_angles)
  end subroutine build_metric_from_angles

  !--------------------------------------------------------------
  subroutine reciprocal_lattice(rv, kv, volume, dim, namespace)
    FLOAT,             intent(in)  :: rv(:,:) !< (1:MAX_DIM, 1:MAX_DIM)
    FLOAT,             intent(out) :: kv(:,:) !< (1:MAX_DIM, 1:MAX_DIM)
    FLOAT,             intent(out) :: volume
    integer,           intent(in)  :: dim
    type(namespace_t), optional, intent(in)  :: namespace

    integer :: ii
    FLOAT :: cross(1:3), rv3(1:3, 1:3)

    PUSH_SUB(reciprocal_lattice)

    kv(:,:) = M_ZERO

    select case (dim)
    case (3)
      cross(1:3) = dcross_product(rv(1:3, 2), rv(1:3, 3)) 
      volume = dot_product(rv(1:3, 1), cross(1:3))

      kv(1:3, 1) = dcross_product(rv(:, 2), rv(:, 3))/volume
      kv(1:3, 2) = dcross_product(rv(:, 3), rv(:, 1))/volume
      kv(1:3, 3) = dcross_product(rv(:, 1), rv(:, 2))/volume    
    case (2)
      rv3(1:3, 1:3) = M_ZERO
      rv3(1:2, 1:2) = rv(1:2, 1:2)
      rv3(3, 3) = M_ONE
      cross(1:3) = dcross_product(rv3(1:3, 1), rv3(1:3, 2)) 
      volume = dot_product(rv3(1:3, 3), cross(1:3))

      kv(1:3, 1) = dcross_product(rv3(:, 2), rv3(:, 3))/volume
      kv(1:3, 2) = dcross_product(rv3(:, 3), rv3(:, 1))/volume
    case (1)
      volume = rv(1, 1)
      kv(1, 1) = M_ONE / rv(1, 1)
    case default ! dim > 3
      message(1) = "Reciprocal lattice for dim > 3 assumes no periodicity."
      call messages_warning(1, namespace=namespace)
      volume = M_ONE
      do ii = 1, dim
        kv(ii, ii) = M_ONE/rv(ii,ii)
        !  At least initialize the thing
        volume = volume * sqrt(sum(rv(:, ii)**2))
      end do
    end select

    if (volume < M_ZERO) then
      message(1) = "Your lattice vectors form a left-handed system."
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(reciprocal_lattice)
  end subroutine reciprocal_lattice

end module lattice_vectors_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

