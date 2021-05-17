#include "global.h"

module atom_oct_m
  use global_oct_m
  use messages_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                               &
    atom_init,                            &
    atom_end,                             &
    atom_get_label,                       &
    atom_set_species,                     &
    atom_get_species,                     &
    atom_same_species,                    &
    atom_classical_init,                  &
    atom_classical_end,                   &
    atom_classical_get_label,             &
    atom_classical_write_xyz,             &
    atom_classical_read_xyz

  type, public :: atom_t
    !private
    character(len=LABEL_LEN)  :: label = ""
    type(species_t), pointer  :: species  =>null() !< pointer to species
    integer, dimension(MAX_DIM) :: c   = 0      !< Constrain on te atom (0 or 1)

    !Components of the force
    FLOAT, dimension(MAX_DIM) :: f_ii     = M_ZERO !< Ion-Ion part
    FLOAT, dimension(MAX_DIM) :: f_vdw    = M_ZERO !< Van der Waals part
    FLOAT, dimension(MAX_DIM) :: f_loc    = M_ZERO !< Local electronic part
    FLOAT, dimension(MAX_DIM) :: f_nl     = M_ZERO !< NL electronic part
    FLOAT, dimension(MAX_DIM) :: f_fields = M_ZERO !< Lasers
    FLOAT, dimension(MAX_DIM) :: f_u      = M_ZERO !< Hubbard forces
    FLOAT, dimension(MAX_DIM) :: f_scf    = M_ZERO !< SCF forces
    FLOAT, dimension(MAX_DIM) :: f_nlcc   = M_ZERO !< NLCC forces
  end type atom_t

  type, public :: atom_classical_t
    !private
    character(len=LABEL_LEN)  :: label  = ""
    FLOAT, dimension(MAX_DIM) :: x      = M_ZERO
    FLOAT, dimension(MAX_DIM) :: v      = M_ZERO
    FLOAT, dimension(MAX_DIM) :: f      = M_ZERO
    FLOAT                     :: charge = M_ZERO
  end type atom_classical_t

  interface atom_same_species
    module procedure atom_same_species_aa
    module procedure atom_same_species_as
  end interface atom_same_species

contains

  ! ---------------------------------------------------------
  subroutine atom_init(this, label, species)
    type(atom_t),                      intent(out) :: this
    character(len=*),                  intent(in)  :: label
    type(species_t), target, optional, intent(in)  :: species

    PUSH_SUB(atom_init)

    this%label = trim(adjustl(label))
    this%species  =>null()
    if(present(species))this%species=>species

    this%f_ii      = M_ZERO
    this%f_vdw     = M_ZERO
    this%f_loc     = M_ZERO
    this%f_nl      = M_ZERO
    this%f_fields  = M_ZERO
    this%f_u       = M_ZERO

    POP_SUB(atom_init)
  end subroutine atom_init

  ! ---------------------------------------------------------
  elemental subroutine atom_end(this)
    type(atom_t), intent(inout) :: this

    this%label = ""
    this%species  =>null()

    this%f_ii      = M_ZERO
    this%f_vdw     = M_ZERO
    this%f_loc     = M_ZERO
    this%f_nl      = M_ZERO
    this%f_fields  = M_ZERO    
    this%f_u       = M_ZERO

  end subroutine atom_end

  ! ---------------------------------------------------------

  pure function atom_get_label(this) result(label)
    type(atom_t), intent(in) :: this

    character(len=len_trim(adjustl(this%label))) :: label

    label=trim(adjustl(this%label))

  end function atom_get_label
  
  ! ---------------------------------------------------------
  subroutine atom_set_species(this, species)
    type(atom_t),            intent(inout) :: this
    type(species_t), target, intent(in)    :: species

    PUSH_SUB(atom_set_species)

    this%species=>species
    POP_SUB(atom_set_species)

  end subroutine atom_set_species
  
  ! ---------------------------------------------------------
  subroutine atom_get_species(this, species)
    type(atom_t),    target,  intent(in)  :: this
    type(species_t), pointer, intent(out) :: species

    ! NO PUSH_SUB, called too often

    species => null()
    if(associated(this%species)) species => this%species

  end subroutine atom_get_species
  
  ! ---------------------------------------------------------
  elemental function atom_same_species_aa(this, that) result(is)
    type(atom_t), intent(in) :: this
    type(atom_t), intent(in) :: that

    logical :: is

    is=(atom_get_label(this)==atom_get_label(that))

  end function atom_same_species_aa

  ! ---------------------------------------------------------
  elemental function atom_same_species_as(this, species) result(is)
    type(atom_t),    intent(in) :: this
    type(species_t), intent(in) :: species

    logical :: is

    is=(atom_get_label(this)==species_label(species))

  end function atom_same_species_as

  ! ---------------------------------------------------------
  pure subroutine atom_classical_init(this, label, x, charge)
    type(atom_classical_t), intent(out) :: this
    character(len=*),       intent(in)  :: label
    FLOAT,    dimension(:), intent(in)  :: x
    FLOAT,                  intent(in)  :: charge

    this%label  = trim(adjustl(label))
    this%x      = x
    this%v      = M_ZERO
    this%f      = M_ZERO
    this%charge = charge

  end subroutine atom_classical_init

  ! ---------------------------------------------------------
  elemental subroutine atom_classical_end(this)
    type(atom_classical_t), intent(inout) :: this

    this%label  = ""
    this%x      = M_ZERO
    this%v      = M_ZERO
    this%f      = M_ZERO
    this%charge = M_ZERO

  end subroutine atom_classical_end

  ! ---------------------------------------------------------
  pure function atom_classical_get_label(this) result(label)
    type(atom_classical_t), intent(in) :: this
    !
    character(len=len_trim(adjustl(this%label))) :: label
    !
    label=trim(adjustl(this%label))
    return
  end function atom_classical_get_label
  
  ! ---------------------------------------------------------
  subroutine atom_classical_write_xyz(this, dim, unit)
    type(atom_classical_t), intent(in) :: this
    integer,      optional, intent(in) :: dim
    integer,                intent(in) :: unit

    character(len=27) :: frmt
    integer           :: i, dim_

    PUSH_SUB(atom_classical_write_xyz)
    dim_=MAX_DIM
    if(present(dim))dim_=dim
    write(unit=frmt, fmt="(a10,i2.2,a15)") "(6x,a1,2x,", dim_, "f12.6,a3,f12.6)"
    write(unit=unit, fmt=frmt) this%label(1:1), &
      (units_from_atomic(units_out%length_xyz_file, this%x(i)), i=1, dim_), " # ", this%charge

    POP_SUB(atom_classical_write_xyz)
  end subroutine atom_classical_write_xyz

  ! ---------------------------------------------------------
  subroutine atom_classical_read_xyz(this, dim, unit)
    type(atom_classical_t), intent(inout) :: this
    integer,      optional, intent(in) :: dim
    integer,                intent(in) :: unit

    character(len=27) :: frmt, dum
    integer           :: i, dim_
    FLOAT, dimension(MAX_DIM) :: tmp

    PUSH_SUB(atom_classical_read_xyz)
    dim_=MAX_DIM
    if(present(dim))dim_=dim
    write(unit=frmt, fmt="(a10,i2.2,a15)") "(6x,a1,2x,", dim_, "f12.6,a3,f12.6)"
    read(unit=unit, fmt=frmt) dum, (tmp, i=1, dim_)

    do i = 1, dim_
      this%x(i) = units_to_atomic(units_out%length_xyz_file, tmp(i))
    end do


    POP_SUB(atom_classical_read_xyz)
  end subroutine atom_classical_read_xyz


end module atom_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
