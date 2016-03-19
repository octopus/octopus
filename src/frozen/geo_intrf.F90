#include "global.h"

module geo_intrf_oct_m

  use atom_oct_m
  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use space_oct_m
  use species_oct_m

  implicit none

  private

  public ::      &
    geo_intrf_t

  public ::               &
    geo_intrf_iterator_t

  public ::          &
    geo_intrf_new,   &
    geo_intrf_del,   &
    geo_intrf_assoc, &
    geo_intrf_init,  &
    geo_intrf_next,  &
    geo_intrf_set,   &
    geo_intrf_get,   &
    geo_intrf_copy,  &
    geo_intrf_end

  public ::          &
    GEO_INTRF_OK,    &
    GEO_INTRF_ERROR

  integer, parameter :: GEO_INTRF_OK    = 0
  integer, parameter :: GEO_INTRF_ERROR =-1

  integer, parameter :: GEO_NULL = 0
  integer, parameter :: GEO_ASSC = 1
  integer, parameter :: GEO_ALLC = 2

  type :: geo_intrf_t
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(geometry_t),    pointer :: pgeo   =>null()
    integer                      :: type   = GEO_NULL
  end type geo_intrf_t

  type :: geo_intrf_iterator_t
    private
    type(geo_intrf_t), pointer :: self =>null()
    type(geometry_t),  pointer :: pgeo =>null()
    integer                    :: natm = 0
    integer                    :: nspc = 0
  end type geo_intrf_iterator_t

  interface geo_intrf_new
    module procedure geo_intrf_new_geometry
    module procedure geo_intrf_new_init
  end interface geo_intrf_new

  interface geo_intrf_init
    module procedure geo_intrf_init_type
    module procedure geo_intrf_init_copy
    module procedure geo_intrf_iterator_init_type
    module procedure geo_intrf_iterator_init_iterator
  end interface geo_intrf_init

  interface geo_intrf_next
    module procedure geo_intrf_iterator_next_atom
    module procedure geo_intrf_iterator_next_species
  end interface geo_intrf_next

  interface geo_intrf_set
    module procedure geo_intrf_set_geometry
  end interface geo_intrf_set

  interface geo_intrf_get
    module procedure geo_intrf_get_config
    module procedure geo_intrf_get_space
    module procedure geo_intrf_get_geometry
  end interface geo_intrf_get

  interface geo_intrf_copy
    module procedure geo_intrf_copy_type
    module procedure geo_intrf_iterator_copy
  end interface geo_intrf_copy

  interface geo_intrf_end
    module procedure geo_intrf_end_type
    module procedure geo_intrf_iterator_end
  end interface geo_intrf_end

contains

  ! ---------------------------------------------------------
  subroutine geo_intrf_new_geometry(this, that)
    type(geo_intrf_t), intent(inout) :: this
    type(geometry_t), pointer        :: that

    PUSH_SUB(geo_intrf_new_geometry)
    
    nullify(that)
    SAFE_ALLOCATE(that)
    call geometry_nullify(that)
    call geo_intrf_set(this, that)
    this%type = GEO_ALLC

    POP_SUB(geo_intrf_new_geometry)
  end subroutine geo_intrf_new_geometry

  ! ---------------------------------------------------------
  subroutine geo_intrf_new_init(this, that, geo_init)
    type(geo_intrf_t), intent(inout) :: this
    type(geometry_t), pointer        :: that

    interface
      subroutine geo_init(this, space, config)
        use geometry_oct_m
        use json_oct_m
        use space_oct_m
        type(geometry_t),    intent(out) :: this
        type(space_t),       intent(in)  :: space
        type(json_object_t), intent(in)  :: config
      end subroutine geo_init
    end interface

    PUSH_SUB(geo_intrf_new_init)
    
    nullify(that)
    call geo_intrf_new(this, that)
    call geo_init(that, this%space, this%config)

    POP_SUB(geo_intrf_new_init)
  end subroutine geo_intrf_new_init

  ! ---------------------------------------------------------
  subroutine geo_intrf_del(this)
    type(geo_intrf_t), intent(inout) :: this

    PUSH_SUB(geo_intrf_del)

    ASSERT(this%type==GEO_ALLC)
    ASSERT(associated(this%pgeo))
    call geometry_end(this%pgeo)
    SAFE_DEALLOCATE_P(this%pgeo)
    nullify(this%pgeo)
    this%type = GEO_NULL
    
    POP_SUB(geo_intrf_del)
  end subroutine geo_intrf_del

  ! ---------------------------------------------------------
  function geo_intrf_assoc(this) result(that)
    type(geo_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(geo_intrf_assoc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    select case(this%type)
    case(GEO_NULL)
      ASSERT(.not.associated(this%pgeo))
      that = .false.
    case(GEO_ASSC,GEO_ALLC)
      ASSERT(associated(this%pgeo))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(geo_intrf_assoc)
  end function geo_intrf_assoc

  ! ---------------------------------------------------------
  subroutine geo_intrf_init_type(this, space, config)
    type(geo_intrf_t),           intent(out) :: this
    type(space_t),       target, intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(geo_intrf_init_type)
    
    this%config => config
    this%space => space
    nullify(this%pgeo)
    this%type = GEO_NULL

    POP_SUB(geo_intrf_init_type)
  end subroutine geo_intrf_init_type

  ! ---------------------------------------------------------
  subroutine geo_intrf_init_copy(this, that)
    type(geo_intrf_t), intent(out) :: this
    type(geo_intrf_t), intent(in)  :: that

    PUSH_SUB(geo_intrf_init_copy)
    
    ASSERT(associated(that%config))
    ASSERT(associated(that%space))
    call geo_intrf_init(this, that%space, that%config)

    POP_SUB(geo_intrf_init_copy)
  end subroutine geo_intrf_init_copy

  ! ---------------------------------------------------------
  subroutine geo_intrf_set_geometry(this, that)
    type(geo_intrf_t),        intent(inout) :: this
    type(geometry_t), target, intent(in)    :: that

    PUSH_SUB(geo_intrf_set_geometry)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(.not.associated(this%pgeo))
    ASSERT(this%type==GEO_NULL)
    this%pgeo => that
    this%type = GEO_ASSC

    POP_SUB(geo_intrf_set_geometry)
  end subroutine geo_intrf_set_geometry

  ! ---------------------------------------------------------
  subroutine geo_intrf_get_config(this, that)
    type(geo_intrf_t),    intent(in) :: this
    type(json_object_t), pointer     :: that

    PUSH_SUB(geo_intrf_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(geo_intrf_get_config)
  end subroutine geo_intrf_get_config

  ! ---------------------------------------------------------
  subroutine geo_intrf_get_space(this, that)
    type(geo_intrf_t), intent(in) :: this
    type(space_t),    pointer     :: that

    PUSH_SUB(geo_intrf_get_space)

    nullify(that)
    if(associated(this%space)) that => this%space

    POP_SUB(geo_intrf_get_space)
  end subroutine geo_intrf_get_space

  ! ---------------------------------------------------------
  subroutine geo_intrf_get_geometry(this, that)
    type(geo_intrf_t), intent(in) :: this
    type(geometry_t), pointer     :: that

    PUSH_SUB(geo_intrf_get_geometry)

    nullify(that)
    if(associated(this%pgeo)) that => this%pgeo

    POP_SUB(geo_intrf_get_geometry)
  end subroutine geo_intrf_get_geometry

  ! ---------------------------------------------------------
  subroutine geo_intrf_copy_type(this, that)
    type(geo_intrf_t), intent(inout) :: this
    type(geo_intrf_t), intent(in)    :: that

    type(geometry_t), pointer :: pgeo

    PUSH_SUB(geo_intrf_copy_type)

    nullify(pgeo)
    call geo_intrf_end(this)
    this%config => that%config
    this%space => that%space
    select case(that%type)
    case(GEO_NULL)
      nullify(this%pgeo)
    case(GEO_ASSC)
      this%pgeo => that%pgeo
    case(GEO_ALLC)
      call geo_intrf_new(this, pgeo)
      call geometry_copy(pgeo, that%pgeo)
      nullify(pgeo)
    case default
      ASSERT(.false.)
    end select
    this%type = that%type

    POP_SUB(geo_intrf_copy_type)
  end subroutine geo_intrf_copy_type

  ! ---------------------------------------------------------
  subroutine geo_intrf_end_type(this)
    type(geo_intrf_t), intent(inout) :: this

    PUSH_SUB(geo_intrf_end_type)

    if(this%type==GEO_ALLC) call geo_intrf_del(this)
    nullify(this%config, this%space, this%pgeo)
    this%type = GEO_NULL

    POP_SUB(geo_intrf_end_type)
  end subroutine geo_intrf_end_type

  ! ---------------------------------------------------------
  subroutine geo_intrf_iterator_init_type(this, that)
    type(geo_intrf_iterator_t), intent(out) :: this
    type(geo_intrf_t),  target, intent(in)  :: that

    PUSH_SUB(geo_intrf_iterator_init_type)

    this%self => that
    this%pgeo => that%pgeo
    this%natm = 0
    this%nspc = 0

    POP_SUB(geo_intrf_iterator_init_type)
  end subroutine geo_intrf_iterator_init_type

  ! ---------------------------------------------------------
  subroutine geo_intrf_iterator_init_iterator(this, that)
    type(geo_intrf_iterator_t), intent(out) :: this
    type(geo_intrf_iterator_t), intent(in)  :: that

    PUSH_SUB(geo_intrf_iterator_init_iterator)

    ASSERT(associated(that%self))
    ASSERT(associated(that%pgeo))
    call geo_intrf_iterator_copy(this, that)

    POP_SUB(geo_intrf_iterator_init_iterator)
  end subroutine geo_intrf_iterator_init_iterator

  ! ---------------------------------------------------------
  subroutine geo_intrf_iterator_next_atom(this, atom, ierr)
    type(geo_intrf_iterator_t), intent(inout) :: this
    type(atom_t),              pointer        :: atom
    integer,          optional, intent(out)   :: ierr

    integer :: jerr

    PUSH_SUB(geo_intrf_iterator_next_atom)

    nullify(atom)
    jerr = GEO_INTRF_ERROR
    if(associated(this%self).and.associated(this%pgeo))then
      this%natm = this%natm + 1
      if(.not.(this%natm>this%pgeo%natoms))then
        atom => this%pgeo%atom(this%natm)
        jerr = GEO_INTRF_OK
      end if
    end if
    if(present(ierr)) ierr = jerr

    POP_SUB(geo_intrf_iterator_next_atom)
  end subroutine geo_intrf_iterator_next_atom

  ! ---------------------------------------------------------
  subroutine geo_intrf_iterator_next_species(this, spec, ierr)
    type(geo_intrf_iterator_t), intent(inout) :: this
    type(species_t),           pointer        :: spec
    integer,          optional, intent(out)   :: ierr

    integer :: jerr

    PUSH_SUB(geo_intrf_iterator_next_species)

    nullify(spec)
    jerr = GEO_INTRF_ERROR
    if(associated(this%self).and.associated(this%pgeo))then
      this%nspc = this%nspc + 1
      if(.not.(this%nspc>this%pgeo%nspecies))then
        spec => this%pgeo%species(this%nspc)
        jerr = GEO_INTRF_OK
      end if
    end if
    if(present(ierr)) ierr = jerr

    POP_SUB(geo_intrf_iterator_next_species)
  end subroutine geo_intrf_iterator_next_species

  ! ---------------------------------------------------------
  subroutine geo_intrf_iterator_copy(this, that)
    type(geo_intrf_iterator_t),         intent(out) :: this
    type(geo_intrf_iterator_t), target, intent(in)  :: that

    PUSH_SUB(geo_intrf_iterator_copy)

    this%self => that%self
    this%pgeo => that%pgeo
    this%natm = that%natm
    this%nspc = that%nspc

    POP_SUB(geo_intrf_iterator_copy)
  end subroutine geo_intrf_iterator_copy

  ! ---------------------------------------------------------
  subroutine geo_intrf_iterator_end(this)
    type(geo_intrf_iterator_t), intent(inout) :: this

    PUSH_SUB(geo_intrf_iterator_end)

    nullify(this%self, this%pgeo)
    this%natm = 0
    this%nspc = 0

    POP_SUB(geo_intrf_iterator_end)
  end subroutine geo_intrf_iterator_end

end module geo_intrf_oct_m

!! Local Variables:
!! mode: f90
!! End:

