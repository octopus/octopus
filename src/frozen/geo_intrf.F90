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

  public ::          &
    GEO_INTRF_OK,    &
    GEO_INTRF_ERROR

  public ::      &
    geo_intrf_t

  public ::               &
    geo_intrf_iterator_t

  public ::          &
    geo_intrf_new,   &
    geo_intrf_del,   &
    geo_intrf_assoc, &
    geo_intrf_alloc, &
    geo_intrf_init,  &
    geo_intrf_next,  &
    geo_intrf_set,   &
    geo_intrf_get,   &
    geo_intrf_copy,  &
    geo_intrf_end

  integer, parameter :: GEO_INTRF_OK    = 0
  integer, parameter :: GEO_INTRF_ERROR =-1

  integer, parameter :: GEO_INTRF_DISA = 0
  integer, parameter :: GEO_INTRF_NULL = 1
  integer, parameter :: GEO_INTRF_ASSC = 2
  integer, parameter :: GEO_INTRF_ALLC = 3

  type :: geo_intrf_t
    private
    type(geometry_t), pointer :: self =>null()
    integer                   :: type = GEO_INTRF_DISA
  end type geo_intrf_t

  type :: geo_intrf_iterator_t
    private
    type(geo_intrf_t), pointer :: self =>null()
    type(geometry_t),  pointer :: pgeo =>null()
    integer                    :: natm = 0
    integer                    :: nspc = 0
  end type geo_intrf_iterator_t

  interface geo_intrf_new
    module procedure geo_intrf_new_type
    module procedure geo_intrf_new_pass
    module procedure geo_intrf_new_copy
  end interface geo_intrf_new

  interface geo_intrf_del
    module procedure geo_intrf_del_type
    module procedure geo_intrf_del_pass
  end interface geo_intrf_del

  interface geo_intrf_init
    module procedure geo_intrf_init_type
    module procedure geo_intrf_iterator_init_type
    module procedure geo_intrf_iterator_init_iterator
  end interface geo_intrf_init

  interface geo_intrf_next
    module procedure geo_intrf_iterator_next_atom
    module procedure geo_intrf_iterator_next_species
  end interface geo_intrf_next

  interface geo_intrf_copy
    module procedure geo_intrf_copy_type
    module procedure geo_intrf_iterator_copy
  end interface geo_intrf_copy

  interface geo_intrf_end
    module procedure geo_intrf_end_type
    module procedure geo_intrf_end_pass
    module procedure geo_intrf_iterator_end
  end interface geo_intrf_end

contains

  ! ---------------------------------------------------------
  subroutine geo_intrf__new__(this, that)
    type(geo_intrf_t), intent(inout) :: this
    type(geometry_t), pointer        :: that

    PUSH_SUB(geo_intrf__new__)

    ASSERT(geo_intrf_isnull(this))
    nullify(that)
    SAFE_ALLOCATE(that)
    call geo_intrf_set(this, that)
    this%type = GEO_INTRF_ALLC

    POP_SUB(geo_intrf__new__)
  end subroutine geo_intrf__new__

  ! ---------------------------------------------------------
  subroutine geo_intrf__del__(this)
    type(geo_intrf_t), intent(inout) :: this

    PUSH_SUB(geo_intrf__del__)

    ASSERT(geo_intrf_assoc(this))
    if(geo_intrf_alloc(this))then
      SAFE_DEALLOCATE_P(this%self)
    end if
    nullify(this%self)
    this%type = GEO_INTRF_NULL

    POP_SUB(geo_intrf__del__)
  end subroutine geo_intrf__del__

  ! ---------------------------------------------------------
  subroutine geo_intrf_new_type(this, space, config)
    type(geo_intrf_t),   intent(inout) :: this
    type(space_t),       intent(in)    :: space
    type(json_object_t), intent(in)    :: config

    PUSH_SUB(geo_intrf_new_type)

    ASSERT(geo_intrf_isnull(this))
    call geo_intrf_new(this, init)

    POP_SUB(geo_intrf_new_type)

  contains

    subroutine init(this)
      type(geometry_t), intent(out) :: this

      PUSH_SUB(geo_intrf_new_type.init)

      call geometry_init_from_data_object(this, space, config)

      POP_SUB(geo_intrf_new_type.init)
    end subroutine init

  end subroutine geo_intrf_new_type

  ! ---------------------------------------------------------
  subroutine geo_intrf_new_pass(this, init)
    type(geo_intrf_t), intent(inout) :: this
    
    interface
      subroutine init(this)
        use geometry_oct_m
        type(geometry_t), intent(out) :: this
      end subroutine init
    end interface

    type(geometry_t), pointer :: self

    PUSH_SUB(geo_intrf_new_pass)
    
    ASSERT(geo_intrf_isnull(this))
    nullify(self)
    call geo_intrf__new__(this, self)
    call init(self)
    nullify(self)

    POP_SUB(geo_intrf_new_pass)
  end subroutine geo_intrf_new_pass

  ! ---------------------------------------------------------
  subroutine geo_intrf_new_copy(this, that)
    type(geo_intrf_t), intent(inout) :: this
    type(geo_intrf_t), intent(in)    :: that

    type(geometry_t), pointer :: self

    PUSH_SUB(geo_intrf_new_copy)

    ASSERT(geo_intrf_isnull(this))
    ASSERT(geo_intrf_assoc(that))
    nullify(self)
    call geo_intrf__new__(this, self)
    call geometry_copy(self, that%self)
    nullify(self)

    POP_SUB(geo_intrf_new_copy)
  end subroutine geo_intrf_new_copy

  ! ---------------------------------------------------------
  subroutine geo_intrf_del_type(this)
    type(geo_intrf_t), intent(inout) :: this

    PUSH_SUB(geo_intrf_del_type)

    ASSERT(geo_intrf_assoc(this))
    call geo_intrf_del(this, geometry_end)

    POP_SUB(geo_intrf_del_type)
  end subroutine geo_intrf_del_type

  ! ---------------------------------------------------------
  subroutine geo_intrf_del_pass(this, finis)
    type(geo_intrf_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use geometry_oct_m
        type(geometry_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(geo_intrf_del_pass)

    ASSERT(geo_intrf_assoc(this))
    if(geo_intrf_alloc(this)) call finis(this%self)
    call geo_intrf__del__(this)

    POP_SUB(geo_intrf_del_pass)
  end subroutine geo_intrf_del_pass

  ! ---------------------------------------------------------
  function geo_intrf_isnull(this) result(that)
    type(geo_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(geo_intrf_isnull)

    select case(this%type)
    case(GEO_INTRF_DISA)
      that = .false.
    case(GEO_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .true.
    case(GEO_INTRF_ASSC,GEO_INTRF_ALLC)
      ASSERT(associated(this%self))
      that = .false.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(geo_intrf_isnull)
  end function geo_intrf_isnull

  ! ---------------------------------------------------------
  function geo_intrf_assoc(this) result(that)
    type(geo_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(geo_intrf_assoc)

    select case(this%type)
    case(GEO_INTRF_DISA)
      that = .false.
    case(GEO_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .false.
    case(GEO_INTRF_ASSC, GEO_INTRF_ALLC)
      ASSERT(associated(this%self))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(geo_intrf_assoc)
  end function geo_intrf_assoc

  ! ---------------------------------------------------------
  function geo_intrf_alloc(this) result(that)
    type(geo_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(geo_intrf_alloc)

    select case(this%type)
    case(GEO_INTRF_DISA)
      that = .false.
    case(GEO_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .false.
    case(GEO_INTRF_ASSC)
      ASSERT(associated(this%self))
      that = .false.
    case(GEO_INTRF_ALLC)
      ASSERT(associated(this%self))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(geo_intrf_alloc)
  end function geo_intrf_alloc

  ! ---------------------------------------------------------
  subroutine geo_intrf_init_type(this)
    type(geo_intrf_t), intent(out) :: this

    PUSH_SUB(geo_intrf_init_type)
    
    nullify(this%self)
    this%type = GEO_INTRF_NULL

    POP_SUB(geo_intrf_init_type)
  end subroutine geo_intrf_init_type

  ! ---------------------------------------------------------
  subroutine geo_intrf_set(this, that)
    type(geo_intrf_t),        intent(inout) :: this
    type(geometry_t), target, intent(in)    :: that

    PUSH_SUB(geo_intrf_set)

    ASSERT(geo_intrf_isnull(this))
    this%self => that
    this%type = GEO_INTRF_ASSC

    POP_SUB(geo_intrf_set)
  end subroutine geo_intrf_set

  ! ---------------------------------------------------------
  subroutine geo_intrf_get(this, that)
    type(geo_intrf_t), intent(in) :: this
    type(geometry_t), pointer     :: that

    PUSH_SUB(geo_intrf_get)

    ASSERT(this%type/=GEO_INTRF_DISA)
    nullify(that)
    if(geo_intrf_assoc(this)) that => this%self

    POP_SUB(geo_intrf_get)
  end subroutine geo_intrf_get

  ! ---------------------------------------------------------
  subroutine geo_intrf_copy_type(this, that)
    type(geo_intrf_t), intent(inout) :: this
    type(geo_intrf_t), intent(in)    :: that

    PUSH_SUB(geo_intrf_copy_type)

    call geo_intrf_end(this)
    select case(that%type)
    case(GEO_INTRF_DISA)
    case(GEO_INTRF_NULL)
    case(GEO_INTRF_ASSC)
      call geo_intrf_set(this, that%self)
    case(GEO_INTRF_ALLC)
      call geo_intrf_new(this, that)
    case default
      ASSERT(.false.)
    end select
    this%type = that%type

    POP_SUB(geo_intrf_copy_type)
  end subroutine geo_intrf_copy_type

  ! ---------------------------------------------------------
  subroutine geo_intrf__end__(this)
    type(geo_intrf_t), intent(inout) :: this

    PUSH_SUB(geo_intrf__end__)

    nullify(this%self)
    this%type = GEO_INTRF_DISA

    POP_SUB(geo_intrf__end__)
  end subroutine geo_intrf__end__

  ! ---------------------------------------------------------
  subroutine geo_intrf_end_type(this)
    type(geo_intrf_t), intent(inout) :: this

    PUSH_SUB(geo_intrf_end_type)

    call geo_intrf_end(this, geometry_end)

    POP_SUB(geo_intrf_end_type)
  end subroutine geo_intrf_end_type

  ! ---------------------------------------------------------
  subroutine geo_intrf_end_pass(this, finis)
    type(geo_intrf_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use geometry_oct_m
        type(geometry_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(geo_intrf_end_pass)

    if(geo_intrf_assoc(this)) call geo_intrf_del(this, finis)
    call geo_intrf__end__(this)

    POP_SUB(geo_intrf_end_pass)
  end subroutine geo_intrf_end_pass

  ! ---------------------------------------------------------
  subroutine geo_intrf_iterator_init_type(this, that)
    type(geo_intrf_iterator_t), intent(out) :: this
    type(geo_intrf_t),  target, intent(in)  :: that

    PUSH_SUB(geo_intrf_iterator_init_type)

    this%self => that
    call geo_intrf_get(that, this%pgeo)
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

