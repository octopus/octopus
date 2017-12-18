#include "global.h"

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

#define BASE_EXTENDED_ITERATOR_TYPE

module base_geometry_oct_m

  use atom_oct_m
  use geo_build_oct_m
  use geo_intrf_oct_m
  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use space_oct_m
  use species_oct_m

#define BASE_TEMPLATE_NAME base_geometry
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private

  public ::          &
    base_geometry_t

  public ::                  &
    base_geometry__init__,   &
    base_geometry__acc__,    &
    base_geometry__update__, &
    base_geometry__reset__,  &
    base_geometry__copy__,   &
    base_geometry__end__

  public ::               &
    base_geometry_new,    &
    base_geometry_del,    &
    base_geometry_init,   &
    base_geometry_acc,    &
    base_geometry_update, &
    base_geometry_reset,  &
    base_geometry_set,    &
    base_geometry_get,    &
    base_geometry_copy,   &
    base_geometry_end

#define BASE_TEMPLATE_NAME base_geometry
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  type :: base_geometry_t
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(refcount_t),    pointer :: rcnt   =>null()
    type(geo_intrf_t)            :: igeo
    type(base_geometry_dict_t)   :: dict
  end type base_geometry_t

  type :: base_geometry_iterator_t
    private
    type(base_geometry_t),      pointer :: self =>null()
    type(geo_intrf_iterator_t)          :: gitr
    type(base_geometry_dict_iterator_t) :: iter
  end type base_geometry_iterator_t

  interface base_geometry__init__
    module procedure base_geometry__init__type
    module procedure base_geometry__init__pass
    module procedure base_geometry__init__copy
  end interface base_geometry__init__

  interface base_geometry_new
    module procedure base_geometry_new_type
    module procedure base_geometry_new_pass
  end interface base_geometry_new

  interface base_geometry_init
    module procedure base_geometry_init_type
  end interface base_geometry_init

  interface base_geometry_set
    module procedure base_geometry_set_geometry
  end interface base_geometry_set

  interface base_geometry_get
    module procedure base_geometry_get_config
    module procedure base_geometry_get_space
    module procedure base_geometry_get_geometry
    module procedure base_geometry_get_sub_geometry
  end interface base_geometry_get

  interface base_geometry__sets__
    module procedure base_geometry__sets__info
    module procedure base_geometry__sets__type
  end interface base_geometry__sets__

  interface base_geometry_next
    module procedure base_geometry_iterator_next_atom
    module procedure base_geometry_iterator_next_species
    module procedure base_geometry_iterator_next_geometry
  end interface base_geometry_next

contains

#define BASE_TEMPLATE_NAME base_geometry
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  function base_geometry_new_type(space, config) result(this)
    type(space_t),       intent(in) :: space
    type(json_object_t), intent(in) :: config

    type(base_geometry_t), pointer :: this

    PUSH_SUB(base_geometry_new_type)

    this => base_geometry_new(space, config, base_geometry_init_type)

    POP_SUB(base_geometry_new_type)
  end function base_geometry_new_type
  
  ! ---------------------------------------------------------
  function base_geometry_new_pass(space, config, init) result(this)
    type(space_t),       intent(in) :: space
    type(json_object_t), intent(in) :: config

    type(base_geometry_t), pointer :: this

    interface
      subroutine init(this, space, config)
        use json_oct_m
        use space_oct_m
        import :: base_geometry_t
        type(base_geometry_t), intent(out) :: this
        type(space_t),         intent(in)  :: space
        type(json_object_t),   intent(in)  :: config
      end subroutine init
    end interface

    PUSH_SUB(base_geometry_new_pass)

    nullify(this)
    SAFE_ALLOCATE(this)
    call init(this, space, config)
    ASSERT(associated(this%rcnt))
    call refcount_set(this%rcnt, dynamic=.true.)
    
    POP_SUB(base_geometry_new_pass)
  end function base_geometry_new_pass

  ! ---------------------------------------------------------
  subroutine base_geometry__init__type(this, space, config)
    type(base_geometry_t),       intent(out) :: this
    type(space_t),       target, intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t),          pointer :: cnfg
    character(len=BASE_GEOMETRY_NAME_LEN) :: name
    integer                               :: ierr
    logical                               :: dflt, refs, accu, allc

    PUSH_SUB(base_geometry__init__type)

    nullify(cnfg)
    this%config => config
    this%space => space
    this%rcnt => refcount_new()
    call geo_intrf_init(this%igeo)
    call json_get(this%config, "default", dflt, ierr)
    if(ierr/=JSON_OK) dflt = .true.
    if(dflt)then
      refs = .false.
      call json_get(this%config, "reference", name, ierr)
      if(ierr==JSON_OK) refs = .true.
      call json_get(this%config, "reduce", accu, ierr)
      if(ierr/=JSON_OK) accu = .false.
      ASSERT(.not.(refs.and.accu))
      call json_get(this%config, "molecule", cnfg, ierr)
      if(ierr==JSON_OK) call geo_intrf_new(this%igeo, this%space, cnfg)
      nullify(cnfg)
      allc = geo_intrf_assoc(this%igeo)
      ASSERT((.not.allc.and.(refs.or.accu)).or.(.not.(refs.or.accu).and.allc))
    end if
    call base_geometry_dict_init(this%dict)

    POP_SUB(base_geometry__init__type)
  end subroutine base_geometry__init__type

  ! ---------------------------------------------------------
  subroutine base_geometry__init__pass(this, init)
    type(base_geometry_t), intent(inout) :: this

    interface
      subroutine init(this)
        use geometry_oct_m
        type(geometry_t), intent(out) :: this
      end subroutine init
    end interface

    integer :: ierr
    logical :: dflt

    PUSH_SUB(base_geometry__init__pass)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(.not.geo_intrf_assoc(this%igeo))
    call json_get(this%config, "default", dflt, ierr)
    if(ierr/=JSON_OK) dflt = .true.
    ASSERT(.not.dflt)
    call geo_intrf_new(this%igeo, init)
    ASSERT(geo_intrf_alloc(this%igeo))

    POP_SUB(base_geometry__init__pass)
  end subroutine base_geometry__init__pass

  ! ---------------------------------------------------------
  subroutine base_geometry__init__copy(this, that)
    type(base_geometry_t), intent(out) :: this
    type(base_geometry_t), intent(in)  :: that

    PUSH_SUB(base_geometry__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%space))
    call base_geometry__init__(this, that%space, that%config)
    call geo_intrf_copy(this%igeo, that%igeo)

    POP_SUB(base_geometry__init__copy)
  end subroutine base_geometry__init__copy

  ! ---------------------------------------------------------
  subroutine base_geometry_init_type(this, space, config)
    type(base_geometry_t), intent(out) :: this
    type(space_t),         intent(in)  :: space
    type(json_object_t),   intent(in)  :: config

    PUSH_SUB(base_geometry_init_type)

    call base_geometry__init__(this, space, config)

    POP_SUB(base_geometry_init_type)
  end subroutine base_geometry_init_type

  ! ---------------------------------------------------------
  subroutine base_geometry__acc__(this, that, config)
    type(base_geometry_t),         intent(inout) :: this
    type(base_geometry_t),         intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    type(geometry_t), pointer :: pgeo
    type(geo_build_t)         :: bgeo

    PUSH_SUB(base_geometry__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(that%config))
    ASSERT(associated(that%space))
    ASSERT(geo_intrf_assoc(that%igeo))
    ASSERT(this%space==that%space)
    nullify(pgeo)
    call geo_build_init(bgeo, this%space)
    if(geo_intrf_assoc(this%igeo))then
      ASSERT(geo_intrf_alloc(this%igeo))
      call base_geometry_get(this, pgeo)
      ASSERT(associated(pgeo))
      call geo_build_extend(bgeo, pgeo)
      nullify(pgeo)
    end if
    call base_geometry__reset__(this)
    call base_geometry_get(that, pgeo)
    ASSERT(associated(pgeo))
    if(present(config))then
      call geo_build_extend(bgeo, pgeo, config)
    else
      call geo_build_extend(bgeo, pgeo)
    end if
    nullify(pgeo)
    call geo_intrf_new(this%igeo, init)
    call geo_build_end(bgeo)
 
    POP_SUB(base_geometry__acc__)
    
  contains

    subroutine init(this)
      type(geometry_t),    intent(out) :: this
      
      PUSH_SUB(base_geometry__acc__.init)
      
      call geo_build_export(bgeo, this)

      POP_SUB(base_geometry__acc__.init)
    end subroutine init

  end subroutine base_geometry__acc__

  ! ---------------------------------------------------------
  subroutine base_geometry__update__(this)
    type(base_geometry_t), intent(inout) :: this

    PUSH_SUB(base_geometry__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))

    POP_SUB(base_geometry__update__)
  end subroutine base_geometry__update__

  ! ---------------------------------------------------------
  subroutine base_geometry__reset__(this)
    type(base_geometry_t), intent(inout) :: this

    PUSH_SUB(base_geometry__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    if(geo_intrf_assoc(this%igeo)) call geo_intrf_del(this%igeo)

    POP_SUB(base_geometry__reset__)
  end subroutine base_geometry__reset__

  ! ---------------------------------------------------------
  subroutine base_geometry_acc(this)
    type(base_geometry_t), intent(inout) :: this

    type(geo_build_t) :: bgeo

    PUSH_SUB(base_geometry_acc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(base_geometry_dict_len(this%dict)>0)
    call geo_build_init(bgeo, this%space)
    call base_geometry__reset__(this)
    call base_geometry__reduce__(this, acc)
    call base_geometry__update__(this)
    call geo_intrf_new(this%igeo, init)
    call geo_build_end(bgeo)

    POP_SUB(base_geometry_acc)

  contains

    subroutine init(this)
      type(geometry_t), intent(out) :: this

      PUSH_SUB(base_geometry_acc.init)

      call geo_build_export(bgeo, this)

      POP_SUB(base_geometry_acc.init)
    end subroutine init

    subroutine acc(this, that, config)
      type(base_geometry_t),         intent(inout) :: this
      type(base_geometry_t),         intent(in)    :: that
      type(json_object_t), optional, intent(in)    :: config

      type(json_array_t), pointer :: list
      type(geometry_t),   pointer :: pgeo
      integer                     :: ierr

      PUSH_SUB(base_geometry_acc.acc)

      ASSERT(associated(that%config))
      ASSERT(associated(that%space))
      ASSERT(geo_intrf_assoc(that%igeo))
      ASSERT(this%space==that%space)
      nullify(list, pgeo)
      call base_geometry_get(that, pgeo)
      ASSERT(associated(pgeo))
      if(present(config))then
        call json_get(config, "positions", list, ierr)
        if(ierr/=JSON_OK) nullify(list)
      end if
      if(associated(list))then
        ASSERT(json_len(list)>0)
        call geo_build_extend(bgeo, pgeo, list)
      else
        call geo_build_extend(bgeo, pgeo)
      end if
      nullify(list, pgeo)

      POP_SUB(base_geometry_acc.acc)
    end subroutine acc

  end subroutine base_geometry_acc

  ! ---------------------------------------------------------
  subroutine base_geometry_update(this)
    type(base_geometry_t), intent(inout) :: this

    PUSH_SUB(base_geometry_update)

    call base_geometry__apply__(this, base_geometry__update__)

    POP_SUB(base_geometry_update)
  end subroutine base_geometry_update

  ! ---------------------------------------------------------
  subroutine base_geometry_reset(this)
    type(base_geometry_t), intent(inout) :: this

    PUSH_SUB(base_geometry_reset)

    call base_geometry__apply__(this, base_geometry__reset__)

    POP_SUB(base_geometry_reset)
  end subroutine base_geometry_reset

  ! ---------------------------------------------------------
  subroutine base_geometry__sets__info(this, name, lock, active)
    type(base_geometry_t), intent(inout) :: this
    character(len=*),      intent(in)    :: name
    logical,     optional, intent(in)    :: lock
    logical,     optional, intent(in)    :: active

    PUSH_SUB(base_geometry__sets__info)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(len_trim(adjustl(name))>0)
    if(present(lock)) continue
    if(present(active)) continue

    POP_SUB(base_geometry__sets__info)
  end subroutine base_geometry__sets__info
  
  ! ---------------------------------------------------------
  subroutine base_geometry__sets__type(this, name, that, config, lock, active)
    type(base_geometry_t), intent(inout) :: this
    character(len=*),      intent(in)    :: name
    type(base_geometry_t), intent(in)    :: that
    type(json_object_t),   intent(in)    :: config
    logical,     optional, intent(in)    :: lock
    logical,     optional, intent(in)    :: active

    type(json_array_iterator_t)           :: iter
    character(len=BASE_GEOMETRY_NAME_LEN) :: rnam
    type(json_object_t),          pointer :: cnfg
    type(json_array_t),           pointer :: list
    type(geometry_t),             pointer :: pgeo
    integer                               :: indx, ierr
    logical                               :: dflt, actv, accu

    PUSH_SUB(base_geometry__sets__type)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(len_trim(adjustl(name))>0)
    ASSERT(associated(that%config))
    ASSERT(associated(that%space))
    ASSERT(this%space==that%space)
    nullify(cnfg, list, pgeo)
    if(present(lock)) continue
    call json_get(this%config, "default", dflt, ierr)
    if(ierr==JSON_OK) dflt = .true.
    if(dflt)then
      actv = .true.
      if(present(active)) actv = active
      call json_get(this%config, "reduce", accu, ierr)
      if(ierr/=JSON_OK) accu = .false.
      if(actv.and.accu)then
        call json_get(config, "positions", list, ierr)
        if(ierr/=JSON_OK) nullify(list)
        if(associated(list))then
          ASSERT(json_len(list)>0)
          call json_init(iter, list)
          do
            nullify(cnfg)
            call json_next(iter, cnfg, ierr)
            if(ierr/=JSON_OK)exit
            call base_geometry__acc__(this, that, cnfg)
          end do
          call json_end(iter)
          nullify(cnfg, list)
        else
          call base_geometry__acc__(this, that)
        end if
        ASSERT(geo_intrf_alloc(this%igeo))
      end if
      call json_get(this%config, "reference", rnam, ierr)
      if(ierr==JSON_OK)then
        ASSERT(.not.accu)
        indx = index(trim(adjustl(rnam)), trim(adjustl(name)))
        if(indx==1)then
          ASSERT(.not.geo_intrf_assoc(this%igeo))
          call base_geometry_get(this, trim(adjustl(rnam)), pgeo)
          ASSERT(associated(pgeo))
          call geo_intrf_set(this%igeo, pgeo)
          nullify(pgeo)
          ASSERT(geo_intrf_assoc(this%igeo))
        end if
      end if
    end if

    POP_SUB(base_geometry__sets__type)
  end subroutine base_geometry__sets__type
  
  ! ---------------------------------------------------------
  subroutine base_geometry__dels__(this, name, that)
    type(base_geometry_t), intent(inout) :: this
    character(len=*),      intent(in)    :: name
    type(base_geometry_t), intent(in)    :: that

    PUSH_SUB(base_geometry__dels__)

    ASSERT(associated(this%config))
    ASSERT(len_trim(adjustl(name))>0)
    ASSERT(associated(that%config))
    ASSERT(associated(that%space))
    ASSERT(this%space==that%space)

    POP_SUB(base_geometry__dels__)
  end subroutine base_geometry__dels__

  ! ---------------------------------------------------------
  subroutine base_geometry_get_sub_geometry(this, name, that)
    type(base_geometry_t),     intent(in)  :: this
    character(len=*),          intent(in)  :: name
    type(geometry_t), pointer, intent(out) :: that

    type(base_geometry_t), pointer :: subs

    PUSH_SUB(base_geometry_get_sub_geometry)

    nullify(subs, that)
    call base_geometry_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_geometry_get(subs, that)
    nullify(subs)
    
    POP_SUB(base_geometry_get_sub_geometry)
  end subroutine base_geometry_get_sub_geometry

  ! ---------------------------------------------------------
  subroutine base_geometry_set_geometry(this, that)
    type(base_geometry_t), intent(inout) :: this
    type(geometry_t),      intent(in)    :: that

    logical :: dflt
    integer :: ierr

    PUSH_SUB(base_geometry_set_geometry)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(.not.geo_intrf_assoc(this%igeo))
    call json_get(this%config, "default", dflt, ierr)
    if(ierr/=JSON_OK) dflt = .true.
    ASSERT(.not.dflt)
    call geo_intrf_set(this%igeo, that)
    ASSERT(geo_intrf_assoc(this%igeo))

    POP_SUB(base_geometry_set_geometry)
  end subroutine base_geometry_set_geometry

  ! ---------------------------------------------------------
  subroutine base_geometry_get_config(this, that)
    type(base_geometry_t), target, intent(in)  :: this
    type(json_object_t),  pointer, intent(out) :: that

    PUSH_SUB(base_geometry_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_geometry_get_config)
  end subroutine base_geometry_get_config

  ! ---------------------------------------------------------
  subroutine base_geometry_get_space(this, that)
    type(base_geometry_t), target, intent(in)  :: this
    type(space_t),        pointer, intent(out) :: that

    PUSH_SUB(base_geometry_get_space)

    nullify(that)
    if(associated(this%space)) that => this%space

    POP_SUB(base_geometry_get_space)
  end subroutine base_geometry_get_space

  ! ---------------------------------------------------------
  subroutine base_geometry_get_geometry(this, that)
    type(base_geometry_t),     intent(in)  :: this
    type(geometry_t), pointer, intent(out) :: that

    PUSH_SUB(base_geometry_get_geometry)

    nullify(that)
    if(geo_intrf_assoc(this%igeo))&
      call geo_intrf_get(this%igeo, that)

    POP_SUB(base_geometry_get_geometry)
  end subroutine base_geometry_get_geometry

  ! ---------------------------------------------------------
  subroutine base_geometry__copy__(this, that)
    type(base_geometry_t), intent(inout) :: this
    type(base_geometry_t), intent(in)    :: that

    type(refcount_t), pointer :: rcnt

    PUSH_SUB(base_geometry__copy__)

    rcnt => this%rcnt
    nullify(this%rcnt)
    call base_geometry__end__(this)
    if(associated(that%config).and.associated(that%space))then
      call base_geometry__init__(this, that%space, that%config)
      call refcount_del(this%rcnt)
      call geo_intrf_copy(this%igeo, that%igeo)
    end if
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(base_geometry__copy__)
  end subroutine base_geometry__copy__

  ! ---------------------------------------------------------
  subroutine base_geometry__end__(this)
    type(base_geometry_t), target, intent(inout) :: this

    PUSH_SUB(base_geometry__end__)

    nullify(this%config, this%space)
    if(associated(this%rcnt)) call refcount_del(this%rcnt)
    call geo_intrf_end(this%igeo)
    ASSERT(base_geometry_dict_len(this%dict)==0)
    call base_geometry_dict_end(this%dict)

    POP_SUB(base_geometry__end__)
  end subroutine base_geometry__end__

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator__init__(this, that)
    type(base_geometry_iterator_t), intent(out) :: this
    type(base_geometry_t),          intent(in)  :: that

    PUSH_SUB(base_geometry_iterator__init__)

    call geo_intrf_init(this%gitr, that%igeo)

    POP_SUB(base_geometry_iterator__init__)
  end subroutine base_geometry_iterator__init__

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator_next_atom(this, atom, ierr)
    type(base_geometry_iterator_t), intent(inout) :: this
    type(atom_t),          pointer, intent(out)   :: atom
    integer,              optional, intent(out)   :: ierr

    PUSH_SUB(base_geometry_iterator_next_atom)

    call geo_intrf_next(this%gitr, atom, ierr)

    POP_SUB(base_geometry_iterator_next_atom)
  end subroutine base_geometry_iterator_next_atom

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator_next_species(this, spec, ierr)
    type(base_geometry_iterator_t), intent(inout) :: this
    type(species_t),       pointer, intent(out)   :: spec
    integer,              optional, intent(out)   :: ierr

    PUSH_SUB(base_geometry_iterator_next_species)

    call geo_intrf_next(this%gitr, spec, ierr)

    POP_SUB(base_geometry_iterator_next_species)
  end subroutine base_geometry_iterator_next_species

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator_next_geometry(this, geometry, ierr)
    type(base_geometry_iterator_t), intent(inout) :: this
    type(geometry_t),      pointer, intent(out)   :: geometry
    integer,              optional, intent(out)   :: ierr

    type(base_geometry_t), pointer :: geom
    integer                        :: jerr

    PUSH_SUB(base_geometry_iterator_next_geometry)

    nullify(geometry, geom)
    call base_geometry_next(this, geom, ierr=jerr)
    if(associated(geom)) call base_geometry_get(geom, geometry)
    if(present(ierr)) ierr = jerr
    nullify(geom)

    POP_SUB(base_geometry_iterator_next_geometry)
  end subroutine base_geometry_iterator_next_geometry

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator__copy__(this, that)
    type(base_geometry_iterator_t), intent(out) :: this
    type(base_geometry_iterator_t), intent(in)  :: that

    PUSH_SUB(base_geometry_iterator__copy__)

    call geo_intrf_copy(this%gitr, that%gitr)

    POP_SUB(base_geometry_iterator__copy__)
  end subroutine base_geometry_iterator__copy__

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator__end__(this)
    type(base_geometry_iterator_t), intent(inout) :: this

    PUSH_SUB(base_geometry_iterator__end__)

    call geo_intrf_end(this%gitr)

    POP_SUB(base_geometry_iterator__end__)
  end subroutine base_geometry_iterator__end__

end module base_geometry_oct_m

#undef BASE_EXTENDED_ITERATOR_TYPE

!! Local Variables:
!! mode: f90
!! End:

