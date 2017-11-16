#include "global.h"

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

module base_states_oct_m

  use base_density_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m

#define BASE_TEMPLATE_NAME base_states
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private

  public ::        &
    base_states_t

  public ::                &
    base_states__init__,   &
    base_states__start__,  &
    base_states__acc__,    &
    base_states__update__, &
    base_states__reset__,  &
    base_states__stop__,   &
    base_states__copy__,   &
    base_states__end__

  public ::             &
    base_states_new,    &
    base_states_del,    &
    base_states_init,   &
    base_states_start,  &
    base_states_acc,    &
    base_states_update, &
    base_states_reset,  &
    base_states_stop,   &
    base_states_get,    &
    base_states_copy,   &
    base_states_end

#define BASE_TEMPLATE_NAME base_states
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  type :: base_states_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(base_states_t), pointer :: prnt   =>null()
    type(base_density_t)         :: density
    type(base_states_dict_t)     :: dict
    type(base_states_list_t)     :: list
  end type base_states_t

  interface base_states__init__
    module procedure base_states__init__type
    module procedure base_states__init__copy
  end interface base_states__init__

  interface base_states_init
    module procedure base_states_init_type
    module procedure base_states_init_copy
  end interface base_states_init

  interface base_states_get
    module procedure base_states_get_info
    module procedure base_states_get_config
    module procedure base_states_get_simulation
    module procedure base_states_get_charge
    module procedure base_states_get_density
    module procedure base_states_get_density_1d
    module procedure base_states_get_density_2d
  end interface base_states_get

  interface base_states_gets
    module procedure base_states_gets_density
    module procedure base_states_gets_density_1d
    module procedure base_states_gets_density_2d
  end interface base_states_gets

  interface base_states_copy
    module procedure base_states_copy_type
  end interface base_states_copy

contains
    
#define BASE_TEMPLATE_NAME base_states
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  subroutine base_states__new__(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states__new__)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(base_states__new__)
  end subroutine base_states__new__

  ! ---------------------------------------------------------
  subroutine base_states__del__(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states__del__)

    if(associated(this))then
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(base_states__del__)
  end subroutine base_states__del__

  ! ---------------------------------------------------------
  subroutine base_states_new(this, that)
    type(base_states_t),  target, intent(inout) :: this
    type(base_states_t), pointer                :: that

    PUSH_SUB(base_states_new)

    nullify(that)
    call base_states__new__(that)
    that%prnt => this
    call base_states_list_push(this%list, that)

    POP_SUB(base_states_new)
  end subroutine base_states_new

  ! ---------------------------------------------------------
  subroutine base_states__iinit__(this, config)
    type(base_states_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(base_states__iinit__)

    this%config => config
    call base_states_dict_init(this%dict)
    call base_states_list_init(this%list)

    POP_SUB(base_states__iinit__)
  end subroutine base_states__iinit__
    
  ! ---------------------------------------------------------
  subroutine base_states__init__type(this, config)
    type(base_states_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_states__init__type)

    nullify(cnfg)
    call base_states__iinit__(this, config)
    call json_get(this%config, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_density__init__(this%density, cnfg)
    nullify(cnfg)

    POP_SUB(base_states__init__type)
  end subroutine base_states__init__type
    
  ! ---------------------------------------------------------
  subroutine base_states__init__copy(this, that, start)
    type(base_states_t), intent(out) :: this
    type(base_states_t), intent(in)  :: that
    logical,   optional, intent(in)  :: start

    logical :: istr

    PUSH_SUB(base_states__init__copy)

    ASSERT(associated(that%config))
    call base_states__iinit__(this, that%config)
    call base_density__init__(this%density, that%density)
    istr = .true.
    if(present(start)) istr = start
    if(istr)then
      if(present(start))then
        ASSERT(associated(that%sim))
      end if
      if(associated(that%sim)) call base_states__start__(this, that%sim)
    end if

    POP_SUB(base_states__init__copy)
  end subroutine base_states__init__copy

  ! ---------------------------------------------------------
  subroutine base_states_init_type(this, config)
    type(base_states_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_states_init_type)

    call base_states__init__(this, config)

    POP_SUB(base_states_init_type)
  end subroutine base_states_init_type
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_init_copy(this, that)
    type(base_states_t), intent(out) :: this
    type(base_states_t), intent(in)  :: that

    PUSH_SUB(base_states_init_copy)

    call base_states__init__(this, that)
    call base_states__build__(this, init)
    
    POP_SUB(base_states_init_copy)

  contains

    recursive subroutine init(this, name, isub)
      type(base_states_t), intent(inout) :: this
      character(len=*),    intent(in)    :: name
      type(base_states_t), intent(in)    :: isub

      type(base_states_t), pointer :: osub
      
      POP_SUB(base_states_init_copy.init)

      nullify(osub)
      if(base_states_list_index(that%list, isub)>0)then
        call base_states_new(this, osub)
        call base_states_init(osub, isub)
        call base_states_sets(this, name, osub)
        nullify(osub)
      else
        call base_states_sets(this, name, isub)
      end if

      PUSH_SUB(base_states_init_copy.init)
    end subroutine init

  end subroutine base_states_init_copy

  ! ---------------------------------------------------------
  subroutine base_states__start__(this, sim)
    type(base_states_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_states__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    call base_density__start__(this%density, sim)

    POP_SUB(base_states__start__)
  end subroutine base_states__start__
    
  ! ---------------------------------------------------------
  subroutine base_states__acc__(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%config))
    ASSERT(associated(that%sim))
    call base_density__acc__(this%density, that%density)

    POP_SUB(base_states__acc__)
  end subroutine base_states__acc__

  ! ---------------------------------------------------------
  subroutine base_states__update__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_density__update__(this%density)

    POP_SUB(base_states__update__)
  end subroutine base_states__update__

  ! ---------------------------------------------------------
  subroutine base_states__reset__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_density__reset__(this%density)

    POP_SUB(base_states__reset__)
  end subroutine base_states__reset__
    
  ! ---------------------------------------------------------
  subroutine base_states__stop__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call base_density__stop__(this%density)

    POP_SUB(base_states__stop__)
  end subroutine base_states__stop__
    
  ! ---------------------------------------------------------
  subroutine base_states_start(this, sim)
    type(base_states_t), intent(inout) :: this
    type(simulation_t),  intent(in)    :: sim

    PUSH_SUB(base_states_start)

    call base_states__apply__(this, start)
    
    POP_SUB(base_states_start)
    
  contains

    subroutine start(this)
      type(base_states_t), intent(inout) :: this

      PUSH_SUB(base_states_start.start)
      
      call base_states__start__(this, sim)

      POP_SUB(base_states_start.start)
    end subroutine start

  end subroutine base_states_start

  ! ---------------------------------------------------------
  recursive subroutine base_states_acc(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states_acc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(base_states_dict_len(this%dict)>0)then
      call base_states__reset__(this)
      call base_states__reduce__(this, base_states__acc__)
      call base_states__update__(this)
    end if
    
    POP_SUB(base_states_acc)
  end subroutine base_states_acc

  ! ---------------------------------------------------------
  subroutine base_states_update(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states_update)

    call base_states__apply__(this, base_states__update__)
    
    POP_SUB(base_states_update)
  end subroutine base_states_update

  ! ---------------------------------------------------------
  subroutine base_states_reset(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states_reset)

    call base_states__apply__(this, base_states__reset__)
    
    POP_SUB(base_states_reset)
  end subroutine base_states_reset

  ! ---------------------------------------------------------
  subroutine base_states_stop(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states_stop)

    call base_states__apply__(this, base_states__stop__)
    
    POP_SUB(base_states_stop)
  end subroutine base_states_stop

  ! ---------------------------------------------------------
  subroutine base_states__sets__(this, name, that)
    type(base_states_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__sets__)

    call base_density_sets(this%density, trim(adjustl(name)), that%density)

    POP_SUB(base_states__sets__)
  end subroutine base_states__sets__
    
  ! ---------------------------------------------------------
  subroutine base_states__dels__(this, name, ierr)
    type(base_states_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    integer,             intent(out)   :: ierr

    PUSH_SUB(base_states__dels__)

    call base_density_dels(this%density, trim(adjustl(name)), ierr)

    POP_SUB(base_states__dels__)
  end subroutine base_states__dels__
    
  ! ---------------------------------------------------------
  subroutine base_states_gets_density(this, name, that)
    type(base_states_t),   intent(in) :: this
    character(len=*),      intent(in) :: name
    type(base_density_t), pointer     :: that

    type(base_states_t), pointer :: subs

    PUSH_SUB(base_states_gets_density)

    nullify(that, subs)
    call base_states_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_states_get(subs, that)

    POP_SUB(base_states_gets_density)
  end subroutine base_states_gets_density
    
  ! ---------------------------------------------------------
  subroutine base_states_gets_density_1d(this, name, that, spin, total)
    type(base_states_t),          intent(in) :: this
    character(len=*),             intent(in) :: name
    real(kind=wp), dimension(:), pointer     :: that
    integer,            optional, intent(in) :: spin
    logical,            optional, intent(in) :: total

    type(base_density_t), pointer :: dnst

    PUSH_SUB(base_states_gets_density_1d)

    nullify(that, dnst)
    call base_states_gets(this, trim(adjustl(name)), dnst)
    if(associated(dnst)) call base_density_get(dnst, that, spin, total)

    POP_SUB(base_states_gets_density_1d)
  end subroutine base_states_gets_density_1d

  ! ---------------------------------------------------------
  subroutine base_states_gets_density_2d(this, name, that, total)
    type(base_states_t),            intent(in) :: this
    character(len=*),               intent(in) :: name
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    type(base_density_t), pointer :: dnst

    PUSH_SUB(base_states_gets_density_2d)

    nullify(that, dnst)
    call base_states_gets(this, trim(adjustl(name)), dnst)
    if(associated(dnst)) call base_density_get(dnst, that, total)

    POP_SUB(base_states_gets_density_2d)
  end subroutine base_states_gets_density_2d

  ! ---------------------------------------------------------
  subroutine base_states_get_info(this, nspin)
    type(base_states_t), intent(in)  :: this
    integer,   optional, intent(out) :: nspin

    PUSH_SUB(base_states_get_info)

    call base_density_get(this%density, nspin=nspin)

    POP_SUB(base_states_get_info)
  end subroutine base_states_get_info
    
  ! ---------------------------------------------------------
  subroutine base_states_get_charge(this, charge, spin, total)
    type(base_states_t), intent(in)  :: this
    real(kind=wp),       intent(out) :: charge
    integer,   optional, intent(out) :: spin
    logical,   optional, intent(out) :: total

    PUSH_SUB(base_states_get_charge)

    call base_density_get(this%density, charge, spin=spin, total=total)

    POP_SUB(base_states_get_charge)
  end subroutine base_states_get_charge
    
  ! ---------------------------------------------------------
  subroutine base_states_get_config(this, that)
    type(base_states_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_states_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_states_get_config)
  end subroutine base_states_get_config
    
  ! ---------------------------------------------------------
  subroutine base_states_get_simulation(this, that)
    type(base_states_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(base_states_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_states_get_simulation)
  end subroutine base_states_get_simulation
    
  ! ---------------------------------------------------------
  subroutine base_states_get_density(this, that)
    type(base_states_t),   target, intent(in) :: this
    type(base_density_t), pointer             :: that

    PUSH_SUB(base_states_get_density)

    that => this%density

    POP_SUB(base_states_get_density)
  end subroutine base_states_get_density
    
  ! ---------------------------------------------------------
  subroutine base_states_get_density_1d(this, that, spin, total)
    type(base_states_t),          intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    integer,            optional, intent(in) :: spin
    logical,            optional, intent(in) :: total

    PUSH_SUB(base_states_get_density_1d)

    nullify(that)
    call base_density_get(this%density, that, spin, total)

    POP_SUB(base_states_get_density_1d)
  end subroutine base_states_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_states_get_density_2d(this, that, total)
    type(base_states_t),            intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    PUSH_SUB(base_states_get_density_2d)

    nullify(that)
    call base_density_get(this%density, that, total)

    POP_SUB(base_states_get_density_2d)
  end subroutine base_states_get_density_2d

  ! ---------------------------------------------------------
  subroutine base_states__copy__(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__copy__)

    call base_states__end__(this)
    if(associated(that%config))then
      call base_states__iinit__(this, that%config)
      this%sim => that%sim
      call base_density__copy__(this%density, that%density)
    end if

    POP_SUB(base_states__copy__)
  end subroutine base_states__copy__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_copy_type(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states_copy_type)

    call base_states_end(this)
    call base_states__copy__(this, that)
    call base_states__build__(this, copy)

    POP_SUB(base_states_copy_type)
    
  contains

    recursive subroutine copy(this, name, isub)
      type(base_states_t), intent(inout) :: this
      character(len=*),    intent(in)    :: name
      type(base_states_t), intent(in)    :: isub

      type(base_states_t), pointer :: osub
      
      POP_SUB(base_states_copy_type.copy)

      nullify(osub)
      if(base_states_list_index(that%list, isub)>0)then
        call base_states_new(this, osub)
        call base_states_copy(osub, isub)
        call base_states_sets(this, name, osub)
        nullify(osub)
      else
        call base_states_sets(this, name, isub)
      end if

      PUSH_SUB(base_states_copy_type.copy)
    end subroutine copy

  end subroutine base_states_copy_type

  ! ---------------------------------------------------------
  subroutine base_states__end__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__end__)

    nullify(this%config, this%sim, this%prnt)
    call base_density__end__(this%density)
    call base_states_dict_end(this%dict)
    ASSERT(base_states_list_len(this%list)==0)
    call base_states_list_end(this%list)

    POP_SUB(base_states__end__)
  end subroutine base_states__end__
    
end module base_states_oct_m

!! Local Variables:
!! mode: f90
!! End:
