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
    base_states__update__, &
    base_states__reset__,  &
    base_states__stop__,   &
    base_states__copy__,   &
    base_states__end__

  public ::             &
    base_states_new,    &
    base_states_del,    &
    base_states_init,   &
    base_states_notify, &
    base_states_start,  &
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
    type(refcount_t),    pointer :: rcnt   =>null()
    type(base_density_t)         :: density
    type(base_states_dict_t)     :: dict
  end type base_states_t

  interface base_states__init__
    module procedure base_states__init__type
    module procedure base_states__init__copy
  end interface base_states__init__

  interface base_states_new
    module procedure base_states_new_type
    module procedure base_states_new_pass
  end interface base_states_new

  interface base_states_init
    module procedure base_states_init_type
  end interface base_states_init

  interface base_states_notify
    module procedure base_states_notify_type
    module procedure base_states_notify_subs
  end interface base_states_notify

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

contains
    
#define BASE_TEMPLATE_NAME base_states
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  function base_states_new_type(config) result(this)
    type(json_object_t), intent(in) :: config

    type(base_states_t), pointer :: this
    
    PUSH_SUB(base_states_new_type)

    this => base_states_new(config, base_states_init_type)

    POP_SUB(base_states_new_type)
  end function base_states_new_type

  ! ---------------------------------------------------------
  function base_states_new_pass(config, init) result(this)
    type(json_object_t), intent(in) :: config

    interface
      subroutine init(this, config)
        use json_oct_m
        import :: base_states_t
        type(base_states_t), intent(out) :: this
        type(json_object_t), intent(in)  :: config
      end subroutine init
    end interface
    
    type(base_states_t), pointer :: this
    
    PUSH_SUB(base_states_new_pass)

    nullify(this)
    SAFE_ALLOCATE(this)
    call init(this, config)
    ASSERT(associated(this%rcnt))
    call refcount_set(this%rcnt, dynamic=.true.)

    POP_SUB(base_states_new_pass)
  end function base_states_new_pass

  ! ---------------------------------------------------------
  subroutine base_states__iinit__(this, config)
    type(base_states_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(base_states__iinit__)

    this%config => config
    this%rcnt => refcount_new()
    call base_states_dict_init(this%dict)

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
  subroutine base_states_notify_type(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states_notify_type)

    call base_density_notify(this%density)
    
    POP_SUB(base_states_notify_type)
  end subroutine base_states_notify_type
    
  ! ---------------------------------------------------------
  subroutine base_states_notify_subs(this, name)
    type(base_states_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name

    type(base_states_t), pointer :: subs
    
    PUSH_SUB(base_states_notify_subs)

    nullify(subs)
    call base_states_gets(this, trim(adjustl(name)), subs)
    ASSERT(associated(subs))
    call base_states_notify(subs)
    nullify(subs)
    
    POP_SUB(base_states_notify_subs)
  end subroutine base_states_notify_subs
    
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
  subroutine base_states__sets__(this, name, that, config, lock, active)
    type(base_states_t),           intent(inout) :: this
    character(len=*),              intent(in)    :: name
    type(base_states_t), optional, intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config
    logical,             optional, intent(in)    :: lock
    logical,             optional, intent(in)    :: active

    PUSH_SUB(base_states__sets__)

    ASSERT(associated(this%config))
    ASSERT(len_trim(adjustl(name))>0)
    if(present(that))then
      ASSERT(associated(that%config))
      call base_density_sets(this%density, trim(adjustl(name)), that%density, config=config, lock=lock, active=active)
    else
      call base_density_sets(this%density, trim(adjustl(name)), config=config, lock=lock, active=active)
    end if

    POP_SUB(base_states__sets__)
  end subroutine base_states__sets__
    
  ! ---------------------------------------------------------
  subroutine base_states__dels__(this, name, that)
    type(base_states_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__dels__)

    ASSERT(associated(this%config))
    ASSERT(len_trim(name)>0)
    ASSERT(associated(that%config))
    call base_density_dels(this%density, trim(adjustl(name)), that%density)

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

    nullify(that)
    if(associated(this%config)) that => this%density

    POP_SUB(base_states_get_density)
  end subroutine base_states_get_density
    
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
  subroutine base_states_get_density_1d(this, that, spin, total)
    type(base_states_t),          intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    integer,            optional, intent(in) :: spin
    logical,            optional, intent(in) :: total

    type(base_density_t), pointer :: dnst

    PUSH_SUB(base_states_get_density_1d)

    nullify(that, dnst)
    call base_states_get(this, dnst)
    if(associated(dnst))&
      call base_density_get(dnst, that, spin, total)
    nullify(dnst)

    POP_SUB(base_states_get_density_1d)
  end subroutine base_states_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_states_get_density_2d(this, that, total)
    type(base_states_t),            intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    type(base_density_t), pointer :: dnst

    PUSH_SUB(base_states_get_density_2d)

    nullify(that, dnst)
    call base_states_get(this, dnst)
    if(associated(dnst))&
      call base_density_get(dnst, that, total)
    nullify(dnst)

    POP_SUB(base_states_get_density_2d)
  end subroutine base_states_get_density_2d

  ! ---------------------------------------------------------
  subroutine base_states__copy__(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    type(refcount_t), pointer :: rcnt

    PUSH_SUB(base_states__copy__)

    rcnt => this%rcnt
    nullify(this%rcnt)
    call base_states__end__(this)
    if(associated(that%config))then
      call base_states__iinit__(this, that%config)
      call refcount_del(this%rcnt)
      this%sim => that%sim
      call base_density__copy__(this%density, that%density)
    end if
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(base_states__copy__)
  end subroutine base_states__copy__
    
  ! ---------------------------------------------------------
  subroutine base_states__end__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__end__)

    nullify(this%config, this%sim)
    if(associated(this%rcnt)) call refcount_del(this%rcnt)
    call base_density__end__(this%density)
    ASSERT(base_states_dict_len(this%dict)==0)
    call base_states_dict_end(this%dict)

    POP_SUB(base_states__end__)
  end subroutine base_states__end__
    
end module base_states_oct_m

!! Local Variables:
!! mode: f90
!! End:
