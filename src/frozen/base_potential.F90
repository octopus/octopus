#include "global.h"

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

#define BASE_LEAF_TYPE

module base_potential_oct_m

  use base_density_oct_m
  use base_states_oct_m
  use base_system_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m

#define BASE_TEMPLATE_NAME base_potential
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private
  
  public ::           &
    base_potential_t

  public ::                   &
    base_potential__init__,   &
    base_potential__start__,  &
    base_potential__acc__,    &
    base_potential__sub__,    &
    base_potential__calc__,   &
    base_potential__update__, &
    base_potential__reset__,  &
    base_potential__stop__,   &
    base_potential__copy__,   &
    base_potential__end__

  public ::                &
    base_potential_new,    &
    base_potential_del,    &
    base_potential_init,   &
    base_potential_start,  &
    base_potential_acc,    &
    base_potential_calc,   &
    base_potential_update, &
    base_potential_reset,  &
    base_potential_stop,   &
    base_potential_set,    &
    base_potential_get,    &
    base_potential_copy,   &
    base_potential_end

#define BASE_TEMPLATE_NAME base_potential
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  integer, parameter :: default_nspin = 1

  type :: base_potential_t
    private
    type(json_object_t),  pointer :: config =>null()
    type(base_system_t),  pointer :: sys    =>null()
    type(simulation_t),   pointer :: sim    =>null()
    type(refcount_t),     pointer :: rcnt   =>null()
    integer                       :: nspin  = default_nspin
    real(kind=wp)                 :: energy = 0.0_wp
    type(storage_t)               :: data
    type(base_potential_dict_t)   :: dict
  end type base_potential_t

  interface base_potential__init__
    module procedure base_potential__init__type
    module procedure base_potential__init__copy
  end interface base_potential__init__

  interface base_potential_new
    module procedure base_potential_new_type
    module procedure base_potential_new_pass
  end interface base_potential_new

  interface base_potential_init
    module procedure base_potential_init_type
  end interface base_potential_init

  interface base_potential_set
    module procedure base_potential_set_info
  end interface base_potential_set

  interface base_potential_get
    module procedure base_potential_get_info
    module procedure base_potential_get_energy
    module procedure base_potential_get_config
    module procedure base_potential_get_system
    module procedure base_potential_get_density
    module procedure base_potential_get_simulation
    module procedure base_potential_get_storage
    module procedure base_potential_get_potential_r1
    module procedure base_potential_get_potential_r2
  end interface base_potential_get

  interface base_potential_gets
    module procedure base_potential_gets_storage
    module procedure base_potential_gets_potential_r1
    module procedure base_potential_gets_potential_r2
  end interface base_potential_gets

contains

#define BASE_TEMPLATE_NAME base_potential
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  function base_potential_new_type(sys, config) result(this)
    type(base_system_t), intent(in)  :: sys
    type(json_object_t), intent(in)  :: config

    type(base_potential_t), pointer :: this

    PUSH_SUB(base_potential_new_type)

    this => base_potential_new(sys, config, base_potential_init_type)

    POP_SUB(base_potential_new_type)
  end function base_potential_new_type

  ! ---------------------------------------------------------
  function base_potential_new_pass(sys, config, init) result(this)
    type(base_system_t), intent(in) :: sys
    type(json_object_t), intent(in) :: config

    type(base_potential_t), pointer :: this
    
    interface
      subroutine init(this, sys, config)
        use json_oct_m
        use base_system_oct_m
        import :: base_potential_t
        type(base_potential_t), intent(out) :: this
        type(base_system_t),    intent(in)  :: sys
        type(json_object_t),    intent(in)  :: config
      end subroutine init
    end interface
    
    PUSH_SUB(base_potential_new_pass)

    nullify(this)
    SAFE_ALLOCATE(this)
    call init(this, sys, config)
    ASSERT(associated(this%rcnt))
    call refcount_set(this%rcnt, dynamic=.true.)

    POP_SUB(base_potential_new_pass)
  end function base_potential_new_pass

  ! ---------------------------------------------------------
  subroutine base_potential__init__type(this, sys, config)
    type(base_potential_t),      intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    logical                      :: uspn
    integer                      :: ierr

    PUSH_SUB(base_potential__init__type)

    nullify(cnfg)
    this%config => config
    this%sys => sys
    this%rcnt => refcount_new()
    this%nspin = default_nspin
    call json_get(this%config, "spin", uspn, ierr)
    if(ierr/=JSON_OK) uspn = .false.
    if(uspn) call base_system_get(this%sys, nspin=this%nspin)
    ASSERT(this%nspin>0)
    ASSERT(this%nspin<3)
    call json_get(this%config, "storage", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(cnfg, "full", .false.)
    call json_set(cnfg, "dimensions", this%nspin)
    call storage_init(this%data, cnfg)
    nullify(cnfg)
    call base_potential_dict_init(this%dict)

    POP_SUB(base_potential__init__type)
  end subroutine base_potential__init__type

  ! ---------------------------------------------------------
  subroutine base_potential__init__copy(this, that)
    type(base_potential_t), intent(out) :: this
    type(base_potential_t), intent(in)  :: that

    PUSH_SUB(base_potential__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call base_potential__init__(this, that%sys, that%config)
    if(associated(that%sim)) call base_potential__start__(this, that%sim)

    POP_SUB(base_potential__init__copy)
  end subroutine base_potential__init__copy

  ! ---------------------------------------------------------
  subroutine base_potential_init_type(this, sys, config)
    type(base_potential_t), intent(out) :: this
    type(base_system_t),    intent(in)  :: sys
    type(json_object_t),    intent(in)  :: config

    PUSH_SUB(base_potential_init_type)

    call base_potential__init__(this, sys, config)

    POP_SUB(base_potential_init_type)
  end subroutine base_potential_init_type

  ! ---------------------------------------------------------
  subroutine base_potential__start__(this, sim)
    type(base_potential_t),     intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_potential__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    call storage_start(this%data, this%sim)

    POP_SUB(base_potential__start__)
  end subroutine base_potential__start__

  ! ---------------------------------------------------------
  subroutine base_potential__acc__(this, that, config)
    type(base_potential_t),        intent(inout) :: this
    type(base_potential_t),        intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    PUSH_SUB(base_potential__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%nspin==that%nspin)
    if(present(config)) continue
    this%energy = this%energy + that%energy
    call storage_add(this%data, that%data)

    POP_SUB(base_potential__acc__)
  end subroutine base_potential__acc__

  ! ---------------------------------------------------------
  subroutine base_potential__sub__(this, that, config)
    type(base_potential_t),        intent(inout) :: this
    type(base_potential_t),        intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    PUSH_SUB(base_potential__sub__)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%nspin==that%nspin)
    if(present(config)) continue
    this%energy = this%energy - that%energy
    call storage_sub(this%data, that%data)

    POP_SUB(base_potential__sub__)
  end subroutine base_potential__sub__

  ! ---------------------------------------------------------
  subroutine base_potential__calc__(this)
    type(base_potential_t), intent(inout) :: this

    type(base_density_t), pointer :: dnst
    type(storage_t),      pointer :: data

    PUSH_SUB(base_potential__calc__)

    nullify(dnst, data)
    this%energy = 0.0_wp
    call base_potential_get(this, dnst)
    ASSERT(associated(dnst))
    call base_density_get(dnst, data)
    ASSERT(associated(data))
    nullify(dnst)
    call storage_integrate(this%data, data, this%energy)
    nullify(data)

    POP_SUB(base_potential__calc__)
  end subroutine base_potential__calc__

  ! ---------------------------------------------------------
  subroutine base_potential__update__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_update(this%data)

    POP_SUB(base_potential__update__)
  end subroutine base_potential__update__

  ! ---------------------------------------------------------
  subroutine base_potential__reset__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%energy = 0.0_wp
    call storage_reset(this%data)

    POP_SUB(base_potential__reset__)
  end subroutine base_potential__reset__

  ! ---------------------------------------------------------
  subroutine base_potential__stop__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call storage_stop(this%data)

    POP_SUB(base_potential__stop__)
  end subroutine base_potential__stop__

  ! ---------------------------------------------------------
  subroutine base_potential_start(this, sim)
    type(base_potential_t), intent(inout) :: this
    type(simulation_t),     intent(in)    :: sim

    PUSH_SUB(base_potential_start)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    call base_potential__apply__(this, start)

    POP_SUB(base_potential_start)

  contains

    subroutine start(this)
      type(base_potential_t), intent(inout) :: this

      PUSH_SUB(base_potential_start.start)

      call base_potential__start__(this, sim)

      POP_SUB(base_potential_start.start)
    end subroutine start

  end subroutine base_potential_start

  ! ---------------------------------------------------------
  subroutine base_potential_acc(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential_acc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    ASSERT(base_potential_dict_len(this%dict)>0)
    call base_potential__reset__(this)
    call base_potential__reduce__(this, base_potential__acc__)
    call base_potential__update__(this)

    POP_SUB(base_potential_acc)
  end subroutine base_potential_acc

  ! ---------------------------------------------------------
  subroutine base_potential_calc(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential_calc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_potential_acc(this)
    call base_potential__apply__(this, base_potential__calc__)

    POP_SUB(base_potential_calc)
  end subroutine base_potential_calc

  ! ---------------------------------------------------------
  recursive subroutine base_potential_update(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential_update)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_potential__apply__(this, base_potential__update__)

    POP_SUB(base_potential_update)
  end subroutine base_potential_update

  ! ---------------------------------------------------------
  subroutine base_potential_reset(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential_reset)

    call base_potential__apply__(this, base_potential__reset__)

    POP_SUB(base_potential_reset)
  end subroutine base_potential_reset

  ! ---------------------------------------------------------
  recursive subroutine base_potential_stop(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential_stop)

    call base_potential__apply__(this, base_potential__stop__)

    POP_SUB(base_potential_stop)
  end subroutine base_potential_stop

  ! ---------------------------------------------------------
  subroutine base_potential_gets_storage(this, name, that)
    type(base_potential_t),   intent(in)  :: this
    character(len=*),         intent(in)  :: name
    type(storage_t), pointer, intent(out) :: that

    type(base_potential_t), pointer :: subs

    PUSH_SUB(base_potential_gets_storage)

    nullify(that, subs)
    call base_potential_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_potential_get(subs, that)
    nullify(subs)

    POP_SUB(base_potential_gets_storage)
  end subroutine base_potential_gets_storage

  ! ---------------------------------------------------------
  subroutine base_potential_gets_potential_r1(this, name, that, spin)
    type(base_potential_t),               intent(in)  :: this
    character(len=*),                     intent(in)  :: name
    real(kind=wp), dimension(:), pointer, intent(out) :: that
    integer,                    optional, intent(in)  :: spin

    type(base_potential_t), pointer :: subs

    PUSH_SUB(base_potential_gets_potential_r1)

    nullify(that, subs)
    call base_potential_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_potential_get(subs, that, spin)

    POP_SUB(base_potential_gets_potential_r1)
  end subroutine base_potential_gets_potential_r1

  ! ---------------------------------------------------------
  subroutine base_potential_gets_potential_r2(this, name, that)
    type(base_potential_t),                 intent(in)  :: this
    character(len=*),                       intent(in)  :: name
    real(kind=wp), dimension(:,:), pointer, intent(out) :: that

    type(base_potential_t), pointer :: subs

    PUSH_SUB(base_potential_gets_potential_r2)

    nullify(that, subs)
    call base_potential_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_potential_get(subs, that)

    POP_SUB(base_potential_gets_potential_r2)
  end subroutine base_potential_gets_potential_r2

  ! ---------------------------------------------------------
  subroutine base_potential_set_info(this, static, energy)
    type(base_potential_t),  intent(inout) :: this
    logical,       optional, intent(in)    :: static
    real(kind=wp), optional, intent(in)    :: energy

    PUSH_SUB(base_potential_set_info)

    ASSERT(associated(this%config))
    call storage_set(this%data, lock=static)
    if(present(energy)) this%energy = energy

    POP_SUB(base_potential_set_info)
  end subroutine base_potential_set_info

  ! ---------------------------------------------------------
  subroutine base_potential_get_info(this, size, nspin, static, fine, use)
    type(base_potential_t), intent(in)  :: this
    integer,      optional, intent(out) :: size
    integer,      optional, intent(out) :: nspin
    logical,      optional, intent(out) :: static
    logical,      optional, intent(out) :: fine
    logical,      optional, intent(out) :: use

    PUSH_SUB(base_potential_get_info)

    ASSERT(associated(this%config))
    if(present(nspin)) nspin = this%nspin
    call storage_get(this%data, size=size, lock=static, fine=fine, alloc=use)

    POP_SUB(base_potential_get_info)
  end subroutine base_potential_get_info

  ! ---------------------------------------------------------
  subroutine base_potential_get_energy(this, energy)
    type(base_potential_t),  intent(in)  :: this
    real(kind=wp),           intent(out) :: energy

    PUSH_SUB(base_potential_get_energy)

    ASSERT(associated(this%config))
    energy = this%energy

    POP_SUB(base_potential_get_energy)
  end subroutine base_potential_get_energy

  ! ---------------------------------------------------------
  subroutine base_potential_get_config(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(json_object_t),   pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_potential_get_config)
  end subroutine base_potential_get_config

  ! ---------------------------------------------------------
  subroutine base_potential_get_system(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(base_system_t),   pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_system)

    nullify(that)
    if(associated(this%sys)) that => this%sys

    POP_SUB(base_potential_get_system)
  end subroutine base_potential_get_system

  ! ---------------------------------------------------------
  subroutine base_potential_get_density(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(base_density_t),  pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_density)

    nullify(that)
    if(associated(this%sys))&
      call base_system_get(this%sys, that)

    POP_SUB(base_potential_get_density)
  end subroutine base_potential_get_density

  ! ---------------------------------------------------------
  subroutine base_potential_get_simulation(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(simulation_t),    pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_potential_get_simulation)
  end subroutine base_potential_get_simulation

  ! ---------------------------------------------------------
  subroutine base_potential_get_storage(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(storage_t),       pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_storage)

    nullify(that)
    if(associated(this%config)) that => this%data

    POP_SUB(base_potential_get_storage)
  end subroutine base_potential_get_storage

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_r1(this, that, spin)
    type(base_potential_t),               intent(in)  :: this
    real(kind=wp), dimension(:), pointer, intent(out) :: that
    integer,                    optional, intent(in)  :: spin

    integer :: ispn

    PUSH_SUB(base_potential_get_potential_r1)

    nullify(that)
    ispn = 1
    if(present(spin)) ispn = spin
    call storage_get(this%data, that, ispn)

    POP_SUB(base_potential_get_potential_r1)
  end subroutine base_potential_get_potential_r1

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_r2(this, that)
    type(base_potential_t),                 intent(in)  :: this
    real(kind=wp), dimension(:,:), pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_potential_r2)

    nullify(that)
    call storage_get(this%data, that)

    POP_SUB(base_potential_get_potential_r2)
  end subroutine base_potential_get_potential_r2

  ! ---------------------------------------------------------
  subroutine base_potential__copy__(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that

    type(refcount_t), pointer :: rcnt

    PUSH_SUB(base_potential__copy__)

    rcnt => this%rcnt
    nullify(this%rcnt)
    call base_potential__end__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_potential__init__(this, that)
      call refcount_del(this%rcnt)
      this%energy = that%energy
      if(associated(that%sim)) call storage_copy(this%data, that%data)
    end if
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(base_potential__copy__)
  end subroutine base_potential__copy__

  ! ---------------------------------------------------------
  subroutine base_potential__end__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__end__)

    nullify(this%config, this%sys, this%sim)
    if(associated(this%rcnt)) call refcount_del(this%rcnt)
    this%nspin = default_nspin
    this%energy = 0.0_wp
    call storage_end(this%data)
    ASSERT(base_potential_dict_len(this%dict)==0)
    call base_potential_dict_end(this%dict)

    POP_SUB(base_potential__end__)
  end subroutine base_potential__end__

end module base_potential_oct_m

#undef BASE_LEAF_TYPE

!! Local Variables:
!! mode: f90
!! End:
