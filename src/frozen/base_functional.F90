#include "global.h"

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

module base_functional_oct_m

  use base_density_oct_m
  use base_states_oct_m
  use base_system_oct_m
  use functional_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use message_oct_m
  use messages_oct_m
  use msgbus_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m

#define BASE_TEMPLATE_NAME base_functional
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private

  public ::            &
    base_functional_t

  public ::                    &
    base_functional__init__,   &
    base_functional__start__,  &
    base_functional__acc__,    &
    base_functional__sub__,    &
    base_functional__calc__,   &
    base_functional__update__, &
    base_functional__reset__,  &
    base_functional__stop__,   &
    base_functional__copy__,   &
    base_functional__end__

  public ::                 &
    base_functional_new,    &
    base_functional_del,    &
    base_functional_init,   &
    base_functional_start,  &
    base_functional_acc,    &
    base_functional_update, &
    base_functional_reset,  &
    base_functional_stop,   &
    base_functional_set,    &
    base_functional_get,    &
    base_functional_copy,   &
    base_functional_end

#define BASE_TEMPLATE_NAME base_functional
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  integer, parameter :: default_nspin = 1

  type :: base_functional_t
    private
    type(json_object_t),  pointer :: config  =>null()
    type(base_system_t),  pointer :: sys     =>null()
    type(base_density_t), pointer :: density =>null()
    type(simulation_t),   pointer :: sim     =>null()
    type(refcount_t),     pointer :: rcnt    =>null()
    integer                       :: nspin   = default_nspin
    real(kind=wp)                 :: factor  = 1.0_wp
    real(kind=wp)                 :: energy  = 0.0_wp
    type(functional_t)            :: funct
    type(storage_t)               :: data
    type(msgbus_t)                :: msgb
    type(base_functional_dict_t)  :: dict
  end type base_functional_t

  interface base_functional__init__
    module procedure base_functional__init__type
    module procedure base_functional__init__copy
  end interface base_functional__init__

  interface base_functional_new
    module procedure base_functional_new_type
    module procedure base_functional_new_pass
  end interface base_functional_new

  interface base_functional_init
    module procedure base_functional_init_type
  end interface base_functional_init

  interface base_functional_set
    module procedure base_functional_set_info
  end interface base_functional_set

  interface base_functional_get
    module procedure base_functional_get_info
    module procedure base_functional_get_config
    module procedure base_functional_get_system
    module procedure base_functional_get_density
    module procedure base_functional_get_simulation
    module procedure base_functional_get_storage
    module procedure base_functional_get_functional_1d
    module procedure base_functional_get_functional_md
  end interface base_functional_get

contains

#define BASE_TEMPLATE_NAME base_functional
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  function base_functional_new_type(sys, config) result(this)
    type(base_system_t), intent(in) :: sys
    type(json_object_t), intent(in) :: config

    type(base_functional_t), pointer :: this

    PUSH_SUB(base_functional_new_type)

    this => base_functional_new(sys, config, base_functional_init_type)

    POP_SUB(base_functional_new_type)
  end function base_functional_new_type

  ! ---------------------------------------------------------
  function base_functional_new_pass(sys, config, init) result(this)
    type(base_system_t), intent(in) :: sys
    type(json_object_t), intent(in) :: config

    interface
      subroutine init(this, sys, config)
        use json_oct_m
        use base_system_oct_m
        import :: base_functional_t
        type(base_functional_t), intent(out) :: this
        type(base_system_t),     intent(in)  :: sys
        type(json_object_t),     intent(in)  :: config
      end subroutine init
    end interface

    type(base_functional_t), pointer :: this

    PUSH_SUB(base_functional_new_pass)

    nullify(this)
    SAFE_ALLOCATE(this)
    call init(this, sys, config)
    ASSERT(associated(this%rcnt))
    call refcount_set(this%rcnt, dynamic=.true.)

    POP_SUB(base_functional_new_pass)
  end function base_functional_new_pass

  ! ---------------------------------------------------------
  subroutine base_functional__init__type(this, sys, config)
    type(base_functional_t),     intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    type(storage_t),     pointer :: dnst
    integer                      :: ierr
    logical                      :: plrz, fuse

    PUSH_SUB(base_functional__init__type)

    nullify(cnfg, dnst)
    this%config => config
    this%sys => sys
    this%rcnt => refcount_new()
    call base_system_get(this%sys, this%density)
    ASSERT(associated(this%density))
    call base_density_get(this%density, nspin=this%nspin)
    ASSERT(this%nspin>0)
    ASSERT(this%nspin<3)
    call json_get(this%config, "factor", this%factor, ierr)
    if(ierr/=JSON_OK) this%factor = 1.0_wp
    plrz = .false.
    if(this%nspin>1)then
      call json_get(this%config, "spin", plrz, ierr)
      if(ierr/=JSON_OK) plrz = .true.
    end if
    call json_get(this%config, "functional", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(cnfg, "polarized", plrz)
    call base_density_get(this%density, dnst, total=(.not.plrz))
    ASSERT(associated(dnst))
    call functional_init(this%funct, dnst, cnfg)
    nullify(dnst, cnfg)
    call json_get(this%config, "storage", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(cnfg, "full", .false.)
    if(plrz) call json_set(cnfg, "dimensions", this%nspin)
    call functional_get(this%funct, use=fuse)
    if(.not.fuse) call json_set(cnfg, "allocate", .false.)
    call storage_init(this%data, cnfg)
    nullify(cnfg)
    call msgbus_init(this%msgb, number=2)
    call base_functional_dict_init(this%dict)

    POP_SUB(base_functional__init__type)
  end subroutine base_functional__init__type

  ! ---------------------------------------------------------
  subroutine base_functional__init__copy(this, that)
    type(base_functional_t), intent(out) :: this
    type(base_functional_t), intent(in)  :: that

    PUSH_SUB(base_functional__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call base_functional__init__(this, that%sys, that%config)
    if(associated(that%sim)) call base_functional__start__(this, that%sim)

    POP_SUB(base_functional__init__copy)
  end subroutine base_functional__init__copy

  ! ---------------------------------------------------------
  subroutine base_functional_init_type(this, sys, config)
    type(base_functional_t), intent(out) :: this
    type(base_system_t),     intent(in)  :: sys
    type(json_object_t),     intent(in)  :: config

    PUSH_SUB(base_functional_init_type)

    call base_functional__init__(this, sys, config)

    POP_SUB(base_functional_init_type)
  end subroutine base_functional_init_type

  ! ---------------------------------------------------------
  subroutine base_functional__start__(this, sim)
    type(base_functional_t),    intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    type(msgbus_t), pointer :: msgb

    PUSH_SUB(base_functional__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    nullify(msgb)
    this%sim => sim
    call functional_start(this%funct, this%sim)
    call storage_start(this%data, this%sim)
    call base_density_get(this%density, msgb)
    ASSERT(associated(msgb))
    call msgbus_attach(this%msgb, msgb, id=2)
    nullify(msgb)

    POP_SUB(base_functional__start__)
  end subroutine base_functional__start__

  ! ---------------------------------------------------------
  subroutine base_functional__acc__(this, that, config)
    type(base_functional_t),       intent(inout) :: this
    type(base_functional_t),       intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    PUSH_SUB(base_functional__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(present(config)) continue
    this%energy = this%energy + that%energy
    call storage_add(this%data, that%data)

    POP_SUB(base_functional__acc__)
  end subroutine base_functional__acc__

  ! ---------------------------------------------------------
  subroutine base_functional__sub__(this, that, config)
    type(base_functional_t),       intent(inout) :: this
    type(base_functional_t),       intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    PUSH_SUB(base_functional__sub__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(present(config)) continue
    this%energy = this%energy - that%energy
    call storage_sub(this%data, that%data)

    POP_SUB(base_functional__sub__)
  end subroutine base_functional__sub__

  ! ---------------------------------------------------------
  subroutine base_functional__calc__(this)
    type(base_functional_t), intent(inout) :: this

    logical :: fuse

    PUSH_SUB(base_functional__calc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_functional__reset__(this)
    call base_functional_get(this, use=fuse)
    if(fuse.and.(abs(this%factor)>tiny(this%factor)))then
      call functional_calc(this%funct, this%energy, this%data)
      if(abs(abs(this%factor)-1.0_wp)>epsilon(this%factor))then
        this%energy = this%factor * this%energy
        call storage_mlt(this%data, this%factor)
      end if
      call base_functional__update__(this)
    end if

    POP_SUB(base_functional__calc__)
  end subroutine base_functional__calc__

  ! ---------------------------------------------------------
  subroutine base_functional__update__(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_update(this%data)

    POP_SUB(base_functional__update__)
  end subroutine base_functional__update__

  ! ---------------------------------------------------------
  subroutine base_functional__reset__(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%energy = 0.0_wp
    call storage_reset(this%data)

    POP_SUB(base_functional__reset__)
  end subroutine base_functional__reset__

  ! ---------------------------------------------------------
  subroutine base_functional__stop__(this)
    type(base_functional_t), intent(inout) :: this

    type(msgbus_t), pointer :: msgb

    PUSH_SUB(base_functional__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim, msgb)
    call functional_stop(this%funct)
    call storage_stop(this%data)
    call base_density_get(this%density, msgb)
    ASSERT(associated(msgb))
    call msgbus_detach(this%msgb, msgb, id=2)
    nullify(msgb)

    POP_SUB(base_functional__stop__)
  end subroutine base_functional__stop__

  ! ---------------------------------------------------------
  subroutine base_functional__recv__(this, update, channel)
    type(base_functional_t), intent(in)  :: this
    logical,                 intent(out) :: update
    integer,       optional, intent(in)  :: channel

    type(msgbus_iterator_t)      :: iter
    type(json_object_t), pointer :: data
    type(message_t),     pointer :: mssg
    logical                      :: updt
    integer                      :: chid, ierr

    PUSH_SUB(base_functional__recv__)

    chid = 1
    if(present(channel)) chid = channel
    update = .false.
    call msgbus_init(iter, this%msgb, id=chid)
    do
      nullify(data, mssg)
      call msgbus_next(iter, mssg)
      if(.not.associated(mssg))exit
      call message_get(mssg, data)
      ASSERT(associated(data))
      call json_get(data, "update", updt, ierr)
      if(ierr==JSON_OK)then
        update = (update .or. updt)
        call msgbus_remove(iter, ierr=ierr)
        ASSERT(ierr==MSGBUS_OK)
      end if
    end do
    nullify(data, mssg)
    
    POP_SUB(base_functional__recv__)
  end subroutine base_functional__recv__
    
  ! ---------------------------------------------------------
  subroutine base_functional_start(this, sim)
    type(base_functional_t), intent(inout) :: this
    type(simulation_t),      intent(in)    :: sim

    PUSH_SUB(base_functional_start)

    call base_functional__apply__(this, start)

    POP_SUB(base_functional_start)
    
  contains

    subroutine start(this)
      type(base_functional_t), intent(inout) :: this

      PUSH_SUB(base_functional_start.start)
      
      call base_functional__start__(this, sim)

      POP_SUB(base_functional_start.start)
    end subroutine start

  end subroutine base_functional_start

  ! ---------------------------------------------------------
  subroutine base_functional_acc(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional_acc)

    ASSERT(base_functional_dict_len(this%dict)>0)
    call base_functional__reset__(this)
    call base_functional__reduce__(this, base_functional__acc__)
    call base_functional__update__(this)
    
    POP_SUB(base_functional_acc)
  end subroutine base_functional_acc

  ! ---------------------------------------------------------
  recursive subroutine base_functional_update(this)
    type(base_functional_t), intent(inout) :: this

    logical :: updt

    PUSH_SUB(base_functional_update)

    ASSERT(associated(this%config))
    call base_functional__recv__(this, updt, channel=2)
    if(updt)then
      call base_density_update(this%density)
      call base_functional__apply__(this, base_functional_update, parent=.false.)
      call base_functional__calc__(this)
    end if

    POP_SUB(base_functional_update)
  end subroutine base_functional_update

  ! ---------------------------------------------------------
  subroutine base_functional_reset(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional_reset)

    call base_functional__apply__(this, base_functional__reset__)

    POP_SUB(base_functional_reset)
  end subroutine base_functional_reset

  ! ---------------------------------------------------------
  subroutine base_functional_stop(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional_stop)

    call base_functional__apply__(this, base_functional__stop__)

    POP_SUB(base_functional_stop)
  end subroutine base_functional_stop

  ! ---------------------------------------------------------
  subroutine base_functional__sets__(this, name, that, config, lock, active)
    type(base_functional_t),           intent(inout) :: this
    character(len=*),                  intent(in)    :: name
    type(base_functional_t), optional, intent(in)    :: that
    type(json_object_t),     optional, intent(in)    :: config
    logical,                 optional, intent(in)    :: lock
    logical,                 optional, intent(in)    :: active

    PUSH_SUB(base_functional__sets__)

    ASSERT(associated(this%config))
    ASSERT(len_trim(name)>0)
    ASSERT(associated(that%config))
    if(present(config)) continue
    if(present(lock)) continue
    if(present(active)) continue
    if(present(that)) call msgbus_attach(this%msgb, that%msgb)
    
    POP_SUB(base_functional__sets__)
  end subroutine base_functional__sets__
    
  ! ---------------------------------------------------------
  subroutine base_functional__dels__(this, name, that)
    type(base_functional_t), intent(inout) :: this
    character(len=*),        intent(in)    :: name
    type(base_functional_t), intent(in)    :: that

    PUSH_SUB(base_functional__dels__)

    ASSERT(associated(this%config))
    ASSERT(len_trim(name)>0)
    ASSERT(associated(that%config))
    call msgbus_detach(this%msgb, that%msgb)

    POP_SUB(base_functional__dels__)
  end subroutine base_functional__dels__

  ! ---------------------------------------------------------
  subroutine base_functional_set_info(this, energy)
    type(base_functional_t), intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: energy

    PUSH_SUB(base_functional_set_info)

    if(present(energy)) this%energy = energy

    POP_SUB(base_functional_set_info)
  end subroutine base_functional_set_info

  ! ---------------------------------------------------------
  subroutine base_functional_get_info(this, id, family, kind, size, nspin, polarized, use, energy)
    type(base_functional_t), intent(in)  :: this
    integer,       optional, intent(out) :: id
    integer,       optional, intent(out) :: family
    integer,       optional, intent(out) :: kind
    integer,       optional, intent(out) :: size
    integer,       optional, intent(out) :: nspin
    logical,       optional, intent(out) :: polarized
    logical,       optional, intent(out) :: use
    real(kind=wp), optional, intent(out) :: energy

    PUSH_SUB(base_functional_get_info)

    ASSERT(associated(this%config))
    if(present(nspin)) nspin = this%nspin
    call functional_get(this%funct, id=id, family=family, kind=kind, polarized=polarized)
    call storage_get(this%data, size=size, alloc=use)
    if(present(energy)) energy = this%energy

    POP_SUB(base_functional_get_info)
  end subroutine base_functional_get_info

  ! ---------------------------------------------------------
  subroutine base_functional_get_config(this, that)
    type(base_functional_t), target, intent(in)  :: this
    type(json_object_t),    pointer, intent(out) :: that

    PUSH_SUB(base_functional_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_functional_get_config)
  end subroutine base_functional_get_config

  ! ---------------------------------------------------------
  subroutine base_functional_get_simulation(this, that)
    type(base_functional_t), target, intent(in)  :: this
    type(simulation_t),     pointer, intent(out) :: that

    PUSH_SUB(base_functional_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_functional_get_simulation)
  end subroutine base_functional_get_simulation

  ! ---------------------------------------------------------
  subroutine base_functional_get_system(this, that)
    type(base_functional_t), target, intent(in)  :: this
    type(base_system_t),    pointer, intent(out) :: that

    PUSH_SUB(base_functional_get_system)

    nullify(that)
    if(associated(this%sys)) that => this%sys

    POP_SUB(base_functional_get_system)
  end subroutine base_functional_get_system

  ! ---------------------------------------------------------
  subroutine base_functional_get_density(this, that)
    type(base_functional_t),         intent(in)  :: this
    type(base_density_t),   pointer, intent(out) :: that

    PUSH_SUB(base_functional_get_density)

    nullify(that)
    if(associated(this%sys))&
      call base_system_get(this%sys, that)

    POP_SUB(base_functional_get_density)
  end subroutine base_functional_get_density

  ! ---------------------------------------------------------
  subroutine base_functional_get_storage(this, that)
    type(base_functional_t), target, intent(in)  :: this
    type(storage_t),        pointer, intent(out) :: that

    PUSH_SUB(base_functional_get_storage)

    nullify(that)
    if(associated(this%config)) that => this%data

    POP_SUB(base_functional_get_storage)
  end subroutine base_functional_get_storage

  ! ---------------------------------------------------------
  subroutine base_functional_get_functional_1d(this, that)
    type(base_functional_t),              intent(inout) :: this
    real(kind=wp), dimension(:), pointer, intent(out)   :: that

    logical :: allc

    PUSH_SUB(base_functional_get_functional_1d)

    nullify(that)
    call storage_get(this%data, alloc=allc)
    if(allc)then
      call base_functional_update(this)
      call storage_get(this%data, that)
    end if

    POP_SUB(base_functional_get_functional_1d)
  end subroutine base_functional_get_functional_1d

  ! ---------------------------------------------------------
  subroutine base_functional_get_functional_md(this, that)
    type(base_functional_t),                intent(inout) :: this
    real(kind=wp), dimension(:,:), pointer, intent(out)   :: that

    logical :: allc

    PUSH_SUB(base_functional_get_functional_md)

    nullify(that)
    call storage_get(this%data, alloc=allc)
    if(allc)then
      call base_functional_update(this)
      call storage_get(this%data, that)
    end if

    POP_SUB(base_functional_get_functional_md)
  end subroutine base_functional_get_functional_md

  ! ---------------------------------------------------------
  subroutine base_functional__copy__(this, that)
    type(base_functional_t), intent(inout) :: this
    type(base_functional_t), intent(in)    :: that

    type(refcount_t), pointer :: rcnt

    PUSH_SUB(base_functional__copy__)

    rcnt => this%rcnt
    nullify(this%rcnt)
    call base_functional__end__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_functional__init__(this, that)
      call refcount_del(this%rcnt)
      this%energy = that%energy
      if(associated(that%sim)) call storage_copy(this%data, that%data)
    end if
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(base_functional__copy__)
  end subroutine base_functional__copy__

  ! ---------------------------------------------------------
  subroutine base_functional__end__(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional__end__)

    nullify(this%config, this%sys, this%sim)
    if(associated(this%rcnt)) call refcount_del(this%rcnt)
    this%nspin = default_nspin
    this%factor = 1.0_wp
    this%energy = 0.0_wp
    call functional_end(this%funct)
    call storage_end(this%data)
    call msgbus_end(this%msgb)
    ASSERT(base_functional_dict_len(this%dict)==0)
    call base_functional_dict_end(this%dict)

    POP_SUB(base_functional__end__)
  end subroutine base_functional__end__

end module base_functional_oct_m

!! Local Variables:
!! mode: f90
!! End:
