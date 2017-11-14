#include "global.h"

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

module base_handle_oct_m

  use base_density_oct_m
  use base_geometry_oct_m
  use base_hamiltonian_oct_m
  use base_model_oct_m
  use base_states_oct_m
  use base_system_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use space_oct_m

#define BASE_TEMPLATE_NAME base_handle
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private

  public ::         &
    HNDL_TYPE_NONE

  public ::        &
    base_handle_t

  public ::                &
    base_handle__init__,   &
    base_handle__start__,  &
    base_handle__update__, &
    base_handle__reset__,  &
    base_handle__stop__,   &
    base_handle__copy__,   &
    base_handle__end__

  public ::             &
    base_handle_new,    &
    base_handle_del,    &
    base_handle_init,   &
    base_handle_start,  &
    base_handle_update, &
    base_handle_reset,  &
    base_handle_stop,   &
    base_handle_get,    &
    base_handle_copy,   &
    base_handle_end

#define BASE_TEMPLATE_NAME base_handle
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  integer, parameter :: HNDL_TYPE_NONE = 0

  type :: base_handle_t
    private
    type(json_object_t), pointer :: config =>null()
    type(base_handle_t), pointer :: prnt   =>null()
    integer                      :: type   = HNDL_TYPE_NONE
    type(simulation_t)           :: sim
    type(base_model_t)           :: model
    type(base_handle_dict_t)     :: dict
    type(base_handle_list_t)     :: list
  end type base_handle_t

  interface base_handle__init__
    module procedure base_handle__init__type
    module procedure base_handle__init__pass
    module procedure base_handle__init__smlt
    module procedure base_handle__init__copy
  end interface base_handle__init__

  interface base_handle_init
    module procedure base_handle_init_type
    module procedure base_handle_init_pass
    module procedure base_handle_init_copy
  end interface base_handle_init

  interface base_handle_get
    module procedure base_handle_get_info
    module procedure base_handle_get_config
    module procedure base_handle_get_simulation
    module procedure base_handle_get_geom
    module procedure base_handle_get_density
    module procedure base_handle_get_states
    module procedure base_handle_get_system
    module procedure base_handle_get_hamiltonian
    module procedure base_handle_get_model
  end interface base_handle_get

  interface base_handle_copy
    module procedure base_handle_copy_type
  end interface base_handle_copy

  interface base_handle_end
    module procedure base_handle_end_type
    module procedure base_handle_end_pass
  end interface base_handle_end

contains

#define BASE_TEMPLATE_NAME base_handle
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  subroutine base_handle__new__(this)
    type(base_handle_t), pointer :: this

    PUSH_SUB(base_handle__new__)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(base_handle__new__)
  end subroutine base_handle__new__

  ! ---------------------------------------------------------
  subroutine base_handle__del__(this)
    type(base_handle_t), pointer :: this

    PUSH_SUB(base_handle__del__)

    if(associated(this))then
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(base_handle__del__)
  end subroutine base_handle__del__

  ! ---------------------------------------------------------
  subroutine base_handle_new(this, that)
    type(base_handle_t),  target, intent(inout) :: this
    type(base_handle_t), pointer                :: that

    PUSH_SUB(base_handle_new)

    nullify(that)
    call base_handle__new__(that)
    that%prnt => this
    call base_handle_list_push(this%list, that)

    POP_SUB(base_handle_new)
  end subroutine base_handle_new

  ! ---------------------------------------------------------
  subroutine base_handle_del(this)
    type(base_handle_t), pointer :: this

    PUSH_SUB(base_handle_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_handle_list_del(this%prnt%list, this)
        call base_handle_end(this)
        call base_handle__del__(this)
      end if
    end if

    POP_SUB(base_handle_del)
  end subroutine base_handle_del

  ! ---------------------------------------------------------
  subroutine base_handle__init__type(this, config)
    type(base_handle_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_handle__init__type)

    nullify(cnfg)
    this%config => config
    call json_get(this%config, "type", this%type, ierr)
    if(ierr/=JSON_OK) this%type = HNDL_TYPE_NONE
    call json_get(this%config, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_model__init__(this%model, cnfg)
    nullify(cnfg)
    call base_handle_dict_init(this%dict)
    call base_handle_list_init(this%list)
    
    POP_SUB(base_handle__init__type)
  end subroutine base_handle__init__type

  ! ---------------------------------------------------------
  recursive subroutine base_handle__init__pass(this, init)
    type(base_handle_t), target, intent(inout) :: this

    interface
      subroutine init(this, config)
        use json_oct_m
        import :: base_handle_t
        type(base_handle_t), intent(out) :: this
        type(json_object_t), intent(in)  :: config
      end subroutine init
    end interface

    type(json_object_iterator_t)        :: iter
    character(len=BASE_HANDLE_NAME_LEN) :: name
    type(json_object_t),        pointer :: dict, cnfg
    type(base_handle_t),        pointer :: hndl
    integer                             :: ierr

    PUSH_SUB(base_handle__init__pass)

    ASSERT(associated(this%config))
    nullify(dict, cnfg, hndl)
    call json_get(this%config, "subsystems", dict, ierr)
    if(ierr==JSON_OK)then
      call json_init(iter, dict)
      do
        nullify(cnfg, hndl)
        call json_next(iter, name, cnfg, ierr)
        if(ierr/=JSON_OK)exit
        call base_handle_new(this, hndl)
        ASSERT(associated(hndl))
        call init(hndl, cnfg)
        hndl%prnt => this
        call base_handle_sets(this, trim(adjustl(name)), hndl)
      end do
      call json_end(iter)
      nullify(cnfg, hndl)
    end if
    nullify(dict)

    POP_SUB(base_handle__init__pass)
  end subroutine base_handle__init__pass

  ! ---------------------------------------------------------
  subroutine base_handle__init__smlt(this)
    type(base_handle_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(geometry_t),    pointer :: geo
    type(space_t),       pointer :: space
    integer                      :: ierr

    PUSH_SUB(base_handle__init__smlt)

    nullify(cnfg, geo, space)
    call base_model_get(this%model, geo)
    ASSERT(associated(geo))
    call base_model_get(this%model, space)
    ASSERT(associated(space))
    call json_get(this%config, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call simulation_init(this%sim, geo, space, cnfg)
    nullify(cnfg, geo, space)

    POP_SUB(base_handle__init__smlt)
  end subroutine base_handle__init__smlt

  ! ---------------------------------------------------------
  subroutine base_handle__init__copy(this, that)
    type(base_handle_t), intent(out) :: this
    type(base_handle_t), intent(in)  :: that

    PUSH_SUB(base_handle__init__copy)

    ASSERT(associated(that%config))
    call base_handle__init__(this, that%config)
    if(simulation_assoc(that%sim))then
      call simulation_copy(this%sim, that%sim)
      call base_model__start__(this%model, this%sim)
    end if

    POP_SUB(base_handle__init__copy)
  end subroutine base_handle__init__copy

  ! ---------------------------------------------------------
  recursive subroutine base_handle_init_type(this, config)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_handle_init_type)

    call base_handle_init(this, config, base_handle_init_type)

    POP_SUB(base_handle_init_type)
  end subroutine base_handle_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_handle_init_pass(this, config, init)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    interface
      subroutine init(this, config)
        use json_oct_m
        import :: base_handle_t
        type(base_handle_t), intent(out) :: this
        type(json_object_t), intent(in)  :: config
      end subroutine init
    end interface

    PUSH_SUB(base_handle_init_pass)

    call base_handle__init__(this, config)
    call base_handle__init__(this, init)
    call base_handle__init__(this)

    POP_SUB(base_handle_init_pass)
  end subroutine base_handle_init_pass

  ! ---------------------------------------------------------
  recursive subroutine base_handle_init_copy(this, that)
    type(base_handle_t), intent(out) :: this
    type(base_handle_t), intent(in)  :: that

    type(base_handle_iterator_t)        :: iter
    character(len=BASE_HANDLE_NAME_LEN) :: name
    type(base_handle_t),        pointer :: osub, isub
    integer                             :: ierr

    PUSH_SUB(base_handle_init_copy)

    call base_handle__init__(this, that)
    call base_handle_init(iter, that)
    do
      nullify(osub, isub)
      call base_handle_next(iter, name, isub, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_new(this, osub)
      call base_handle_init(osub, isub)
      call base_handle_sets(this, name, osub)
    end do
    call base_handle_end(iter)
    nullify(osub, isub)
    call base_handle__init__(this)

    POP_SUB(base_handle_init_copy)
  end subroutine base_handle_init_copy

  ! ---------------------------------------------------------
  subroutine base_handle__start__(this, grid)
    type(base_handle_t),    intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid

    PUSH_SUB(base_handle__start__)

    ASSERT(associated(this%config))
    !call base_handle__reduce__(this, extend)
    if(present(grid)) call simulation_start(this%sim, grid)
    ASSERT(simulation_assoc(this%sim))
    call base_model__start__(this%model, this%sim)

    POP_SUB(base_handle__start__)
    
  contains

    subroutine extend(this, that)
      type(base_handle_t), intent(inout) :: this
      type(base_handle_t), intent(in)    :: that

      PUSH_SUB(base_handle__start__.extend)
      
      call simulation_extend(this%sim, that%sim, that%config)

      POP_SUB(base_handle__start__.extend)
    end subroutine extend

  end subroutine base_handle__start__

  ! ---------------------------------------------------------
  subroutine base_handle__update__(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle__update__)

    ASSERT(associated(this%config))
    call base_model__update__(this%model)

    POP_SUB(base_handle__update__)
  end subroutine base_handle__update__

  ! ---------------------------------------------------------
  subroutine base_handle__reset__(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle__reset__)

    ASSERT(associated(this%config))
    call base_model__reset__(this%model)

    POP_SUB(base_handle__reset__)
  end subroutine base_handle__reset__

  ! ---------------------------------------------------------
  subroutine base_handle__stop__(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle__stop__)

    ASSERT(associated(this%config))
    call base_model__stop__(this%model)

    POP_SUB(base_handle__stop__)
  end subroutine base_handle__stop__

  ! ---------------------------------------------------------
  subroutine base_handle_start(this, grid)
    type(base_handle_t),    intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid

    PUSH_SUB(base_handle_start)

    call base_handle__apply__(this, start)

    POP_SUB(base_handle_start)
    
  contains

    subroutine start(this)
      type(base_handle_t), intent(inout) :: this

      PUSH_SUB(base_handle_start.start)
      
      call base_handle__start__(this, grid)

      POP_SUB(base_handle_start.start)
    end subroutine start

  end subroutine base_handle_start

  ! ---------------------------------------------------------
  subroutine base_handle_update(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle_update)

    call base_handle__apply__(this, base_handle__update__)

    POP_SUB(base_handle_update)
  end subroutine base_handle_update

  ! ---------------------------------------------------------
  subroutine base_handle_reset(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle_reset)

    call base_handle__apply__(this, base_handle__reset__)

    POP_SUB(base_handle_reset)
  end subroutine base_handle_reset

  ! ---------------------------------------------------------
  subroutine base_handle_stop(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle_stop)

    call base_handle__apply__(this, base_handle__stop__)

    POP_SUB(base_handle_stop)
  end subroutine base_handle_stop

  ! ---------------------------------------------------------
  subroutine base_handle__sets__(this, name, that)
    type(base_handle_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_handle_t), intent(in)    :: that

    PUSH_SUB(base_handle__sets__)

    call base_model_sets(this%model, trim(adjustl(name)), that%model)

    POP_SUB(base_handle__sets__)
  end subroutine base_handle__sets__

  ! ---------------------------------------------------------
  subroutine base_handle__dels__(this, name, ierr)
    type(base_handle_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    integer,             intent(out)   :: ierr

    PUSH_SUB(base_handle__dels__)

    call base_model_dels(this%model, trim(adjustl(name)), ierr)

    POP_SUB(base_handle__dels__)
  end subroutine base_handle__dels__

  ! ---------------------------------------------------------
  subroutine base_handle_get_info(this, type)
    type(base_handle_t), intent(in)  :: this
    integer,  optional,  intent(out) :: type

    PUSH_SUB(base_handle_get_info)

    if(present(type)) type = this%type

    POP_SUB(base_handle_get_info)
  end subroutine base_handle_get_info

  ! ---------------------------------------------------------
  subroutine base_handle_get_config(this, that)
    type(base_handle_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_handle_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_handle_get_config)
  end subroutine base_handle_get_config

  ! ---------------------------------------------------------
  subroutine base_handle_get_simulation(this, that)
    type(base_handle_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(base_handle_get_simulation)

    that => this%sim

    POP_SUB(base_handle_get_simulation)
  end subroutine base_handle_get_simulation

  ! ---------------------------------------------------------
  subroutine base_handle_get_geom(this, that)
    type(base_handle_t),    intent(in) :: this
    type(base_geometry_t), pointer     :: that

    PUSH_SUB(base_handle_get_geom)

    call base_model_get(this%model, that)
    
    POP_SUB(base_handle_get_geom)
  end subroutine base_handle_get_geom

  ! ---------------------------------------------------------
  subroutine base_handle_get_density(this, that)
    type(base_handle_t),   intent(in) :: this
    type(base_density_t), pointer     :: that

    PUSH_SUB(base_handle_get_density)

    call base_model_get(this%model, that)
    
    POP_SUB(base_handle_get_density)
  end subroutine base_handle_get_density

  ! ---------------------------------------------------------
  subroutine base_handle_get_states(this, that)
    type(base_handle_t),  intent(in) :: this
    type(base_states_t), pointer     :: that

    PUSH_SUB(base_handle_get_states)

    call base_model_get(this%model, that)
    
    POP_SUB(base_handle_get_states)
  end subroutine base_handle_get_states

  ! ---------------------------------------------------------
  subroutine base_handle_get_system(this, that)
    type(base_handle_t),  intent(in) :: this
    type(base_system_t), pointer     :: that

    PUSH_SUB(base_handle_get_system)

    call base_model_get(this%model, that)
    
    POP_SUB(base_handle_get_system)
  end subroutine base_handle_get_system

  ! ---------------------------------------------------------
  subroutine base_handle_get_hamiltonian(this, that)
    type(base_handle_t),       intent(in) :: this
    type(base_hamiltonian_t), pointer     :: that

    PUSH_SUB(base_handle_get_hamiltonian)

    call base_model_get(this%model, that)
    
    POP_SUB(base_handle_get_hamiltonian)
  end subroutine base_handle_get_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_handle_get_model(this, that)
    type(base_handle_t), target, intent(in) :: this
    type(base_model_t), pointer             :: that

    PUSH_SUB(base_handle_get_model)

    nullify(that)
    if(associated(this%config)) that => this%model

    POP_SUB(base_handle_get_model)
  end subroutine base_handle_get_model

  ! ---------------------------------------------------------
  subroutine base_handle__copy__(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    PUSH_SUB(base_handle__copy__)

    call base_handle__end__(this)
    if(associated(that%config))then
      call base_handle__init__(this, that)
      if(simulation_assoc(that%sim)) call base_model__copy__(this%model, that%model)
    end if

    POP_SUB(base_handle__copy__)
  end subroutine base_handle__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_copy_type(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    type(base_handle_iterator_t)        :: iter
    character(len=BASE_HANDLE_NAME_LEN) :: name
    type(base_handle_t),        pointer :: osub, isub
    integer                             :: ierr

    PUSH_SUB(base_handle_copy_type)

    nullify(osub, isub)
    call base_handle_end(this)
    call base_handle__copy__(this, that)
    call base_handle_init(iter, that)
    do
      nullify(osub, isub)
      call base_handle_next(iter, name, isub, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_new(this, osub)
      call base_handle_copy(osub, isub)
      call base_handle_sets(this, name, osub)
    end do
    call base_handle_end(iter)
    nullify(osub, isub)

    POP_SUB(base_handle_copy_type)
  end subroutine base_handle_copy_type

  ! ---------------------------------------------------------
  subroutine base_handle__end__(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle__end__)

    nullify(this%config, this%prnt)
    this%type = HNDL_TYPE_NONE
    call simulation_end(this%sim)
    call base_model__end__(this%model)
    call base_handle_dict_end(this%dict)
    ASSERT(base_handle_list_len(this%list)==0)
    call base_handle_list_end(this%list)

    POP_SUB(base_handle__end__)
  end subroutine base_handle__end__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_end_type(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle_end_type)

    call base_handle_end(this, base_handle_end_type)
    
    POP_SUB(base_handle_end_type)
  end subroutine base_handle_end_type

  ! ---------------------------------------------------------
  recursive subroutine base_handle_end_pass(this, finis)
    type(base_handle_t), intent(inout) :: this

    interface
      subroutine finis(this)
        import :: base_handle_t
        type(base_handle_t), intent(inout) :: this
      end subroutine finis
    end interface
    
    type(base_handle_t), pointer :: subs

    PUSH_SUB(base_handle_end_pass)

    do
      nullify(subs)
      call base_handle_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call finis(subs)
      call base_handle__del__(subs)
    end do
    nullify(subs)
    call base_handle__end__(this)

    POP_SUB(base_handle_end_pass)
  end subroutine base_handle_end_pass

end module base_handle_oct_m

!! Local Variables:
!! mode: f90
!! End:
