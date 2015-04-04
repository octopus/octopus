#include "global.h"

module base_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: json_object_t

  use simulation_m, only: &
    simulation_t

  use base_system_m, only: &
    base_system_t

  use root_hamiltonian_m, only: &
    root_hamiltonian__rpush__,  &
    root_hamiltonian__rpop__,   &
    root_hamiltonian__rdel__,   &
    root_hamiltonian__rset__,   &
    root_hamiltonian__rget__

  use root_hamiltonian_m, only: &
    root_hamiltonian__init__,   &
    root_hamiltonian__start__,  &
    root_hamiltonian__update__, &
    root_hamiltonian__stop__,   &
    root_hamiltonian__reset__,  &
    root_hamiltonian__acc__,    &
    root_hamiltonian__add__,    &
    root_hamiltonian__copy__,   &
    root_hamiltonian__end__

  use root_hamiltonian_m, only:                   &
    base_hamiltonian_t   => root_hamiltonian_t,   &
    base_hamiltonian_set => root_hamiltonian_set, &
    base_hamiltonian_get => root_hamiltonian_get

  use root_hamiltonian_m, only: &
    HMLT_TYPE_NONE,             &
    HMLT_TYPE_TERM,             &
    HMLT_TYPE_POTN,             &
    HMLT_TYPE_FNCT,             &
    HMLT_TYPE_HMLT

  use root_hamiltonian_m, only:                                   &
    BASE_HAMILTONIAN_OK          => ROOT_HAMILTONIAN_OK,          &
    BASE_HAMILTONIAN_KEY_ERROR   => ROOT_HAMILTONIAN_KEY_ERROR,   &
    BASE_HAMILTONIAN_EMPTY_ERROR => ROOT_HAMILTONIAN_EMPTY_ERROR

#define TEMPLATE_NAME hamiltonian
#define INTERNAL_ACCESS
#define INCLUDE_PREFIX
#include "iterator_code.F90"
#undef INCLUDE_PREFIX
#undef INTERNAL_ACCESS
#undef TEMPLATE_NAME

  implicit none

  private

  public ::                  &
    base_hamiltonian_t,      &
    base_hamiltonian_new,    &
    base_hamiltonian_del,    &
    base_hamiltonian_init,   &
    base_hamiltonian_start,  &
    base_hamiltonian_update, &
    base_hamiltonian_stop,   &
    base_hamiltonian_next,   &
    base_hamiltonian_set,    &
    base_hamiltonian_get,    &
    base_hamiltonian_copy,   &
    base_hamiltonian_end

  public ::         &
    HMLT_TYPE_NONE, &
    HMLT_TYPE_TERM, &
    HMLT_TYPE_POTN, &
    HMLT_TYPE_FNCT, &
    HMLT_TYPE_HMLT

  public ::                       &
    BASE_HAMILTONIAN_OK,          &
    BASE_HAMILTONIAN_KEY_ERROR,   &
    BASE_HAMILTONIAN_EMPTY_ERROR

  interface base_hamiltonian_init
    module procedure base_hamiltonian_init_hamiltonian
    module procedure base_hamiltonian_init_copy
  end interface base_hamiltonian_init

  interface base_hamiltonian_update
    module procedure base_hamiltonian_update_hamiltonian
    module procedure base_hamiltonian_update_pass
  end interface base_hamiltonian_update

  interface base_hamiltonian_copy
    module procedure base_hamiltonian_copy_hamiltonian
  end interface base_hamiltonian_copy

  interface base_hamiltonian_end
    module procedure base_hamiltonian_end_hamiltonian
  end interface base_hamiltonian_end

#define TEMPLATE_NAME hamiltonian
#define INTERNAL_ACCESS
#define INCLUDE_HEADER
#include "iterator_code.F90"
#undef INCLUDE_HEADER
#undef INTERNAL_ACCESS
#undef TEMPLATE_NAME

contains

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_new(this, that)
    type(base_hamiltonian_t),  target, intent(inout) :: this
    type(base_hamiltonian_t), pointer                :: that
    !
    PUSH_SUB(base_hamiltonian_new)
    nullify(that)
    SAFE_ALLOCATE(that)
    call root_hamiltonian__rset__(that, this)
    call root_hamiltonian__rpush__(this, that)
    POP_SUB(base_hamiltonian_new)
    return
  end subroutine base_hamiltonian_new

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__idel__(this)
    type(base_hamiltonian_t), pointer :: this
    !
    PUSH_SUB(base_hamiltonian__idel__)
    SAFE_DEALLOCATE_P(this)
    nullify(this)
    POP_SUB(base_hamiltonian__idel__)
    return
  end subroutine base_hamiltonian__idel__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del(this)
    type(base_hamiltonian_t), pointer :: this
    !
    type(base_hamiltonian_t), pointer :: prnt
    !
    PUSH_SUB(base_hamiltonian_del)
    nullify(prnt)
    if(associated(this))then
      call root_hamiltonian__rget__(this, prnt)
      if(associated(prnt))then
        call root_hamiltonian__rdel__(prnt, this)
        call base_hamiltonian_end(this)
        call base_hamiltonian__idel__(this)
      end if
      nullify(prnt)
    end if
    POP_SUB(base_hamiltonian_del)
    return
  end subroutine base_hamiltonian_del

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_init_hamiltonian(this, sys, config)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_system_t),      intent(in)  :: sys
    type(json_object_t),      intent(in)  :: config
    !
    PUSH_SUB(base_hamiltonian_init_hamiltonian)
    call root_hamiltonian__init__(this, sys, config)
    POP_SUB(base_hamiltonian_init_hamiltonian)
    return
  end subroutine base_hamiltonian_init_hamiltonian

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_init_copy(this, that)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_hamiltonian_t), intent(in)  :: that
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: osub, isub
    type(json_object_t),      pointer :: cnfg
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_init_copy)
    nullify(cnfg, osub, isub)
    call root_hamiltonian__init__(this, that)
    call base_hamiltonian_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_hamiltonian_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_new(this, osub)
      call base_hamiltonian_init(osub, isub)
      call root_hamiltonian__add__(this, osub, cnfg)
    end do
    call base_hamiltonian_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(base_hamiltonian_init_copy)
    return
  end subroutine base_hamiltonian_init_copy

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_start(this, sim)
    type(base_hamiltonian_t), intent(inout) :: this
    type(simulation_t),       intent(in)    :: sim
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: subs
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_start)
    nullify(subs)
    call base_hamiltonian_init(iter, this)
    do
      nullify(subs)
      call base_hamiltonian_next(iter, subs, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_start(subs, sim)
    end do
    call base_hamiltonian_end(iter)
    nullify(subs)
    call root_hamiltonian__start__(this, sim)
    POP_SUB(base_hamiltonian_start)
    return
  end subroutine base_hamiltonian_start

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_update_hamiltonian(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    PUSH_SUB(base_hamiltonian_update_hamiltonian)
    call base_hamiltonian_update(this, root_hamiltonian__update__)
    POP_SUB(base_hamiltonian_update_hamiltonian)
    return
  end subroutine base_hamiltonian_update_hamiltonian

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_update_pass(this, hamiltonian_update)
    type(base_hamiltonian_t), intent(inout) :: this
    interface
      subroutine hamiltonian_update(this)
        import :: base_hamiltonian_t
        type(base_hamiltonian_t), intent(inout) :: this
      end subroutine hamiltonian_update
    end interface
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: subs
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_update_pass)
    nullify(subs)
    call root_hamiltonian__reset__(this)
    call hamiltonian_update(this)
    call base_hamiltonian_init(iter, this)
    do
      nullify(subs)
      call base_hamiltonian_next(iter, subs, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_update(subs, hamiltonian_update)
      call root_hamiltonian__acc__(this, subs)
    end do
    call base_hamiltonian_end(iter)
    nullify(subs)
    POP_SUB(base_hamiltonian_update_pass)
    return
  end subroutine base_hamiltonian_update_pass

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_stop(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: subs
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_stop)
    nullify(subs)
    call base_hamiltonian_init(iter, this)
    do
      nullify(subs)
      call base_hamiltonian_next(iter, subs, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_stop(subs)
    end do
    call base_hamiltonian_end(iter)
    nullify(subs)
    call root_hamiltonian__stop__(this)
    POP_SUB(base_hamiltonian_stop)
    return
  end subroutine base_hamiltonian_stop

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_copy_hamiltonian(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: osub, isub
    type(json_object_t),      pointer :: cnfg
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_copy_hamiltonian)
    nullify(cnfg, osub, isub)
    call base_hamiltonian_end(this)
    call root_hamiltonian__copy__(this, that)
    call base_hamiltonian_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_hamiltonian_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_new(this, osub)
      call base_hamiltonian_copy(osub, isub)
      call root_hamiltonian__add__(this, osub, cnfg)
    end do
    call base_hamiltonian_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(base_hamiltonian_copy_hamiltonian)
    return
  end subroutine base_hamiltonian_copy_hamiltonian

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_end_hamiltonian(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(base_hamiltonian_t), pointer :: subs
    !
    PUSH_SUB(base_hamiltonian_end_hamiltonian)
    do
      nullify(subs)
      call root_hamiltonian__rpop__(this, subs)
      if(.not.associated(subs))exit
      call base_hamiltonian_end(subs)
      call base_hamiltonian__idel__(subs)
    end do
    nullify(subs)
    call root_hamiltonian__end__(this)
    POP_SUB(base_hamiltonian_end_hamiltonian)
    return
  end subroutine base_hamiltonian_end_hamiltonian

#define TEMPLATE_NAME hamiltonian
#define INTERNAL_ACCESS
#define INCLUDE_BODY
#include "iterator_code.F90"
#undef INCLUDE_BODY
#undef INTERNAL_ACCESS
#undef TEMPLATE_NAME

end module base_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
