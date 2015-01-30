#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME
#undef LIST_INCLUDE_PREFIX
#undef LIST_INCLUDE_HEADER
#undef LIST_INCLUDE_BODY

#define LIST_TEMPLATE_NAME domain
#define LIST_INCLUDE_PREFIX
#include "tlist.F90"
#undef LIST_INCLUDE_PREFIX

module domain_m

  use global_m
  use messages_m
  use profiling_m

  use kinds_m,     only: wp
  use geometry_m,  only: geometry_t
  use simul_box_m, only: simul_box_t, simul_box_in_box

  implicit none

  private
  public ::           &
    domain_init,      &
    domain_start,     &
    domain_in_domain, &
    domain_copy,      &
    domain_end

#define LIST_INCLUDE_HEADER
#include "tlist.F90"
#undef LIST_INCLUDE_HEADER

  type, public :: domain_t
    private
    type(simul_box_t), pointer :: sb  =>null()
    type(geometry_t),  pointer :: geo =>null()
    type(domain_list_t)        :: list
  end type domain_t

  interface domain_init
    module procedure domain_init_domain
    module procedure domain_init_build
  end interface domain_init

contains
  
#define LIST_INCLUDE_BODY
#include "tlist.F90"
#undef LIST_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine domain_init_domain(this)
    type(domain_t), intent(out) :: this
    !
    PUSH_SUB(domain_init_domain)
    call domain_list_init(this%list)
    POP_SUB(domain_init_domain)
    return
  end subroutine domain_init_domain

  ! ---------------------------------------------------------
  subroutine domain_init_build(this, that)
    type(domain_t), intent(inout) :: this
    type(domain_t), intent(in)    :: that
    !
    PUSH_SUB(domain_init_build)
    call domain_list_push(this%list, that)
    POP_SUB(domain_init_build)
    return
  end subroutine domain_init_build

  ! ---------------------------------------------------------
  subroutine domain_start(this, sb, geo)
    type(domain_t),            intent(out) :: this
    type(simul_box_t), target, intent(in)  :: sb
    type(geometry_t),  target, intent(in)  :: geo
    !
    PUSH_SUB(domain_start)
    ASSERT(.not.associated(this%sb))
    ASSERT(.not.associated(this%geo))
    this%sb=>sb
    this%geo=>geo
    POP_SUB(domain_start)
    return
  end subroutine domain_start

  ! ---------------------------------------------------------
  function domain_in_domain_aux(this, x) result(in)
    type(domain_t),              intent(in) :: this
    real(kind=wp), dimension(:), intent(in) :: x
    !
    logical :: in
    !
    PUSH_SUB(domain_in_domain_aux)
    ASSERT(associated(this%sb))
    ASSERT(associated(this%geo))
    in=simul_box_in_box(this%sb, this%geo, x)
    POP_SUB(domain_in_domain_aux)
    return
  end function domain_in_domain_aux

  ! ---------------------------------------------------------
  function domain_in_domain(this, x) result(in)
    type(domain_t),              intent(in) :: this
    real(kind=wp), dimension(:), intent(in) :: x
    !
    logical :: in
    !
    type(domain_list_iterator_t) :: iter
    type(domain_t),      pointer :: domain
    integer                      :: ierr
    !
    PUSH_SUB(domain_in_domain)
    nullify(domain)
    in=domain_in_domain_aux(this, x)
    if(.not.in)then
      call domain_list_init(iter, this%list)
      do
        nullify(domain)
        call domain_list_next(iter, domain, ierr)
        if(ierr/=DOMAIN_LIST_OK)exit
        in=domain_in_domain_aux(domain, x)
        if(in)exit
      end do
      call domain_list_end(iter)
    end if
    POP_SUB(domain_in_domain)
    return
  end function domain_in_domain

  ! ---------------------------------------------------------
  subroutine domain_copy(this, that)
    type(domain_t), intent(out) :: this
    type(domain_t), intent(in)  :: that
    !
    PUSH_SUB(domain_copy)
    this%sb=>that%sb
    this%geo=>that%geo
    call domain_list_copy(this%list, that%list)
    POP_SUB(domain_copy)
    return
  end subroutine domain_copy

  ! ---------------------------------------------------------
  subroutine domain_end(this)
    type(domain_t), intent(inout) :: this
    !
    nullify(this%sb, this%geo)
    call domain_list_end(this%list)
    return
  end subroutine domain_end

end module domain_m

#undef LIST_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
