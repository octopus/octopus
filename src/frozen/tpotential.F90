#include "global.h"
#include "template.h"

module TEMPLATE(potential_m)

  use global_m
  use messages_m
  use profiling_m

  use json_m,          only: JSON_OK, json_object_t, json_get
  use kinds_m,         only: wp

  use storage_m, only:     &
    storage_t,             &
    storage_init,          &
    storage_update,        &
    storage_get_size,      &
    storage_get_dimension, &
    storage_get_storage,   &
    storage_copy,          &
    storage_end

  use storage_m, only:          &
    storage_interpolation_t,    &
    storage_interpolation_init, &
    storage_interpolation_eval, &
    storage_interpolation_copy, &
    storage_interpolation_end

  use storage_m, only:                              &
    TEMPLATE(POTENTIAL_INTRP_OK)=>STORAGE_INTRP_OK, &
    TEMPLATE(POTENTIAL_INTRP_OD)=>STORAGE_INTRP_OD, &
    TEMPLATE(POTENTIAL_INTRP_NI)=>STORAGE_INTRP_NI

  use TEMPLATE(density_m), only:      &
    density_t => TEMPLATE(density_t)

  use TEMPLATE(geometry_m), only:       &
    geometry_t => TEMPLATE(geometry_t)

  use TEMPLATE(simulation_m), only:         &
    simulation_t => TEMPLATE(simulation_t)

  use TEMPLATE(system_m), only:         &
    system_t   => TEMPLATE(system_t),   &
    system_get => TEMPLATE(system_get)

  use TEMPLATE(term_m), only: &
    term_t          => TEMPLATE(term_t),          &
    term_init       => TEMPLATE(term_init),       &
#ifdef SUBTEMPLATE_NAME
    term_extend     => TEMPLATE(term_extend),     &
#endif
    term_get        => TEMPLATE(term_get),        &
    term_set_energy => TEMPLATE(term_set_energy), &
    term_get_energy => TEMPLATE(term_get_energy), &
    term_copy       => TEMPLATE(term_copy),       &
    term_end        => TEMPLATE(term_end)

#ifdef SUBTEMPLATE_NAME
  use SUBTEMPLATE(m), only:  &
    sub_t => SUBTEMPLATE(t)
#endif

  implicit none

  private
  public ::                            &
    TEMPLATE(potential_init),          &
    TEMPLATE(potential_start),         &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(potential_extend),        &
#endif
    TEMPLATE(potential_update),        &
    TEMPLATE(potential_get_size),      &
    TEMPLATE(potential_get_nspin),     &
    TEMPLATE(potential_get),           &
    TEMPLATE(potential_set_energy),    &
    TEMPLATE(potential_get_energy),    &
    TEMPLATE(potential_get_potential), &
    TEMPLATE(potential_copy),          &
    TEMPLATE(potential_end)

  public ::                                 &
    TEMPLATE(POTENTIAL_INTRP_OK),           & 
    TEMPLATE(POTENTIAL_INTRP_OD),           &
    TEMPLATE(POTENTIAL_INTRP_NI),           &
    TEMPLATE(potential_interpolation_init), & 
    TEMPLATE(potential_interpolation_eval), & 
    TEMPLATE(potential_interpolation_copy), & 
    TEMPLATE(potential_interpolation_end)

  type, public :: TEMPLATE(potential_t)
    private
    type(simulation_t),  pointer :: sim     =>null()
    integer                      :: nspin   = 0
    logical                      :: alloc   = .false.
    real(kind=wp)                :: default = 0.0_wp
    type(term_t)                 :: term
    type(storage_t)              :: potential
  end type TEMPLATE(potential_t)

  type, public :: TEMPLATE(potential_interpolation_t)
    private
    type(TEMPLATE(potential_t)), pointer :: self =>null()
    type(storage_interpolation_t)        :: intrp
  end type TEMPLATE(potential_interpolation_t)

  interface TEMPLATE(potential_get)
    module procedure TEMPLATE(potential_get_config)
    module procedure TEMPLATE(potential_get_simulation)
    module procedure TEMPLATE(potential_get_system)
    module procedure TEMPLATE(potential_get_geometry)
    module procedure TEMPLATE(potential_get_density)
  end interface TEMPLATE(potential_get)

  interface TEMPLATE(potential_get_potential)
    module procedure TEMPLATE(potential_get_potential_1d)
    module procedure TEMPLATE(potential_get_potential_2d)
  end interface TEMPLATE(potential_get_potential)

  interface TEMPLATE(potential_interpolation_eval)
    module procedure TEMPLATE(potential_interpolation_eval_1d)
    module procedure TEMPLATE(potential_interpolation_eval_2d)
  end interface TEMPLATE(potential_interpolation_eval)

contains

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_init)(this, sys, config)
    type(TEMPLATE(potential_t)), intent(out) :: this
    type(system_t),      target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: ierr
    !
    this%sim=>null()
    call term_init(this%term, sys, config)
    call json_get(config, "SpinComponents", this%nspin, ierr)
    if(ierr/=JSON_OK)this%nspin=1
    call json_get(config, "allocate", this%alloc, ierr)
    if(ierr/=JSON_OK)this%alloc=.true.
    this%default=0.0_wp
    return
  end subroutine TEMPLATE(potential_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_start)(this, sim)
    type(TEMPLATE(potential_t)), intent(inout) :: this
    type(simulation_t),  target, intent(in)    :: sim
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call term_get(this%term, cnfg)
    ASSERT(associated(cnfg))
    call json_get(cnfg, "default", this%default, ierr)
    if(ierr/=JSON_OK)this%default=0.0_wp
    if(this%alloc)&
      call storage_init(this%potential, this%sim, this%nspin, this%default)
    return
  end subroutine TEMPLATE(potential_start)

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_extend)(this, that)
    type(TEMPLATE(potential_t)), intent(inout) :: this
    type(sub_t),       optional, intent(in)    :: that
    !
    call term_extend(this%term, that)
    return
  end subroutine TEMPLATE(potential_extend)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_update)(this)
    type(TEMPLATE(potential_t)), intent(inout) :: this
    !
    ASSERT(associated(this%sim))
    if(this%alloc)&
      call storage_update(this%potential)
    return
  end subroutine TEMPLATE(potential_update)

  ! ---------------------------------------------------------
  elemental subroutine TEMPLATE(potential_set_energy)(this, that)
    type(TEMPLATE(potential_t)), intent(inout) :: this
    real(kind=wp),               intent(in)    :: that
    !
    call term_set_energy(this%term, that)
    return
  end subroutine TEMPLATE(potential_set_energy)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(potential_get_size)(this) result(np)
    type(TEMPLATE(potential_t)), intent(in) :: this
    !
    integer :: np
    !
    np=storage_get_size(this%potential)
    return
  end function TEMPLATE(potential_get_size)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(potential_get_nspin)(this) result(that)
    type(TEMPLATE(potential_t)), intent(in) :: this
    !
    integer :: that
    !
    that=storage_get_dimension(this%potential)
    return
  end function TEMPLATE(potential_get_nspin)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(potential_get_energy)(this) result(that)
    type(TEMPLATE(potential_t)), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=term_get_energy(this%term)
    return
  end function TEMPLATE(potential_get_energy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_get_config)(this, that)
    type(TEMPLATE(potential_t)), target, intent(in) :: this
    type(json_object_t),        pointer             :: that
    !
    call term_get(this%term, that)
    return
  end subroutine TEMPLATE(potential_get_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_get_simulation)(this, that)
    type(TEMPLATE(potential_t)), target, intent(in) :: this
    type(simulation_t),         pointer             :: that
    !
    that=>null()
    if(associated(this%sim))&
      that=>this%sim
    return
  end subroutine TEMPLATE(potential_get_simulation)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_get_system)(this, that)
    type(TEMPLATE(potential_t)), target, intent(in) :: this
    type(system_t),             pointer             :: that
    !
    call term_get(this%term, that)
    return
  end subroutine TEMPLATE(potential_get_system)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_get_geometry)(this, that)
    type(TEMPLATE(potential_t)), target, intent(in) :: this
    type(geometry_t),           pointer             :: that
    !
    call term_get(this%term, that)
    return
  end subroutine TEMPLATE(potential_get_geometry)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_get_density)(this, that)
    type(TEMPLATE(potential_t)), target, intent(in) :: this
    type(density_t),            pointer             :: that
    !
    call term_get(this%term, that)
    return
  end subroutine TEMPLATE(potential_get_density)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_get_potential_1d)(this, that)
    type(TEMPLATE(potential_t)),  intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    !
    that=>null()
    if(this%alloc)&
      call storage_get_storage(this%potential, that)
    return
  end subroutine TEMPLATE(potential_get_potential_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_get_potential_2d)(this, that)
    type(TEMPLATE(potential_t)),    intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    !
    that=>null()
    if(this%alloc)&
      call storage_get_storage(this%potential, that)
    return
  end subroutine TEMPLATE(potential_get_potential_2d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_copy)(this, that)
    type(TEMPLATE(potential_t)), intent(out) :: this
    type(TEMPLATE(potential_t)), intent(in)  :: that
    !
    this%sim=>that%sim
    this%nspin=that%nspin
    this%alloc=that%alloc
    this%default=that%default
    call term_copy(this%term, that%term)
    if(this%alloc)&
      call storage_copy(this%potential, that%potential)
    return
  end subroutine TEMPLATE(potential_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_end)(this)
    type(TEMPLATE(potential_t)), intent(inout) :: this
    !
    if(this%alloc)&
      call storage_end(this%potential)
    call term_end(this%term)
    this%default=0.0_wp
    this%alloc=.false.
    this%nspin=0
    this%sim=>null()
    return
  end subroutine TEMPLATE(potential_end)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_interpolation_init)(this, that)
    type(TEMPLATE(potential_interpolation_t)), intent(out) :: this
    type(TEMPLATE(potential_t)),       target, intent(in)  :: that
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: type, ierr
    !
    cnfg=>null()
    this%self=>that
    if(that%alloc)then
      call term_get(that%term, cnfg)
      ASSERT(associated(cnfg))
      call json_get(cnfg, "interpolation", type, ierr)
      if(ierr==JSON_OK)then
        call storage_interpolation_init(this%intrp, that%potential, type)
      else
        call storage_interpolation_init(this%intrp, that%potential)
      end if
      cnfg=>null()
    end if
    return
  end subroutine TEMPLATE(potential_interpolation_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_interpolation_eval_1d)(this, x, val, ierr)
    type(TEMPLATE(potential_interpolation_t)), intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val
    integer,                     intent(out) :: ierr
    !
    if(this%self%alloc)then
      call storage_interpolation_eval(this%intrp, x, val, ierr)
    else
      val=this%self%default
    end if
    return
  end subroutine TEMPLATE(potential_interpolation_eval_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_interpolation_eval_2d)(this, x, val, ierr)
    type(TEMPLATE(potential_interpolation_t)), intent(in)  :: this
    real(kind=wp),               dimension(:), intent(in)  :: x
    real(kind=wp),               dimension(:), intent(out) :: val
    integer,                                   intent(out) :: ierr
    !
    if(this%self%alloc)then
      call storage_interpolation_eval(this%intrp, x, val, ierr)
    else
      val=this%self%default
    end if
    return
  end subroutine TEMPLATE(potential_interpolation_eval_2d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_interpolation_copy)(this, that)
    type(TEMPLATE(potential_interpolation_t)), intent(out) :: this
    type(TEMPLATE(potential_interpolation_t)), intent(in)  :: that
    !
    this%self=>that%self
    call storage_interpolation_copy(this%intrp, that%intrp)
    return
  end subroutine TEMPLATE(potential_interpolation_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(potential_interpolation_end)(this)
    type(TEMPLATE(potential_interpolation_t)), intent(inout) :: this
    !
    call storage_interpolation_end(this%intrp)
    this%self=>null()
    return
  end subroutine TEMPLATE(potential_interpolation_end)

end module TEMPLATE(potential_m)

!! Local Variables:
!! mode: f90
!! End:
