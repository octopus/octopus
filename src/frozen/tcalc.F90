#include "global.h"
#include "template.h"

module TEMPLATE(calc_m)

  use global_m
  use messages_m
  use profiling_m

  use basis_m,    only: basis_t
  use geometry_m, only: geometry_t
  use json_m,     only: JSON_OK, json_object_t, json_get
  use kinds_m,    only: wp
  use space_m,    only: space_t

  use TEMPLATE(density_m), only:      &
    density_t => TEMPLATE(density_t)

  use TEMPLATE(hamiltonian_m), only:                    &
    hamiltonian_t      => TEMPLATE(hamiltonian_t),      &
    hamiltonian_init   => TEMPLATE(hamiltonian_init),   &
    hamiltonian_start  => TEMPLATE(hamiltonian_start),  &
    hamiltonian_update => TEMPLATE(hamiltonian_update), &
    hamiltonian_get    => TEMPLATE(hamiltonian_get),    &
    hamiltonian_copy   => TEMPLATE(hamiltonian_copy),   &
    hamiltonian_end    => TEMPLATE(hamiltonian_end)

  use TEMPLATE(hamiltonian_m), only:                                            &
    hamiltonian_interpolation_t    => TEMPLATE(hamiltonian_interpolation_t),    &
    hamiltonian_interpolation_init => TEMPLATE(hamiltonian_interpolation_init), &
    hamiltonian_interpolation_eval => TEMPLATE(hamiltonian_interpolation_eval), &
    hamiltonian_interpolation_copy => TEMPLATE(hamiltonian_interpolation_copy), &
    hamiltonian_interpolation_end  => TEMPLATE(hamiltonian_interpolation_end)

#ifdef EXTERNAL
  use TEMPLATE(external_m), only:       &
    external_t => TEMPLATE(external_t)
#endif

#ifdef HARTREE
  use TEMPLATE(hartree_m), only:      &
    hartree_t => TEMPLATE(hartree_t)
#endif

#ifdef IONIC
  use TEMPLATE(ionic_m), only:    &
    ionic_t => TEMPLATE(ionic_t)
#endif

#ifdef TNADD
  use TEMPLATE(tnadd_m), only:    &
    tnadd_t => TEMPLATE(tnadd_t)
#endif

  use TEMPLATE(simulation_m), only:                 &
    simulation_t, &     !=> TEMPLATE(simulation_t),     &
    simulation_init  => TEMPLATE(simulation_init),  &
    simulation_start => TEMPLATE(simulation_start), &
    simulation_copy  => TEMPLATE(simulation_copy),  &
    simulation_end   => TEMPLATE(simulation_end)

  use TEMPLATE(states_m), only:     &
    states_t => TEMPLATE(states_t)

  use TEMPLATE(system_m), only:               &
    system_t      => TEMPLATE(system_t),      &
    system_init   => TEMPLATE(system_init),   &
    system_start  => TEMPLATE(system_start),  &
    system_update => TEMPLATE(system_update), &
    system_get    => TEMPLATE(system_get),    &
    system_copy   => TEMPLATE(system_copy),   &
    system_end    => TEMPLATE(system_end)

  use TEMPLATE(system_m), only:                             &
    system_iterator_init => TEMPLATE(system_iterator_init)

  use TEMPLATE(system_m), only:                                     &
    TEMPLATE(CALC_ITER_OK)       => TEMPLATE(SYSTEM_ITER_OK),       &
    TEMPLATE(CALC_ITER_END)      => TEMPLATE(SYSTEM_ITER_END),      & 
    TEMPLATE(calc_iterator_t)    => TEMPLATE(system_iterator_t),    &
    TEMPLATE(calc_iterator_next) => TEMPLATE(system_iterator_next), &
    TEMPLATE(calc_iterator_copy) => TEMPLATE(system_iterator_copy), &
    TEMPLATE(calc_iterator_end)  => TEMPLATE(system_iterator_end)

  use TEMPLATE(system_m), only:                                       &
    system_interpolation_t    => TEMPLATE(system_interpolation_t),    &
    system_interpolation_init => TEMPLATE(system_interpolation_init), &
    system_interpolation_eval => TEMPLATE(system_interpolation_eval), &
    system_interpolation_copy => TEMPLATE(system_interpolation_copy), &
    system_interpolation_end  => TEMPLATE(system_interpolation_end)

#ifdef SUBTEMPLATE_NAME
  use TEMPLATE(hamiltonian_m), only:                    &
    hamiltonian_extend => TEMPLATE(hamiltonian_extend)

  use TEMPLATE(simulation_m), only:                   &
    simulation_extend => TEMPLATE(simulation_extend)

  use TEMPLATE(system_m), only:               &
    system_extend => TEMPLATE(system_extend)

  use SUBTEMPLATE(m), only:  &
    sub_t => SUBTEMPLATE(t)
#endif

  implicit none

  private
  public ::                &
    TEMPLATE(calc_init),   &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(calc_extend), &
#endif
    TEMPLATE(calc_update), &
    TEMPLATE(calc_start),  &
    TEMPLATE(calc_get),    &
    TEMPLATE(calc_copy),   &
    TEMPLATE(calc_end)

  public ::                       &
    TEMPLATE(CALC_ITER_OK),       &
    TEMPLATE(CALC_ITER_END),      & 
    TEMPLATE(calc_iterator_t),    &
    TEMPLATE(calc_iterator_init), &
    TEMPLATE(calc_iterator_next), &
    TEMPLATE(calc_iterator_copy), &
    TEMPLATE(calc_iterator_end)

  public ::                            &
    TEMPLATE(calc_interpolation_init), &
    TEMPLATE(calc_interpolation_eval), &
    TEMPLATE(calc_interpolation_copy), &
    TEMPLATE(calc_interpolation_end)

  integer, parameter :: NONE = 0
  integer, parameter :: ISYS = 1
  integer, parameter :: IHAM = 2

  type, public :: TEMPLATE(calc_t)
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(simulation_t)           :: ism
    type(system_t)               :: sys
    type(hamiltonian_t)          :: hm
  end type TEMPLATE(calc_t)

  type, public :: TEMPLATE(calc_interpolation_t)
    private
    type(TEMPLATE(calc_t)),            pointer :: self  =>null()
    type(system_interpolation_t),      pointer :: sntrp =>null()
    type(hamiltonian_interpolation_t), pointer :: hntrp =>null()
    integer                                    :: type  = NONE
  end type TEMPLATE(calc_interpolation_t)

  interface TEMPLATE(calc_get)
    module procedure TEMPLATE(calc_get_config)
    module procedure TEMPLATE(calc_get_simulation)
    module procedure TEMPLATE(calc_get_system)
    module procedure TEMPLATE(calc_get_space)
    module procedure TEMPLATE(calc_get_geometry)
    module procedure TEMPLATE(calc_get_states)
    module procedure TEMPLATE(calc_get_density)
    module procedure TEMPLATE(calc_get_hamiltonian)
#ifdef EXTERNAL
    module procedure TEMPLATE(calc_get_external)
#endif
#ifdef HARTREE
    module procedure TEMPLATE(calc_get_hartree)
#endif
#ifdef IONIC
    module procedure TEMPLATE(calc_get_ionic)
#endif
#ifdef TNADD
    module procedure TEMPLATE(calc_get_tnadd)
#endif
  end interface TEMPLATE(calc_get)

  interface TEMPLATE(calc_update)
#ifdef SUBTEMPLATE_NAME
    module procedure TEMPLATE(calc_update_build)
#endif
    module procedure TEMPLATE(calc_update_finalize)
  end interface TEMPLATE(calc_update)

  ! ---------------------------------------------------------
  interface TEMPLATE(calc_interpolation_init)
    module procedure TEMPLATE(calc_interpolation_init_density)
#ifdef EXTERNAL
    module procedure TEMPLATE(calc_interpolation_init_external)
#endif
  end interface TEMPLATE(calc_interpolation_init)

  interface TEMPLATE(calc_interpolation_eval)
    module procedure TEMPLATE(calc_interpolation_eval_1d)
    module procedure TEMPLATE(calc_interpolation_eval_2d)
  end interface TEMPLATE(calc_interpolation_eval)

contains

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_init)(this, space, config)
    type(TEMPLATE(calc_t)), target, intent(out) :: this
    type(space_t),          target, intent(in)  :: space
    type(json_object_t),    target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    type(geometry_t),    pointer :: geo
    integer                      :: ierr
    !
    this%config=>config
    this%space=>space
    this%sim=>this%ism
    nullify(cnfg)
    call json_get(this%config, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    print *, "system_init"
    call system_init(this%sys, space, cnfg)
    nullify(cnfg, geo)
    call system_get(this%sys, geo)
    call json_get(config, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    print *, "simulation_init"
    call simulation_init(this%sim, geo, space, cnfg)
    nullify(cnfg, geo)
    call json_get(this%config, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    print *, "hamiltonian_init"
    call hamiltonian_init(this%hm, this%sys, cnfg)
    nullify(cnfg)
    return
  end subroutine TEMPLATE(calc_init)

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_extend)(this, that)
    type(TEMPLATE(calc_t)), intent(inout) :: this
    type(sub_t),  optional, intent(in)    :: that
    !
    print *, "simulation_extend"
    call simulation_extend(this%sim, that)
    print *, "system_extend"
    call system_extend(this%sys, that)
    print *, "hamiltonian_extend"
    call hamiltonian_extend(this%hm, that)
    return
  end subroutine TEMPLATE(calc_extend)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_update_build)(this, that)
    type(TEMPLATE(calc_t)), intent(inout) :: this
    type(sub_t),            intent(in)    :: that
    !
    print *, "***: system_update"
    call system_update(this%sys, that)
    print *, "***: hamiltonian_update"
    call hamiltonian_update(this%hm, that)
    return
  end subroutine TEMPLATE(calc_update_build)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_update_finalize)(this)
    type(TEMPLATE(calc_t)), intent(inout) :: this
    !
    print *, "***: system_update_finalize"
    call system_update(this%sys)
    print *, "***: hamiltonian_update_finalize"
    call hamiltonian_update(this%hm)
    return
  end subroutine TEMPLATE(calc_update_finalize)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_start)(this, sim)
    type(TEMPLATE(calc_t)),               intent(inout) :: this
    type(simulation_t), optional, target, intent(in)    :: sim
    !
    if(present(sim))then
      this%sim=>sim
      print *, "***: simulation_end"
      call simulation_end(this%ism)
    else
      print *, "***: simulation_start"
      call simulation_start(this%sim)
    end if
    print *, "***: system_start"
    call system_start(this%sys, this%sim)
    print *, "***: hamiltonian_start"
    call hamiltonian_start(this%hm, this%sim)
    return
  end subroutine TEMPLATE(calc_start)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_config)(this, that)
    type(TEMPLATE(calc_t)), target, intent(in) :: this
    type(json_object_t),   pointer             :: that
    !
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    return
  end subroutine TEMPLATE(calc_get_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_simulation)(this, that)
    type(TEMPLATE(calc_t)), target, intent(in) :: this
    type(simulation_t),    pointer             :: that
    !
    that=>this%sim
    return
  end subroutine TEMPLATE(calc_get_simulation)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_system)(this, that)
    type(TEMPLATE(calc_t)), target, intent(in) :: this
    type(system_t),        pointer             :: that
    !
    that=>this%sys
    return
  end subroutine TEMPLATE(calc_get_system)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_space)(this, that)
    type(TEMPLATE(calc_t)), target, intent(in) :: this
    type(space_t),         pointer             :: that
    !
    call system_get(this%sys, that)
    return
  end subroutine TEMPLATE(calc_get_space)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_geometry)(this, that)
    type(TEMPLATE(calc_t)), target, intent(in) :: this
    type(geometry_t),      pointer             :: that
    !
    call system_get(this%sys, that)
    return
  end subroutine TEMPLATE(calc_get_geometry)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_states)(this, that)
    type(TEMPLATE(calc_t)),  target, intent(in) :: this
    type(states_t),         pointer             :: that
    !
    call system_get(this%sys, that)
    return
  end subroutine TEMPLATE(calc_get_states)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_density)(this, that)
    type(TEMPLATE(calc_t)),  target, intent(in) :: this
    type(density_t),        pointer             :: that
    !
    call system_get(this%sys, that)
    return
  end subroutine TEMPLATE(calc_get_density)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_hamiltonian)(this, that)
    type(TEMPLATE(calc_t)), target, intent(in) :: this
    type(hamiltonian_t),   pointer             :: that
    !
    that=>this%hm
    return
  end subroutine TEMPLATE(calc_get_hamiltonian)

#ifdef EXTERNAL
  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_external)(this, that)
    type(TEMPLATE(calc_t)), intent(in) :: this
    type(external_t),      pointer     :: that
    !
    call hamiltonian_get(this%hm, that)
    return
  end subroutine TEMPLATE(calc_get_external)
#endif

#ifdef HARTREE
  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_hartree)(this, that)
    type(TEMPLATE(calc_t)), intent(in) :: this
    type(hartree_t),       pointer     :: that
    !
    call hamiltonian_get(this%hm, that)
    return
  end subroutine TEMPLATE(calc_get_hartree)
#endif

#ifdef IONIC
  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_ionic)(this, that)
    type(TEMPLATE(calc_t)), intent(in) :: this
    type(ionic_t),         pointer     :: that
    !
    call hamiltonian_get(this%hm, that)
    return
  end subroutine TEMPLATE(calc_get_ionic)
#endif

#ifdef TNADD
  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_get_tnadd)(this, that)
    type(TEMPLATE(calc_t)), intent(in) :: this
    type(tnadd_t),         pointer     :: that
    !
    call hamiltonian_get(this%hm, that)
    return
  end subroutine TEMPLATE(calc_get_tnadd)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_copy)(this_out, this_in)
    type(TEMPLATE(calc_t)), intent(out) :: this_out
    type(TEMPLATE(calc_t)), intent(in)  :: this_in
    !
    this_out%config=>this_in%config
    this_out%sim=>this_in%sim
    call simulation_copy(this_out%ism, this_in%ism)
    call system_copy(this_out%sys, this_in%sys)
    call hamiltonian_copy(this_out%hm, this_in%hm)
    return
  end subroutine TEMPLATE(calc_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_end)(this)
    type(TEMPLATE(calc_t)), intent(inout) :: this
    !
    call hamiltonian_end(this%hm)
    call system_end(this%sys)
    call simulation_end(this%ism)
    nullify(this%config, this%sim)
    return
  end subroutine TEMPLATE(calc_end)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_iterator_init)(this, that)
    type(TEMPLATE(calc_iterator_t)), intent(out) :: this
    type(TEMPLATE(calc_t)),          intent(in)  :: that
    !
    call system_iterator_init(this, that%sys)
    return
  end subroutine TEMPLATE(calc_iterator_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_interpolation_init_density)(this, that, what)
    type(TEMPLATE(calc_interpolation_t)), intent(out) :: this
    type(TEMPLATE(calc_t)),       target, intent(in)  :: that
    type(density_t),                      intent(in)  :: what
    !
    this%self=>that
    SAFE_ALLOCATE(this%sntrp)
    call system_interpolation_init(this%sntrp, that%sys, what)
    this%hntrp=>null()
    this%type=ISYS
    print *, "***: calc_interpolation_init_density: ", this%type
    return
  end subroutine TEMPLATE(calc_interpolation_init_density)

#ifdef EXTERNAL
  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_interpolation_init_external)(this, that, what)
    type(TEMPLATE(calc_interpolation_t)), intent(out) :: this
    type(TEMPLATE(calc_t)),       target, intent(in)  :: that
    type(external_t),                     intent(in)  :: what
    !
    this%self=>that
    this%sntrp=>null()
    SAFE_ALLOCATE(this%hntrp)
    call hamiltonian_interpolation_init(this%hntrp, that%hm, what)
    this%type=IHAM
    print *, "***: calc_interpolation_init_external: ", this%type
    return
  end subroutine TEMPLATE(calc_interpolation_init_external)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_interpolation_eval_1d)(this, x, val)
    type(TEMPLATE(calc_interpolation_t)), intent(in)  :: this
    real(kind=wp),          dimension(:), intent(in)  :: x
    real(kind=wp),                        intent(out) :: val
    !
    select case(this%type)
    case(ISYS)
      call system_interpolation_eval(this%sntrp, x, val)
    case(IHAM)
      call hamiltonian_interpolation_eval(this%hntrp, x, val)
    case default
      print *, "***: calc_interpolation_eval_1d: ", this%type
      ASSERT(.false.)
    end select
    return
  end subroutine TEMPLATE(calc_interpolation_eval_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_interpolation_eval_2d)(this, x, val)
    type(TEMPLATE(calc_interpolation_t)), intent(in)  :: this
    real(kind=wp),          dimension(:), intent(in)  :: x
    real(kind=wp),          dimension(:), intent(out) :: val
    !
    select case(this%type)
    case(ISYS)
      call system_interpolation_eval(this%sntrp, x, val)
    case default
      print *, "***: calc_interpolation_eval_2d: ", this%type
      ASSERT(.false.)
    end select
    return
  end subroutine TEMPLATE(calc_interpolation_eval_2d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_interpolation_copy)(this, that)
    type(TEMPLATE(calc_interpolation_t)),         intent(out) :: this
    type(TEMPLATE(calc_interpolation_t)), target, intent(in)  :: that
    !
    this%self=>that%self
    select case(that%type)
    case(ISYS)
      call system_interpolation_copy(this%sntrp, that%sntrp)
    case(IHAM)
      call hamiltonian_interpolation_copy(this%hntrp, that%hntrp)
    case default
      ASSERT(.false.)
    end select
    this%type=that%type
    return
  end subroutine TEMPLATE(calc_interpolation_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(calc_interpolation_end)(this)
    type(TEMPLATE(calc_interpolation_t)), intent(inout) :: this
    !
    select case(this%type)
    case(ISYS)
      call system_interpolation_end(this%sntrp)
      SAFE_DEALLOCATE_P(this%sntrp)
    case(IHAM)
      call hamiltonian_interpolation_end(this%hntrp)
      SAFE_DEALLOCATE_P(this%hntrp)
    case default
      ASSERT(.false.)
    end select
    this%type=NONE
    nullify(this%sntrp, this%hntrp, this%self)
    return
  end subroutine TEMPLATE(calc_interpolation_end)

end module TEMPLATE(calc_m)

!! Local Variables:
!! mode: f90
!! End:
