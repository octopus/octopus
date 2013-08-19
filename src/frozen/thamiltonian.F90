#include "global.h"
#include "template.h"

module TEMPLATE(hamiltonian_m)

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp

#ifdef EXTERNAL
  use TEMPLATE(external_m), only:                               &
    external_t             => TEMPLATE(external_t),             &
    external_init          => TEMPLATE(external_init),          &
    external_start         => TEMPLATE(external_start),         &
    external_update        => TEMPLATE(external_update),        &
    external_get_energy    => TEMPLATE(external_get_energy),    &
    external_get_potential => TEMPLATE(external_get_potential), &
    external_copy          => TEMPLATE(external_copy),          &
    external_end           => TEMPLATE(external_end)

  use TEMPLATE(external_m), only:                                         &
    external_interpolation_t    => TEMPLATE(external_interpolation_t),    &
    external_interpolation_init => TEMPLATE(external_interpolation_init), &
    external_interpolation_eval => TEMPLATE(external_interpolation_eval), &
    external_interpolation_copy => TEMPLATE(external_interpolation_copy), &
    external_interpolation_end  => TEMPLATE(external_interpolation_end)
#endif

#ifdef IONIC
  use TEMPLATE(ionic_m), only:                                &
    ionic_t               => TEMPLATE(ionic_t),               &
    ionic_init            => TEMPLATE(ionic_init),            &
    ionic_update          => TEMPLATE(ionic_update),          &
    ionic_get_energy      => TEMPLATE(ionic_get_energy),      &
    ionic_get_interaction => TEMPLATE(ionic_get_interaction), &
    ionic_copy            => TEMPLATE(ionic_copy),            &
    ionic_end             => TEMPLATE(ionic_end)
#endif

#ifdef TNADD
  use TEMPLATE(tnadd_m), only:                            &
    tnadd_t             => TEMPLATE(tnadd_t),             &
    tnadd_init          => TEMPLATE(tnadd_init),          &
    tnadd_start         => TEMPLATE(tnadd_start),         &
    tnadd_update        => TEMPLATE(tnadd_update),        &
    tnadd_calc          => TEMPLATE(tnadd_calc),          &
    tnadd_get_energy    => TEMPLATE(tnadd_get_energy),    &
    tnadd_get_potential => TEMPLATE(tnadd_get_potential), &
    tnadd_copy          => TEMPLATE(tnadd_copy),          &
    tnadd_end           => TEMPLATE(tnadd_end)
#endif

  use TEMPLATE(simulation_m), only:         &
    simulation_t !=> TEMPLATE(simulation_t)

  use TEMPLATE(system_m), only:     &
    system_t => TEMPLATE(system_t)

#ifdef SUBTEMPLATE_NAME
  use SUBTEMPLATE(m), only:  &
    sub_t => SUBTEMPLATE(t)
#ifdef EXTERNAL
  use TEMPLATE(external_m), only:                 &
    external_extend => TEMPLATE(external_extend)
#endif
#endif

  implicit none

  private
  public ::                                &
    TEMPLATE(hamiltonian_init),            &
    TEMPLATE(hamiltonian_start),           &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(hamiltonian_extend),          &
#endif
    TEMPLATE(hamiltonian_update),          &
    TEMPLATE(hamiltonian_get),             &
    TEMPLATE(hamiltonian_get_energy),      &
    TEMPLATE(hamiltonian_get_potential),   &
    TEMPLATE(hamiltonian_get_interaction), &
    TEMPLATE(hamiltonian_copy),            &
    TEMPLATE(hamiltonian_end)
  
  public ::                                   &
    TEMPLATE(hamiltonian_interpolation_init), &
    TEMPLATE(hamiltonian_interpolation_eval), &
    TEMPLATE(hamiltonian_interpolation_copy), &
    TEMPLATE(hamiltonian_interpolation_end)

  type, public :: TEMPLATE(hamiltonian_t)
    private
    type(json_object_t), pointer :: config => null()
    type(system_t),      pointer :: sys    => null()
    type(simulation_t),  pointer :: sim    => null()
#ifdef EXTERNAL
    type(external_t)             :: ep
#endif
#ifdef IONIC
    type(ionic_t)                :: ionic
#endif
#ifdef TNADD
    type(tnadd_t)                :: tnadd
#endif
  end type TEMPLATE(hamiltonian_t)

  integer, parameter :: NONE = 0
  integer, parameter :: IEXT = 1

  type, public :: TEMPLATE(hamiltonian_interpolation_t)
    private
    type(TEMPLATE(hamiltonian_t)),  pointer :: self  =>null()
#ifdef EXTERNAL
    type(external_interpolation_t), pointer :: xntrp =>null()
#endif
    integer                                 :: type  = NONE
  end type TEMPLATE(hamiltonian_interpolation_t)

  interface TEMPLATE(hamiltonian_get)
#ifdef EXTERNAL
    module procedure TEMPLATE(hamiltonian_get_external)
#endif
#ifdef IONIC
    module procedure TEMPLATE(hamiltonian_get_ionic)
#endif
#ifdef TNADD
    module procedure TEMPLATE(hamiltonian_get_tnadd)
#endif
  end interface TEMPLATE(hamiltonian_get)

  interface TEMPLATE(hamiltonian_get_energy)
    module procedure TEMPLATE(hamiltonian_get_energy_total)
#ifdef EXTERNAL
    module procedure external_get_energy
#endif
#ifdef IONIC
    module procedure ionic_get_energy
#endif
#ifdef TNADD
    module procedure tnadd_get_energy
#endif
  end interface TEMPLATE(hamiltonian_get_energy)

  interface TEMPLATE(hamiltonian_get_potential)
    !module procedure TEMPLATE(hamiltonian_get_potential_total)
#ifdef EXTERNAL
    module procedure external_get_potential
#endif
#ifdef TNADD
    module procedure tnadd_get_potential
#endif
  end interface TEMPLATE(hamiltonian_get_potential)

  interface TEMPLATE(hamiltonian_get_interaction)
    !module procedure TEMPLATE(hamiltonian_get_potential_total)
#ifdef IONIC
    module procedure ionic_get_interaction
#endif
  end interface TEMPLATE(hamiltonian_get_interaction)

  interface TEMPLATE(hamiltonian_update)
#ifdef SUBTEMPLATE_NAME
    module procedure TEMPLATE(hamiltonian_update_build)
#endif
    module procedure TEMPLATE(hamiltonian_update_finish)
  end interface TEMPLATE(hamiltonian_update)

  interface TEMPLATE(hamiltonian_interpolation_init)
#ifdef EXTERNAL
    module procedure TEMPLATE(hamiltonian_interpolation_init_external)
#endif
  end interface TEMPLATE(hamiltonian_interpolation_init)

contains

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_init)(this, sys, config)
    type(TEMPLATE(hamiltonian_t)), intent(out) :: this
    type(system_t),        target, intent(in)  :: sys
    type(json_object_t),   target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    this%config=>config
    this%sys=>sys
    this%sim=>null()
    nullify(cnfg)
#ifdef EXTERNAL
    call json_get(this%config, "external", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    print *, "***: external_init"
    call external_init(this%ep, sys, cnfg)
    nullify(cnfg)
#endif
#ifdef IONIC
    call json_get(this%config, "ionic", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    print *, "***: ionic_init"
    call ionic_init(this%ionic, sys, cnfg)
    nullify(cnfg)
#endif
#ifdef TNADD
    call json_get(this%config, "tnadd", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    print *, "***: tnadd_init"
    call tnadd_init(this%tnadd, sys, cnfg)
    nullify(cnfg)
#endif
    return
  end subroutine TEMPLATE(hamiltonian_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_start)(this, sim)
    type(TEMPLATE(hamiltonian_t)), intent(inout) :: this
    type(simulation_t),    target, intent(in)    :: sim
    !
    type(json_object_t) :: config
    !
    this%sim=>sim
#ifdef EXTERNAL
    call external_start(this%ep, sim)
#endif
#ifdef TNADD
    call tnadd_start(this%tnadd, sim)
#endif
    return
  end subroutine TEMPLATE(hamiltonian_start)

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_extend)(this, that)
    type(TEMPLATE(hamiltonian_t)), intent(inout) :: this
    type(sub_t),                   intent(in)    :: that
    !
#ifdef EXTERNAL
    print *, "***: external_extend"
    call external_extend(this%ep, that)
#endif
    return
  end subroutine TEMPLATE(hamiltonian_extend)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_update_build)(this, that)
    type(TEMPLATE(hamiltonian_t)), intent(inout) :: this
    type(sub_t),                   intent(in)    :: that
    !
#ifdef EXTERNAL
    call external_update(this%ep, that)
#endif
    return
  end subroutine TEMPLATE(hamiltonian_update_build)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_update_finish)(this)
    type(TEMPLATE(hamiltonian_t)), intent(inout) :: this
    !
#ifdef EXTERNAL
    print *, "***: hamiltonian_update_finish: external_update"
    call external_update(this%ep)
#endif
#ifdef IONIC
    print *, "***: hamiltonian_update_finish: ionic_update"
    call ionic_update(this%ionic)
#endif
#ifdef TNADD
    print *, "***: hamiltonian_update_finish: tnadd_update"
    call tnadd_update(this%tnadd)
#endif
    return
  end subroutine TEMPLATE(hamiltonian_update_finish)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(hamiltonian_get_energy_total)(this) result(that)
    type(TEMPLATE(hamiltonian_t)), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=0.0_wp
    return
  end function TEMPLATE(hamiltonian_get_energy_total)

#ifdef EXTERNAL
  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_get_external)(this, that)
    type(TEMPLATE(hamiltonian_t)), target, intent(in) :: this
    type(external_t),             pointer             :: that
    !
    that=>this%ep
    return
  end subroutine TEMPLATE(hamiltonian_get_external)
#endif

#ifdef IONIC
  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_get_ionic)(this, that)
    type(TEMPLATE(hamiltonian_t)), target, intent(in) :: this
    type(ionic_t),                pointer             :: that
    !
    that=>this%ionic
    return
  end subroutine TEMPLATE(hamiltonian_get_ionic)
#endif

#ifdef TNADD
  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_get_tnadd)(this, that)
    type(TEMPLATE(hamiltonian_t)), target, intent(in) :: this
    type(tnadd_t),                pointer             :: that
    !
    that=>this%tnadd
    return
  end subroutine TEMPLATE(hamiltonian_get_tnadd)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_copy)(this_out, this_in)
    type(TEMPLATE(hamiltonian_t)), intent(out) :: this_out
    type(TEMPLATE(hamiltonian_t)), intent(in)  :: this_in
    !
    this_out%config=>this_in%config
    this_out%sys=>this_in%sys
    this_out%sim=>this_in%sim
#ifdef EXTERNAL
    call external_copy(this_out%ep, this_in%ep)
#endif
#ifdef IONIC
    call ionic_copy(this_out%ionic, this_in%ionic)
#endif
#ifdef TNADD
    call tnadd_copy(this_out%tnadd, this_in%tnadd)
#endif
    return
  end subroutine TEMPLATE(hamiltonian_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_end)(this)
    type(TEMPLATE(hamiltonian_t)), intent(inout) :: this
    !
#ifdef TNADD
    call tnadd_end(this%tnadd)
#endif
#ifdef IONIC
    call ionic_end(this%ionic)
#endif
#ifdef EXTERNAL
    call external_end(this%ep)
#endif
    nullify(this%sim, this%sys, this%config)
    return
  end subroutine TEMPLATE(hamiltonian_end)

#ifdef EXTERNAL
  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_interpolation_init_external)(this, that, what)
    type(TEMPLATE(hamiltonian_interpolation_t)), intent(out) :: this
    type(TEMPLATE(hamiltonian_t)),               intent(in)  :: that
    type(external_t),                            intent(in)  :: what
    !
    nullify(this%self, this%xntrp)
#ifdef EXTERNAL
    SAFE_ALLOCATE(this%xntrp)
    call external_interpolation_init(this%xntrp, that%ep)
    this%type=IEXT
#else
    this%type=NONE
#endif
    return
  end subroutine TEMPLATE(hamiltonian_interpolation_init_external)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_interpolation_eval)(this, x, val)
    type(TEMPLATE(hamiltonian_interpolation_t)), intent(in)  :: this
    real(kind=wp),                 dimension(:), intent(in)  :: x
    real(kind=wp),                               intent(out) :: val
    !
    select case(this%type)
#ifdef EXTERNAL
    case(IEXT)
      call external_interpolation_eval(this%xntrp, x, val)
#endif
    case default
      ASSERT(.false.)
    end select
    return
  end subroutine TEMPLATE(hamiltonian_interpolation_eval)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_interpolation_copy)(this, that)
    type(TEMPLATE(hamiltonian_interpolation_t)),         intent(out) :: this
    type(TEMPLATE(hamiltonian_interpolation_t)), target, intent(in)  :: that
    !
    this%self=>that%self
    select case(that%type)
#ifdef EXTERNAL
    case(IEXT)
      call external_interpolation_copy(this%xntrp, that%xntrp)
#endif
    case default
      ASSERT(.false.)
    end select
    this%type=that%type
    return
  end subroutine TEMPLATE(hamiltonian_interpolation_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hamiltonian_interpolation_end)(this)
    type(TEMPLATE(hamiltonian_interpolation_t)), intent(inout) :: this
    !
    select case(this%type)
#ifdef EXTERNAL
    case(IEXT)
      call external_interpolation_end(this%xntrp)
      SAFE_DEALLOCATE_P(this%xntrp)
#endif
    case default
      ASSERT(.false.)
    end select
    this%type=NONE
    nullify(this%xntrp, this%self)
    return
  end subroutine TEMPLATE(hamiltonian_interpolation_end)

end module TEMPLATE(hamiltonian_m)

!! Local Variables:
!! mode: f90
!! End:
