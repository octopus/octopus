#include "global.h"
#include "template.h"

module TEMPLATE(m)

  use global_m
  use messages_m
  use profiling_m

  use atom_m,    only: atom_t
  use basis_m,   only: basis_t, basis_init, basis_copy, basis_end, basis_to_external, basis_to_internal
  use json_m,    only: JSON_OK, json_object_t, json_get
  use kinds_m,   only: wp
  use space_m,   only: space_t
  use species_m, only: species_t

  use TEMPLATE(calc_m), only:             &
    calc_t      => TEMPLATE(calc_t),      &
    calc_init   => TEMPLATE(calc_init),   &
    calc_start  => TEMPLATE(calc_start),  &
    calc_update => TEMPLATE(calc_update), &
    calc_get    => TEMPLATE(calc_get),    &
    calc_copy   => TEMPLATE(calc_copy),   &
    calc_end    => TEMPLATE(calc_end)

  use TEMPLATE(density_m), only:      &
    density_t => TEMPLATE(density_t)

  use TEMPLATE(geometry_m), only:       &
    geometry_t !=> TEMPLATE(geometry_t)

  use TEMPLATE(hamiltonian_m), only:          &
    hamiltonian_t => TEMPLATE(hamiltonian_t)

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

  use TEMPLATE(simulation_m), only:         &
    simulation_t !=> TEMPLATE(simulation_t)

  use TEMPLATE(states_m), only:     &
    states_t => TEMPLATE(states_t)

  use TEMPLATE(system_m), only:     &
    system_t => TEMPLATE(system_t)

  use TEMPLATE(hamiltonian_m), only:                      &
    get_energy    => TEMPLATE(hamiltonian_get_energy),    &
    get_potential => TEMPLATE(hamiltonian_get_potential)

  use TEMPLATE(calc_m), only:                           &
    calc_iterator_t    => TEMPLATE(calc_iterator_t),    &
    calc_iterator_init => TEMPLATE(calc_iterator_init), &
    calc_iterator_next => TEMPLATE(calc_iterator_next), &
    calc_iterator_copy => TEMPLATE(calc_iterator_copy), &
    calc_iterator_end  => TEMPLATE(calc_iterator_end)

  use TEMPLATE(calc_m), only:                      &
    TEMPLATE(ITER_OK)  => TEMPLATE(CALC_ITER_OK),  &
    TEMPLATE(ITER_END) => TEMPLATE(CALC_ITER_END)

  use TEMPLATE(calc_m), only:                                     &
    calc_interpolation_t    => TEMPLATE(calc_interpolation_t),    &
    calc_interpolation_init => TEMPLATE(calc_interpolation_init), &
    calc_interpolation_eval => TEMPLATE(calc_interpolation_eval), &
    calc_interpolation_copy => TEMPLATE(calc_interpolation_copy), &
    calc_interpolation_end  => TEMPLATE(calc_interpolation_end)

#ifdef SUBTEMPLATE_NAME
  use json_m,  only: json_array_t, json_array_iterator_t, json_init, json_next, json_len, json_end

  use TEMPLATE(calc_m), only:             &
    calc_extend => TEMPLATE(calc_extend)

  use SUBTEMPLATE(m), only:          &
    sub_t     => SUBTEMPLATE(t),     &
    sub_init  => SUBTEMPLATE(init),  &
    sub_start => SUBTEMPLATE(start), &
    sub_get   => SUBTEMPLATE(get),   &
    sub_end   => SUBTEMPLATE(end)

  use SUBTEMPLATE(list_m), only:         &
    list_t    => SUBTEMPLATE(list_t),    &
    list_len  => SUBTEMPLATE(list_len),  &
    list_init => SUBTEMPLATE(list_init), &
    list_push => SUBTEMPLATE(list_push), &
    list_pop  => SUBTEMPLATE(list_pop),  &
    list_copy => SUBTEMPLATE(list_copy), &
    list_end  => SUBTEMPLATE(list_end)

  use SUBTEMPLATE(list_m), only:                           &
    list_iterator_t    => SUBTEMPLATE(list_iterator_t),    &
    list_iterator_init => SUBTEMPLATE(list_iterator_init), &
    list_iterator_next => SUBTEMPLATE(list_iterator_next), &
    list_iterator_end  => SUBTEMPLATE(list_iterator_end)
#endif

  implicit none

  private

  public ::           &
    TEMPLATE(init),   &
    TEMPLATE(start),  &
    TEMPLATE(use),    &
    TEMPLATE(get),    &
    TEMPLATE(copy),   &
    TEMPLATE(end)

  public ::                  &
    TEMPLATE(ITER_OK),       &
    TEMPLATE(ITER_END),      & 
    TEMPLATE(iterator_init), &
    TEMPLATE(iterator_next), &
    TEMPLATE(iterator_copy), &
    TEMPLATE(iterator_end)

  public ::                       &
    TEMPLATE(interpolation_init), &
    TEMPLATE(interpolation_eval), &
    TEMPLATE(interpolation_copy), &
    TEMPLATE(interpolation_end)

  type, public :: TEMPLATE(t)
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(basis_t)                :: basis
    type(calc_t)                 :: calc
#ifdef SUBTEMPLATE_NAME
    type(list_t)                 :: list
#endif
  end type TEMPLATE(t)

  type, public :: TEMPLATE(iterator_t)
    private
    type(TEMPLATE(t)), pointer :: self
    type(calc_iterator_t)      :: iter
  end type TEMPLATE(iterator_t)

  type, public :: TEMPLATE(interpolation_t)
    private
    type(TEMPLATE(t)), pointer :: self =>null()
    type(calc_interpolation_t) :: intrp
  end type TEMPLATE(interpolation_t)

  interface TEMPLATE(get)
    module procedure TEMPLATE(get_config)
    module procedure TEMPLATE(get_basis)
    module procedure TEMPLATE(get_calc)
    module procedure TEMPLATE(get_simulation)
    module procedure TEMPLATE(get_system)
    module procedure TEMPLATE(get_space)
    module procedure TEMPLATE(get_geometry)
    module procedure TEMPLATE(get_states)
    module procedure TEMPLATE(get_density)
    module procedure TEMPLATE(get_hamiltonian)
#ifdef EXTERNAL
    module procedure TEMPLATE(get_external)
#endif
#ifdef HARTREE
    module procedure TEMPLATE(get_hartree)
#endif
#ifdef IONIC
    module procedure TEMPLATE(get_ionic)
#endif
#ifdef TNADD
    module procedure TEMPLATE(get_tnadd)
#endif
  end interface TEMPLATE(get)

  interface TEMPLATE(interpolation_init)
    module procedure TEMPLATE(interpolation_init_density)
#ifdef EXTERNAL
    module procedure TEMPLATE(interpolation_init_external)
#endif
  end interface TEMPLATE(interpolation_init)

  interface TEMPLATE(interpolation_eval)
    module procedure TEMPLATE(interpolation_eval_1d)
    module procedure TEMPLATE(interpolation_eval_2d)
  end interface TEMPLATE(interpolation_eval)

contains

  ! ---------------------------------------------------------
  subroutine TEMPLATE(init_common)(this, space, config)
    type(TEMPLATE(t)),           intent(out) :: this
    type(space_t),       target, intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    nullify(cnfg)
    this%config=>config
    this%space=>space
    call json_get(this%config, "basis", cnfg, ierr)
    if(ierr==JSON_OK)call basis_init(this%basis, cnfg)
    nullify(cnfg)
    call json_get(this%config, "calculation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call calc_init(this%calc, space, cnfg)
    nullify(cnfg)
    return
  end subroutine TEMPLATE(init_common)

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(init)(this, space, config)
    type(TEMPLATE(t)),   intent(out) :: this
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config
    !
    type(json_array_t),  pointer :: list
    integer                      :: ierr
    !
    call TEMPLATE(init_common)(this, space, config)
    nullify(list)
    call json_get(this%config, "list", list, ierr)
    ASSERT(ierr==JSON_OK)
    call list_init(this%list)
    if(json_len(list)>0)&
      call TEMPLATE(list_init)(this, list)
    nullify(list)
    return
  end subroutine TEMPLATE(init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_init)(this, list)
    type(TEMPLATE(t)),  intent(inout) :: this
    type(json_array_t), intent(in)    :: list
    !
    type(sub_t),         pointer :: that
    type(json_object_t), pointer :: cnfg
    type(json_array_iterator_t)  :: iter
    integer                      :: ierr
    !
    nullify(cnfg, that)
    call json_init(iter, list)
    do
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      SAFE_ALLOCATE(that)
      call sub_init(that, this%space, cnfg)
      call list_push(this%list, that)
      call TEMPLATE(extend)(this, that)
      nullify(cnfg, that)
    end do
    call TEMPLATE(extend)(this)
    call json_end(iter)
    return
  end subroutine TEMPLATE(list_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(extend)(this, that)
    type(TEMPLATE(t)),     intent(inout) :: this
    type(sub_t), optional, intent(in)    :: that
    !
    call calc_extend(this%calc, that)
    return
  end subroutine TEMPLATE(extend)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(update)(this, sim)
    type(TEMPLATE(t)),  intent(inout) :: this
    type(simulation_t), intent(in)    :: sim
    !
    type(sub_t),         pointer :: sub
    type(json_object_t), pointer :: cnfg
    type(list_iterator_t)        :: iter
    logical                      :: temp, pass
    integer                      :: ierr
    !
    nullify(sub)
    call list_iterator_init(iter, this%list)
    do
      call list_iterator_next(iter, sub)
      if(.not.associated(sub))exit
      nullify(cnfg)
      call sub_get(sub, cnfg)
      call json_get(cnfg, "temporary", temp, ierr)
      if(ierr/=JSON_OK)temp=.false.
      call json_get(cnfg, "pass_simulation", pass, ierr)
      if(ierr/=JSON_OK)pass=.false.
      if(pass)then
        print *, "***: sub_start(sub, sim)"
        call sub_start(sub, sim)
      else
        print *, "***: sub_start(sub)"
        call sub_start(sub)
      end if
      print *, "***: calc_update(this%calc, sub)"
      call calc_update(this%calc, sub)
      if(temp)then
        print *, "***: sub_end(sub)"
        call sub_end(sub)
      end if
      nullify(sub)
    end do
    call list_iterator_end(iter)
    return
  end subroutine TEMPLATE(update)

  ! ---------------------------------------------------------
  function TEMPLATE(use)(this) result(is)
    type(TEMPLATE(t)), intent(in) :: this
    !
    logical :: is
    !
    is=(list_len(this%list)>0)
    return
  end function TEMPLATE(use)
#else

  ! ---------------------------------------------------------
  subroutine TEMPLATE(init)(this, space, config)
    type(TEMPLATE(t)),   intent(out) :: this
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config
    !
    call TEMPLATE(init_common)(this, space, config)
    return
  end subroutine TEMPLATE(init)

  ! ---------------------------------------------------------
  function TEMPLATE(use)(this) result(is)
    type(TEMPLATE(t)), intent(in) :: this
    !
    logical :: is
    !
    is=.false. !(list_len(this%list)>0)
    return
  end function TEMPLATE(use)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(start)(this, that)
    type(TEMPLATE(t)),            intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: that
    !
    type(simulation_t), pointer :: sim
    !
    nullify(sim)
    print *, "***: calc_start"
    call calc_start(this%calc, that)
#ifdef SUBTEMPLATE_NAME
    print *, "***: calc_get"
    call calc_get(this%calc, sim)
    ASSERT(associated(sim))
    print *, "***: update"
    call TEMPLATE(update)(this, sim)
    nullify(sim)
#endif
    print *, "***: calc_update(this%calc)"
    call calc_update(this%calc)
    return
  end subroutine TEMPLATE(start)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_config)(this, that)
    type(TEMPLATE(t)),    target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    return
  end subroutine TEMPLATE(get_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_basis)(this, that)
    type(TEMPLATE(t)), target, intent(in) :: this
    type(basis_t),    pointer             :: that
    !
    that=>this%basis
    return
  end subroutine TEMPLATE(get_basis)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_calc)(this, that)
    type(TEMPLATE(t)), target, intent(in) :: this
    type(calc_t),     pointer             :: that
    !
    that=>this%calc
    return
  end subroutine TEMPLATE(get_calc)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_simulation)(this, that)
    type(TEMPLATE(t)),   intent(in) :: this
    type(simulation_t), pointer     :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_simulation)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_system)(this, that)
    type(TEMPLATE(t)), intent(in) :: this
    type(system_t),   pointer     :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_system)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_space)(this, that)
    type(TEMPLATE(t)), intent(in) :: this
    type(space_t),    pointer     :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_space)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_geometry)(this, that)
    type(TEMPLATE(t)), intent(in) :: this
    type(geometry_t), pointer     :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_geometry)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_states)(this, that)
    type(TEMPLATE(t)), target, intent(in) :: this
    type(states_t),   pointer             :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_states)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_density)(this, that)
    type(TEMPLATE(t)), target, intent(in) :: this
    type(density_t),  pointer             :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_density)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_hamiltonian)(this, that)
    type(TEMPLATE(t)),    intent(in) :: this
    type(hamiltonian_t), pointer     :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_hamiltonian)

#ifdef EXTERNAL
  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_external)(this, that)
    type(TEMPLATE(t)), intent(in) :: this
    type(external_t), pointer     :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_external)
#endif

#ifdef HARTREE
  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_hartree)(this, that)
    type(TEMPLATE(t)), intent(in) :: this
    type(hartree_t),  pointer     :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_hartree)
#endif

#ifdef IONIC
  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_ionic)(this, that)
    type(TEMPLATE(t)), intent(in) :: this
    type(ionic_t),    pointer     :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_ionic)
#endif

#ifdef TNADD
  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_tnadd)(this, that)
    type(TEMPLATE(t)), intent(in) :: this
    type(tnadd_t),    pointer     :: that
    !
    call calc_get(this%calc, that)
    return
  end subroutine TEMPLATE(get_tnadd)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(copy)(this, that)
    type(TEMPLATE(t)), intent(out) :: this
    type(TEMPLATE(t)), intent(in)  :: that
    !
    this%config=>that%config
    this%space=>that%space
    call basis_copy(this%basis, that%basis)
    call calc_copy(this%calc, that%calc)
#ifdef SUBTEMPLATE_NAME
    call list_copy(this%list, that%list)
#endif
    return
  end subroutine TEMPLATE(copy)

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_end)(this)
    type(TEMPLATE(t)), intent(inout) :: this
    !
    type(sub_t), pointer :: that
    !
    nullify(that)
    do
      call list_pop(this%list, that)
      if(.not.associated(that))exit
      call sub_end(that)
      SAFE_DEALLOCATE_P(that)
      nullify(that)
    end do
    call list_end(this%list)
    return
  end subroutine TEMPLATE(list_end)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(end)(this)
    type(TEMPLATE(t)), intent(inout) :: this
    !
#ifdef SUBTEMPLATE_NAME
    call TEMPLATE(list_end)(this)
#endif
    call calc_end(this%calc)
    call basis_end(this%basis)
    nullify(this%space, this%config)
    return
  end subroutine TEMPLATE(end)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_init)(this, that)
    type(TEMPLATE(iterator_t)), intent(out) :: this
    type(TEMPLATE(t)),  target, intent(in)  :: that
    !
    this%self=>that
    call calc_iterator_init(this%iter, that%calc)
    return
  end subroutine TEMPLATE(iterator_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next)(this, atom, species, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    type(atom_t),               intent(out)   :: atom
    type(species_t),            intent(out)   :: species
    integer,                    intent(out)   :: ierr
    !
    call calc_iterator_next(this%iter, atom, species, ierr)
    if(ierr==TEMPLATE(ITER_OK))&
      call basis_to_external(this%self%basis, atom%x, atom%x)
    return
  end subroutine TEMPLATE(iterator_next)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_copy)(this, that)
    type(TEMPLATE(iterator_t)), intent(out) :: this
    type(TEMPLATE(iterator_t)), intent(in)  :: that
    !
    this%self=>that%self
    call calc_iterator_copy(this%iter, that%iter)
    return
  end subroutine TEMPLATE(iterator_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_end)(this)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    !
    call calc_iterator_end(this%iter)
    this%self=>null()
    return
  end subroutine TEMPLATE(iterator_end)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(interpolation_init_density)(this, that, what)
    type(TEMPLATE(interpolation_t)), intent(out) :: this
    type(TEMPLATE(t)),       target, intent(in)  :: that
    type(density_t),                 intent(in)  :: what
    !
    this%self=>that
    call calc_interpolation_init(this%intrp, that%calc, what)
    return
  end subroutine TEMPLATE(interpolation_init_density)

#ifdef EXTERNAL
  ! ---------------------------------------------------------
  subroutine TEMPLATE(interpolation_init_external)(this, that, what)
    type(TEMPLATE(interpolation_t)), intent(out) :: this
    type(TEMPLATE(t)),       target, intent(in)  :: that
    type(external_t),                intent(in)  :: what
    !
    this%self=>that
    call calc_interpolation_init(this%intrp, that%calc, what)
    return
  end subroutine TEMPLATE(interpolation_init_external)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(interpolation_eval_1d)(this, x, val)
    type(TEMPLATE(interpolation_t)), intent(in)  :: this
    real(kind=wp),     dimension(:), intent(in)  :: x
    real(kind=wp),                   intent(out) :: val
    !
    real(kind=wp), dimension(size(x)) :: y
    !
    call basis_to_internal(this%self%basis, x, y)
    call calc_interpolation_eval(this%intrp, y, val)
    return
  end subroutine TEMPLATE(interpolation_eval_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(interpolation_eval_2d)(this, x, val)
    type(TEMPLATE(interpolation_t)), intent(in)  :: this
    real(kind=wp),     dimension(:), intent(in)  :: x
    real(kind=wp),     dimension(:), intent(out) :: val
    !
    real(kind=wp), dimension(size(x)) :: y
    !
    call basis_to_internal(this%self%basis, x, y)
    call calc_interpolation_eval(this%intrp, y, val)
    return
  end subroutine TEMPLATE(interpolation_eval_2d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(interpolation_copy)(this, that)
    type(TEMPLATE(interpolation_t)),         intent(out) :: this
    type(TEMPLATE(interpolation_t)), target, intent(in)  :: that
    !
    this%self=>that%self
    call calc_interpolation_copy(this%intrp, that%intrp)
    return
  end subroutine TEMPLATE(interpolation_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(interpolation_end)(this)
    type(TEMPLATE(interpolation_t)), intent(inout) :: this
    !
    call calc_interpolation_end(this%intrp)
    this%self=>null()
    return
  end subroutine TEMPLATE(interpolation_end)

end module TEMPLATE(m)

#define MODULE_TYPE TEMPLATE(m)
#define TYPE TEMPLATE(t)
#include "tlist.F90"
#undef MODULE_TYPE
#undef TYPE

!! Local Variables:
!! mode: f90
!! End:
