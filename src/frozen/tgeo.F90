#include "global.h"
#include "template.h"

module TEMPLATE(geo_m)

  use global_m
  use messages_m
  use profiling_m

  use atom_m,        only: atom_t, atom_end, atom_set_species, atom_get_species
  use distributed_m, only: distributed_nullify
  use json_m !,        only: json_object_t
  use space_m,       only: space_t
  use species_m,     only: species_t, species_end, species_set_index, species_index

  use TEMPLATE(geometry_m), only:             &
    geometry_t    !=> TEMPLATE(geometry_t),    &

  use TEMPLATE(geometry_m), only:             &
    geometry_init => TEMPLATE(geometry_init), &
    geometry_copy => TEMPLATE(geometry_copy), &
    geometry_end  => TEMPLATE(geometry_end)

#ifdef SUBTEMPLATE_NAME
  use atom_list_m,  only:        &
    list_t    => atom_list_t,    &
    list_len  => atom_list_len,  &
    list_init => atom_list_init, &
    list_push => atom_list_push, &
    list_pop  => atom_list_pop,  &
    list_copy => atom_list_copy, &
    list_end  => atom_list_end

  use atom_list_m,  only:                          &
    list_iterator_t    => atom_list_iterator_t,    &
    list_iterator_init => atom_list_iterator_init, &
    list_iterator_next => atom_list_iterator_next, &
    list_iterator_end  => atom_list_iterator_end

  use atomspec_table_m, only:          &
    table_t    => atomspec_table_t,    &
    table_len  => atomspec_table_len,  &
    table_init => atomspec_table_init, &
    table_pop  => atomspec_table_pop,  &
    table_set  => atomspec_table_set,  &
    table_get  => atomspec_table_get,  &
    table_copy => atomspec_table_copy, &
    table_end  => atomspec_table_end

  use SUBTEMPLATE(m), only:      &
    sub_t   => SUBTEMPLATE(t),   &
    sub_get => SUBTEMPLATE(get)

  use SUBTEMPLATE(m), only:                          &
    SUB_ITER_OK       => SUBTEMPLATE(ITER_OK),       &
    sub_iterator_t    => SUBTEMPLATE(iterator_t),    &
    sub_iterator_init => SUBTEMPLATE(iterator_init), &
    sub_iterator_next => SUBTEMPLATE(iterator_next), &
    sub_iterator_end  => SUBTEMPLATE(iterator_end)
#endif

  implicit none

  private
  public ::               &
    TEMPLATE(geo_init),   &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(geo_extend), &
#endif
    TEMPLATE(geo_get),    &
    TEMPLATE(geo_copy),   &
    TEMPLATE(geo_end)

  public ::                      &
    TEMPLATE(geo_iterator_init), &
    TEMPLATE(geo_iterator_next), &
    TEMPLATE(geo_iterator_copy), &
    TEMPLATE(geo_iterator_end)

  integer, public :: TEMPLATE(GEO_ITER_OK)  = 0
  integer, public :: TEMPLATE(GEO_ITER_END) = 1

  integer, parameter :: TABLE_INIT_LEN = 11

  type, public :: TEMPLATE(geo_t)
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    logical                      :: init   = .false.
    type(geometry_t)             :: geo
#ifdef SUBTEMPLATE_NAME
    type(list_t)                 :: list
    type(table_t)                :: table
#endif
  end type TEMPLATE(geo_t)

  type, public :: TEMPLATE(geo_iterator_t)
    private
    integer                                :: n        = 0
    integer                                :: natoms   = 0
    integer                                :: nspecies = 0
    type(atom_t),    dimension(:), pointer :: atom     =>null()
    type(species_t), dimension(:), pointer :: species  =>null()
  end type TEMPLATE(geo_iterator_t)

contains

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_init)(this, space, config)
    type(TEMPLATE(geo_t)),                 intent(out) :: this
    type(space_t),                 target, intent(in)  :: space
    type(json_object_t), optional, target, intent(in)  :: config
    !
    this%config=>null()
    if(present(config))this%config=>config
    this%space=>space
    this%init=.false.
    call list_init(this%list)
    call table_init(this%table, TABLE_INIT_LEN)
    return
  end subroutine TEMPLATE(geo_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_extend_build)(this, that)
    type(TEMPLATE(geo_t)), intent(inout) :: this
    type(sub_t),           intent(in)    :: that
    !
    type(atom_t),    pointer :: atom
    type(species_t), pointer :: spec
    type(sub_iterator_t)     :: iter
    integer                  :: ierr
    !
    nullify(atom, spec)
    call sub_iterator_init(iter, that)
    do
      SAFE_ALLOCATE(atom)
      SAFE_ALLOCATE(spec)
      call sub_iterator_next(iter, atom, spec, ierr)
      if(ierr/=SUB_ITER_OK)then
        SAFE_DEALLOCATE_P(atom)
        SAFE_DEALLOCATE_P(spec)
        exit
      end if
      call list_push(this%list, atom)
      call table_set(this%table, atom, spec)
    end do
    call sub_iterator_end(iter)
    nullify(atom, spec)
    return
  end subroutine TEMPLATE(geo_extend_build)

  ! ---------------------------------------------------------
  subroutine geometry_extend_finalize(this, space, list, table)
    type(geometry_t),      intent(inout) :: this
    type(space_t), target, intent(in)    :: space
    type(list_t),          intent(inout) :: list
    type(table_t),         intent(inout) :: table
    !
    type(list_iterator_t)    :: iter
    type(atom_t),    pointer :: atom
    type(species_t), pointer :: spec
    integer                  :: j, icnt, jcnt
    !
    print *, "this%natoms=list_len(list)"
    this%natoms=list_len(list)
    print *, "A!SSERT(this%natoms>0)", this%natoms
    ASSERT(this%natoms>0)
    print *, "S!AFE_ALLOCATE(this%atom(this%natoms))"
    SAFE_ALLOCATE(this%atom(this%natoms))
    print *, "this%space=>space"
    this%space=>space
    print *, "this%nspecies=table_len(table)"
    this%nspecies=table_len(table)
    print *, "A!SSERT(this%nspecies>0)", this%nspecies
    ASSERT(this%nspecies>0)
    print *, "S!AFE_ALLOCATE(this%species(this%nspecies))"
    SAFE_ALLOCATE(this%species(this%nspecies))
    print *, "icnt=0"
    icnt=0
    jcnt=0
    nullify(atom, spec)
    print *, "call list_iterator_init(iter, list)"
    call list_iterator_init(iter, list)
    do
      call list_iterator_next(iter, atom)
      if(.not.associated(atom))exit
      icnt=icnt+1
      call table_get(table, atom, spec)
      this%atom(icnt)=atom
      nullify(atom)
      j=species_index(spec)
      if(j==0)then
        jcnt=jcnt+1
        j=jcnt
        call species_set_index(spec, j)
        this%species(j)=spec
      end if
      nullify(spec)
      call atom_set_species(this%atom(icnt), this%species(j))
    end do
    call list_iterator_end(iter)
    ASSERT(icnt==this%natoms)
    ASSERT(jcnt==this%nspecies)
    this%ncatoms=0
    this%catom=>null()
    this%only_user_def=.false.
    this%species_time_dependent=.false.
    this%kinetic_energy=M_ZERO
    this%nlpp=.false.
    this%nlcc=.false.
    call distributed_nullify(this%atoms_dist, this%natoms)
    nullify(this%ionic_interaction_type)
    nullify(this%ionic_interaction_parameter)
    this%reduced_coordinates=.false.
    return
  end subroutine geometry_extend_finalize

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_extend_finalize)(this)
    type(TEMPLATE(geo_t)), intent(inout) :: this
    !
    call geometry_extend_finalize(this%geo, this%space, this%list, this%table)
    call TEMPLATE(geo_table_end)(this%table)
    call TEMPLATE(geo_list_end)(this%list)
    this%init=.true.
    return
  end subroutine TEMPLATE(geo_extend_finalize)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_extend)(this, that)
    type(TEMPLATE(geo_t)), intent(inout) :: this
    type(sub_t), optional, intent(in)    :: that
    !
    if(present(that))then
      call TEMPLATE(geo_extend_build)(this, that)
    else
      call TEMPLATE(geo_extend_finalize)(this)
    end if
    return
  end subroutine TEMPLATE(geo_extend)

#else
  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_init)(this, space, config)
    type(TEMPLATE(geo_t)),       intent(out) :: this
    type(space_t),       target, intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config
    !
    this%config=>config
    this%space=>space
    call geometry_init(this%geo, space, config)
    this%init=.true.
    return
  end subroutine TEMPLATE(geo_init)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_get)(this, that)
    type(TEMPLATE(geo_t)), target, intent(in) :: this
    type(geometry_t),     pointer             :: that
    !
    that=>this%geo
    return
  end subroutine TEMPLATE(geo_get)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_copy)(this, that)
    type(TEMPLATE(geo_t)), intent(out) :: this
    type(TEMPLATE(geo_t)), intent(in)  :: that
    !
    this%config=>that%config
    this%space=>that%space
    call geometry_copy(this%geo, that%geo)
#ifdef SUBTEMPLATE_NAME
    call list_copy(this%list, that%list)
    call table_copy(this%table, that%table)
#endif
    return
  end subroutine TEMPLATE(geo_copy)

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_list_end)(this)
    type(list_t), intent(inout) :: this
    !
    type(atom_t), pointer :: atom
    !
    nullify(atom)
    do
      call list_pop(this, atom)
      if(.not.associated(atom))exit
      call atom_end(atom)
      SAFE_DEALLOCATE_P(atom)
      nullify(atom)
    end do
    call list_end(this)
    return
  end subroutine TEMPLATE(geo_list_end)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_table_end)(this)
    type(table_t), intent(inout) :: this
    !
    type(atom_t),    pointer :: atom
    type(species_t), pointer :: spec
    !
    nullify(atom, spec)
    do
      call table_pop(this, atom, spec)
      if(.not.associated(spec))exit
      call species_end(spec)
      SAFE_DEALLOCATE_P(spec)
      nullify(atom, spec)
    end do
    call table_end(this)
    return
  end subroutine TEMPLATE(geo_table_end)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_end)(this)
    type(TEMPLATE(geo_t)), intent(inout) :: this
    !
#ifdef SUBTEMPLATE_NAME
    call TEMPLATE(geo_table_end)(this%table)
    call TEMPLATE(geo_list_end)(this%list)
#endif
    call geometry_end(this%geo)
    nullify(this%space, this%config)
    return
  end subroutine TEMPLATE(geo_end)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_iterator_init)(this, that)
    type(TEMPLATE(geo_iterator_t)), intent(out) :: this
    type(TEMPLATE(geo_t)),  target, intent(in)  :: that
    !
    print *, "this%n=0", associated(that%config), that%init
    this%n=0
    print *, "this%natoms=that%geo%natoms"
    this%natoms=that%geo%natoms
    print *, "this%nspecies=that%geo%nspecies"
    this%nspecies=that%geo%nspecies
    print *, "A!SSERT(associated(that%geo%atom))"
    ASSERT(associated(that%geo%atom))
    print *, "this%atom=>that%geo%atom"
    this%atom=>that%geo%atom
    print *, "A!SSERT(associated(that%geo%species))"
    ASSERT(associated(that%geo%species))
    print *, "this%species=>that%geo%species"
    this%species=>that%geo%species
    print *, "return"
    return
  end subroutine TEMPLATE(geo_iterator_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_iterator_next)(this, atom, species, ierr)
    type(TEMPLATE(geo_iterator_t)), intent(inout) :: this
    type(atom_t),                   intent(out)   :: atom
    type(species_t),                intent(out)   :: species
    integer,                        intent(out)   :: ierr
    !
    type(species_t), pointer :: spec
    !
    nullify(spec)
    this%n=this%n+1
    ierr=TEMPLATE(GEO_ITER_END)
    if(this%n<=this%natoms)then
      ierr=TEMPLATE(GEO_ITER_OK)
      atom=this%atom(this%n)
      call atom_get_species(atom, spec)
      species=spec
      nullify(spec)
      call species_set_index(species, 0)
      call atom_set_species(atom, species)
    end if
    return
  end subroutine TEMPLATE(geo_iterator_next)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_iterator_copy)(this, that)
    type(TEMPLATE(geo_iterator_t)), intent(out) :: this
    type(TEMPLATE(geo_iterator_t)), intent(in)  :: that
    !
    this%n=that%n
    this%natoms=that%natoms
    this%nspecies=that%nspecies
    this%atom=>that%atom
    this%species=>that%species
    return
  end subroutine TEMPLATE(geo_iterator_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(geo_iterator_end)(this)
    type(TEMPLATE(geo_iterator_t)), intent(inout) :: this
    !
    this%n=0
    this%natoms=0
    this%nspecies=0
    nullify(this%atom, this%species)
    return
  end subroutine TEMPLATE(geo_iterator_end)

end module TEMPLATE(geo_m)

!! Local Variables:
!! mode: f90
!! End:
