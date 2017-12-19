#include "global.h"

module dnst_oct_m

  use global_oct_m
  use intrpl_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m

  implicit none

  private

  public :: &
    dnst_t

  public ::      &
    dnst_init,   &
    dnst_start,  &
    dnst_acc,    &
    dnst_update, &
    dnst_reset,  &
    dnst_stop,   &
    dnst_set,    &
    dnst_get,    &
    dnst_copy,   &
    dnst_end

  integer, parameter :: default_nspin = 1

  type :: dnst_t
    private
    type(json_object_t),             pointer :: config =>null()
    type(simulation_t),              pointer :: sim    =>null()
    type(storage_t),                 pointer :: total  =>null()
    logical                                  :: xtrnl  = .false.
    integer                                  :: nspin  = default_nspin
    real(kind=wp), dimension(:), allocatable :: charge
    type(storage_t)                          :: data
  end type dnst_t

  interface dnst__acc__
    module procedure dnst__acc__type
    module procedure dnst__acc__cnfg
    module procedure dnst__acc__list
  end interface dnst__acc__

  interface dnst_init
    module procedure dnst_init_type
    module procedure dnst_init_copy
  end interface dnst_init

  interface dnst_set
    module procedure dnst_set_info
    module procedure dnst_set_charge
  end interface dnst_set

  interface dnst_get
    module procedure dnst_get_info
    module procedure dnst_get_charge
    module procedure dnst_get_config
    module procedure dnst_get_simulation
    module procedure dnst_get_storage
    module procedure dnst_get_density_1d
    module procedure dnst_get_density_2d
  end interface dnst_get

contains

  ! ---------------------------------------------------------
  subroutine dnst_init_type(this, config)
    type(dnst_t),        target, intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(dnst_init_type)

    nullify(cnfg)
    this%config => config
    call json_get(this%config, "external", this%xtrnl, ierr)
    if(ierr/=JSON_OK) this%xtrnl = .false.
    call json_get(this%config, "nspin", this%nspin, ierr)
    if(ierr/=JSON_OK) this%nspin = default_nspin
    ASSERT(this%nspin>0)
    ASSERT(this%nspin<3)
    SAFE_ALLOCATE(this%charge(1:this%nspin))
    call json_get(this%config, "charge", this%charge, ierr)
    if(ierr/=JSON_OK) this%charge = 0.0_wp
    call json_get(this%config, "storage", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(cnfg, "fine", .true.)
    call json_set(cnfg, "dimensions", this%nspin)
    call json_set(cnfg, "default", 0.0_wp)
    if(this%nspin>1)then
      SAFE_ALLOCATE(this%total)
      call storage_init(this%total, cnfg)
    else
      this%total => this%data
    end if
    call storage_init(this%data, cnfg)
    nullify(cnfg)

    POP_SUB(dnst_init_type)
  end subroutine dnst_init_type

  ! ---------------------------------------------------------
  subroutine dnst_init_copy(this, that, start)
    type(dnst_t),      intent(out) :: this
    type(dnst_t),      intent(in)  :: that
    logical, optional, intent(in)  :: start

    logical :: istr

    PUSH_SUB(dnst_init_copy)

    ASSERT(associated(that%config))
    call dnst_init(this, that%config)
    istr = .true.
    if(present(start)) istr = start
    if(istr)then
      if(present(start))then
        ASSERT(associated(that%sim))
      end if
      if(associated(that%sim)) call dnst_start(this, that%sim)
    end if

    POP_SUB(dnst_init_copy)
  end subroutine dnst_init_copy

  ! ---------------------------------------------------------
  subroutine dnst_start(this, sim)
    type(dnst_t),               intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(dnst_start)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    if(this%nspin>1) call storage_start(this%total, this%sim, ndim=1)
    call storage_start(this%data, this%sim)

    POP_SUB(dnst_start)
  end subroutine dnst_start

  ! ---------------------------------------------------------
  pure subroutine dnst_adjust_spin_1_n(this, that)
    real(kind=wp),               intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    select case(size(that))
    case(1)
      this = that(1)
    case(2)
      this = sum(that)
    case default
      this = -1.0_wp
    end select

  end subroutine dnst_adjust_spin_1_n

  ! ---------------------------------------------------------
  pure subroutine dnst_adjust_spin_2_n(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    select case(size(that))
    case(1)
      this = 0.5_wp*that(1)
    case(2)
      this = that
    case default
      this = -1.0_wp
    end select

  end subroutine dnst_adjust_spin_2_n

  ! ---------------------------------------------------------
  pure subroutine dnst_adjust_spin(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    select case(size(this))
    case(1)
      call dnst_adjust_spin_1_n(this(1), that)
    case(2)
      call dnst_adjust_spin_2_n(this, that)
    case default
      this = -1.0_wp
    end select

  end subroutine dnst_adjust_spin

  ! ---------------------------------------------------------
  subroutine dnst__load__(this, that, intp)
    type(dnst_t),   intent(inout) :: this
    type(dnst_t),   intent(in)    :: that
    type(intrpl_t), intent(in)    :: intp

    real(kind=wp), dimension(:), allocatable :: irho
    integer                                  :: ospn, ispn

    PUSH_SUB(dnst__load__)

    call dnst__reset__(this)
    call dnst_adjust_spin(this%charge, that%charge)
    call dnst_get(this, nspin=ospn)
    call dnst_get(that, nspin=ispn)
    if(ispn==ospn)then
      call storage_load(this%data, fnc1)
    else
      SAFE_ALLOCATE(irho(ispn))
      call storage_load(this%data, fnc2)
      SAFE_DEALLOCATE_A(irho)
    end if
    call dnst__update__(this)

    POP_SUB(dnst__load__)

  contains

    subroutine fnc1(x, f)
      real(kind=wp), dimension(:), intent(in)  :: x
      real(kind=wp), dimension(:), intent(out) :: f

      PUSH_SUB(dnst__load__.fnc1)

      call intrpl_eval(intp, x, f)
      
      POP_SUB(dnst__load__.fnc1)
    end subroutine fnc1

    subroutine fnc2(x, f)
      real(kind=wp), dimension(:), intent(in)  :: x
      real(kind=wp), dimension(:), intent(out) :: f

      integer :: ierr

      PUSH_SUB(dnst__load__.fnc2)

      f = 0.0_wp
      call intrpl_eval(intp, x, irho, ierr=ierr)
      if(ierr==INTRPL_OK) call dnst_adjust_spin(f, irho)
      
      POP_SUB(dnst__load__.fnc2)
    end subroutine fnc2

  end subroutine dnst__load__

  ! ---------------------------------------------------------
  subroutine dnst__acc__type(this, that)
    type(dnst_t), intent(inout) :: this
    type(dnst_t), intent(in)    :: that

    integer :: ispn

    PUSH_SUB(dnst__acc__type)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%nspin==that%nspin)
    do ispn = 1, this%nspin
      if(that%charge(ispn)>tiny(that%charge))then
        this%charge(ispn) = this%charge(ispn) + that%charge(ispn)
        call storage_add(this%data, that%data, ispn)
      end if
    end do

    POP_SUB(dnst__acc__type)
  end subroutine dnst__acc__type

  ! ---------------------------------------------------------
  subroutine dnst__acc__cnfg(this, that, type, config)
    type(dnst_t),        intent(inout) :: this
    type(dnst_t),        intent(in)    :: that
    integer,             intent(in)    :: type
    type(json_object_t), intent(in)    :: config

    type(dnst_t)   :: dnst
    type(intrpl_t) :: intp
    
    PUSH_SUB(dnst__acc__cnfg)

    ASSERT(.false.)
    
    call dnst_init(dnst, this, start=.true.)
    call dnst_set(dnst, static=.false., external=.false.)
    call intrpl_init(intp, that%data, type=type)
    call intrpl_attach(intp, config)
    call dnst__load__(dnst, that, intp)
    call intrpl_detach(intp)
    call dnst__acc__(this, dnst)
    call intrpl_end(intp)
    call dnst_end(dnst)

    POP_SUB(dnst__acc__cnfg)
  end subroutine dnst__acc__cnfg

  ! ---------------------------------------------------------
  subroutine dnst__acc__list(this, that, type, list)
    type(dnst_t),       intent(inout) :: this
    type(dnst_t),       intent(in)    :: that
    integer,            intent(in)    :: type
    type(json_array_t), intent(in)    :: list

    type(json_object_t), pointer :: cnfg
    type(json_array_iterator_t)  :: iter
    type(dnst_t)                 :: dnst
    type(intrpl_t)               :: intp
    integer                      :: ierr

    PUSH_SUB(dnst__acc__list)

    call dnst_init(dnst, this, start=.true.)
    call dnst_set(dnst, static=.false., external=.false.)
    call intrpl_init(intp, that%data, type=type)
    ASSERT(json_len(list)>0)
    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call intrpl_attach(intp, cnfg)
      call dnst__load__(dnst, that, intp)
      call intrpl_detach(intp)
      call dnst__acc__(this, dnst)
    end do
    call json_end(iter)
    nullify(cnfg)
    call intrpl_end(intp)
    call dnst_end(dnst)

    POP_SUB(dnst__acc__list)
  end subroutine dnst__acc__list

  ! ---------------------------------------------------------
  subroutine dnst_acc(this, that, config)
    type(dnst_t),                  intent(inout) :: this
    type(dnst_t),                  intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    type(json_array_t), pointer :: list
    integer                     :: type, ierr

    PUSH_SUB(dnst_acc)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    nullify(list)
    if(present(config))then
      call json_get(config, "interpolation", type, ierr)
      if(ierr/=JSON_OK) type = INTRPL_DEFAULT
      call json_get(config, "positions", list, ierr)
      if(ierr/=JSON_OK) nullify(list)
    end if
    if(associated(list))then
      call dnst__acc__(this, that, type, list)
    else
      call dnst__acc__(this, that)
    end if
    nullify(list)

    POP_SUB(dnst_acc)
  end subroutine dnst_acc

  ! ---------------------------------------------------------
  subroutine dnst__update__(this)
    type(dnst_t), intent(inout) :: this

    real(kind=wp) :: chrg
    integer       :: ispn

    PUSH_SUB(dnst__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_update(this%data)
    do ispn = 1, this%nspin
      chrg = this%charge(ispn)
      if(chrg>tiny(chrg))then
        call storage_norm(this%data, ispn, norm=chrg)
      else
        this%charge(ispn) = 0.0_wp
        call storage_reset(this%data, ispn)
      end if
    end do

    POP_SUB(dnst__update__)
  end subroutine dnst__update__

  ! ---------------------------------------------------------
  subroutine dnst_update(this)
    type(dnst_t), intent(inout) :: this

    logical :: static

    PUSH_SUB(dnst_update)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_get(this%data, lock=static)
    call storage_update(this%data)
    if(.not.(static.or.this%xtrnl)) call dnst__update__(this)
    if(this%nspin>1) call storage_reduce(this%total, this%data)

    POP_SUB(dnst_update)
  end subroutine dnst_update

  ! ---------------------------------------------------------
  subroutine dnst__reset__(this)
    type(dnst_t), intent(inout) :: this

    PUSH_SUB(dnst__reset__)

    ASSERT(associated(this%config))
    this%charge = 0.0_wp
    if(associated(this%sim)) call storage_reset(this%data)

    POP_SUB(dnst__reset__)
  end subroutine dnst__reset__

  ! ---------------------------------------------------------
  subroutine dnst_reset(this)
    type(dnst_t), intent(inout) :: this

    logical :: static

    PUSH_SUB(dnst_reset)

    ASSERT(associated(this%config))
    call storage_get(this%data, lock=static)
    if(.not.(static.or.this%xtrnl))then
      if(this%nspin>1)then
        if(associated(this%sim)) call storage_reset(this%total)
      end if
      call dnst__reset__(this)
    end if

    POP_SUB(dnst_reset)
  end subroutine dnst_reset

  ! ---------------------------------------------------------
  subroutine dnst_stop(this)
    type(dnst_t), intent(inout) :: this

    PUSH_SUB(dnst_stop)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    if(this%nspin>1) call storage_stop(this%total)
    call storage_stop(this%data)

    POP_SUB(dnst_stop)
  end subroutine dnst_stop

  ! ---------------------------------------------------------
  subroutine dnst_set_info(this, static, external)
    type(dnst_t),      intent(inout) :: this
    logical, optional, intent(in)    :: external
    logical, optional, intent(in)    :: static

    PUSH_SUB(dnst_set_info)

    ASSERT(associated(this%config))
    if(present(external)) this%xtrnl = external
    call storage_set(this%data, lock=static)
    if(this%nspin>1) call storage_set(this%total, lock=static)

    POP_SUB(dnst_set_info)
  end subroutine dnst_set_info

  ! ---------------------------------------------------------
  subroutine dnst_set_charge(this, charge, spin, total)
    type(dnst_t),      intent(inout) :: this
    real(kind=wp),     intent(in)    :: charge
    integer, optional, intent(in)    :: spin
    logical, optional, intent(in)    :: total

    logical :: ittl

    PUSH_SUB(dnst_set_charge)

    ASSERT(associated(this%config))
    ittl = .false.
    if(present(total)) ittl = total
    if(ittl)then
      ASSERT(.not.present(spin))
      this%charge = charge / real(this%nspin, kind=wp)
    else
      if(present(spin))then
        ASSERT(spin<=this%nspin)
        this%charge(spin) = charge
      else
        ASSERT(this%nspin==1)
        this%charge(1) = charge
      end if
    end if

    POP_SUB(dnst_set_charge)
  end subroutine dnst_set_charge

  ! ---------------------------------------------------------
  subroutine dnst_get_info(this, size, nspin, static, external, fine, use)
    type(dnst_t),      intent(in)  :: this
    integer, optional, intent(out) :: size
    integer, optional, intent(out) :: nspin
    logical, optional, intent(out) :: external
    logical, optional, intent(out) :: static
    logical, optional, intent(out) :: fine
    logical, optional, intent(out) :: use

    PUSH_SUB(dnst_get_info)

    ASSERT(associated(this%config))
    if(present(nspin)) nspin = this%nspin
    if(present(external)) external = this%xtrnl
    call storage_get(this%data, size=size, lock=static, fine=fine, alloc=use)

    POP_SUB(dnst_get_info)
  end subroutine dnst_get_info

  ! ---------------------------------------------------------
  subroutine dnst_get_charge(this, charge, spin, total)
    type(dnst_t),      intent(in)  :: this
    real(kind=wp),     intent(out) :: charge
    integer, optional, intent(in)  :: spin
    logical, optional, intent(in)  :: total

    logical :: ittl

    PUSH_SUB(dnst_get_charge)

    ASSERT(associated(this%config))
    ittl = .false.
    if(present(total)) ittl = total
    if(ittl)then
      ASSERT(.not.present(spin))
      charge = sum(this%charge)
    else
      if(present(spin))then
        ASSERT(spin<=this%nspin)
        charge = this%charge(spin)
      else
        ASSERT(this%nspin==1)
        charge = this%charge(1)
      end if
    end if

    POP_SUB(dnst_get_charge)
  end subroutine dnst_get_charge

  ! ---------------------------------------------------------
  subroutine dnst_get_config(this, that)
    type(dnst_t),         target, intent(in)  :: this
    type(json_object_t), pointer, intent(out) :: that

    PUSH_SUB(dnst_get_config)

    nullify(that)
    if(associated(this%config)) that=>this%config

    POP_SUB(dnst_get_config)
  end subroutine dnst_get_config

  ! ---------------------------------------------------------
  subroutine dnst_get_simulation(this, that)
    type(dnst_t),        target, intent(in)  :: this
    type(simulation_t), pointer, intent(out) :: that

    PUSH_SUB(dnst_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(dnst_get_simulation)
  end subroutine dnst_get_simulation

  ! ---------------------------------------------------------
  subroutine dnst_get_storage(this, that, total)
    type(dnst_t),     target, intent(in)  :: this
    type(storage_t), pointer, intent(out) :: that
    logical,        optional, intent(in)  :: total

    logical :: ittl

    PUSH_SUB(dnst_get_storage)

    nullify(that)
    ittl = .false.
    if(present(total)) ittl = total
    if(ittl)then
      that => this%total
    else
      that => this%data
    end if

    POP_SUB(dnst_get_storage)
  end subroutine dnst_get_storage

  ! ---------------------------------------------------------
  subroutine dnst_get_density_1d(this, that, spin, total)
    type(dnst_t),                         intent(in)  :: this
    real(kind=wp), dimension(:), pointer, intent(out) :: that
    integer,                    optional, intent(in)  :: spin
    logical,                    optional, intent(in)  :: total

    integer :: ispn
    logical :: ittl

    PUSH_SUB(dnst_get_density_1d)

    nullify(that)
    ispn = 1
    if(present(spin)) ispn = spin
    ittl = .false.
    if(present(total)) ittl = total
    if(ittl)then
      ASSERT(.not.present(spin))
      call storage_get(this%total, that)
    else
      call storage_get(this%data, that, ispn)
    end if

    POP_SUB(dnst_get_density_1d)
  end subroutine dnst_get_density_1d

  ! ---------------------------------------------------------
  subroutine dnst_get_density_2d(this, that, total)
    type(dnst_t),                           intent(in)  :: this
    real(kind=wp), dimension(:,:), pointer, intent(out) :: that
    logical,                      optional, intent(in)  :: total

    logical :: ittl

    PUSH_SUB(dnst_get_density_2d)

    nullify(that)
    ittl = .false.
    if(present(total)) ittl = total
    if(ittl)then
      call storage_get(this%total, that)
    else
      call storage_get(this%data, that)
    end if
    
    POP_SUB(dnst_get_density_2d)
  end subroutine dnst_get_density_2d

  ! ---------------------------------------------------------
  subroutine dnst_copy(this, that)
    type(dnst_t), intent(inout) :: this
    type(dnst_t), intent(in)    :: that

    PUSH_SUB(dnst_copy)

    call dnst_end(this)
    if(associated(that%config))then
      call dnst_init(this, that)
      this%charge(:) = that%charge(:)
      if(associated(that%sim)) then
        call storage_copy(this%data, that%data)
        call dnst_update(this)
      end if
    end if

    POP_SUB(dnst_copy)
  end subroutine dnst_copy

  ! ---------------------------------------------------------
  subroutine dnst_end(this)
    type(dnst_t), intent(inout) :: this

    PUSH_SUB(dnst_end)

    nullify(this%config, this%sim)
    if(this%nspin>1)then
      call storage_end(this%total)
      SAFE_DEALLOCATE_P(this%total)
    end if
    nullify(this%total)
    this%xtrnl = .false.
    this%nspin = 0
    SAFE_DEALLOCATE_A(this%charge)
    call storage_end(this%data)

    POP_SUB(dnst_end)
  end subroutine dnst_end

end module dnst_oct_m

!! Local Variables:
!! mode: f90
!! End:

