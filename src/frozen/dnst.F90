#include "global.h"

module dnst_oct_m

  use global_oct_m
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
    logical                                  :: static = .false.
    logical                                  :: xtrnl  = .false.
    logical                                  :: updt   = .false.
    integer                                  :: nspin  = 0
    real(kind=wp), dimension(:), allocatable :: charge
    type(storage_t)                          :: data
  end type dnst_t

  interface dnst_init
    module procedure dnst_init_type
    module procedure dnst_init_copy
  end interface dnst_init

  interface dnst_set
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
    call json_get(this%config, "static", this%static, ierr)
    if(ierr/=JSON_OK) this%static = .false.
    call json_get(this%config, "external", this%xtrnl, ierr)
    if(ierr/=JSON_OK) this%xtrnl = .false.
    this%updt = .false.
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
  subroutine dnst_init_copy(this, that)
    type(dnst_t), intent(out) :: this
    type(dnst_t), intent(in)  :: that

    PUSH_SUB(dnst_init_copy)

    ASSERT(associated(that%config))
    call dnst_init(this, that%config)
    if(associated(that%sim)) call dnst_start(this, that%sim)

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
  subroutine dnst_acc(this, that)
    type(dnst_t), intent(inout) :: this
    type(dnst_t), intent(in)    :: that

    integer :: ispn

    PUSH_SUB(dnst_acc)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%nspin==that%nspin)
    do ispn = 1, this%nspin
      if(that%charge(ispn)>0.0_wp)then
        this%charge(ispn) = this%charge(ispn) + that%charge(ispn)
        call storage_add(this%data, that%data, ispn)
        this%updt = .true.
      end if
    end do

    POP_SUB(dnst_acc)
  end subroutine dnst_acc

  ! ---------------------------------------------------------
  subroutine dnst_update(this)
    type(dnst_t), intent(inout) :: this

    real(kind=wp) :: chrg, intg
    integer       :: ispn

    PUSH_SUB(dnst_update)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(.not.this%static)then
      call storage_update(this%data)
      if(.not.this%xtrnl)then
        do ispn = 1, this%nspin
          chrg = this%charge(ispn)
          if(chrg>0.0_wp)then
            call storage_integrate(this%data, ispn, intg)
            ASSERT(.not.(intg<0.0_wp))
            if(intg>0.0_wp) call storage_mlt(this%data, ispn, chrg/intg)
          else
            this%charge(ispn) = 0.0_wp
            call storage_reset(this%data, ispn)
          end if
        end do
      end if
    end if
    if(.not.this%static.or.this%updt)then
      this%updt = .false.
      if(this%nspin>1)then
        call storage_reduce(this%total, this%data)
        if(.not.this%xtrnl)then
          chrg = sum(this%charge)
          if(chrg>0.0_wp)then
            call storage_integrate(this%total, intg)
            ASSERT(.not.(intg<0.0_wp))
            if(intg>0.0_wp) call storage_mlt(this%total, chrg/intg)
          else
            call storage_reset(this%total)
          end if
        end if
      end if
    end if

    POP_SUB(dnst_update)
  end subroutine dnst_update

  ! ---------------------------------------------------------
  subroutine dnst_reset(this)
    type(dnst_t), intent(inout) :: this

    PUSH_SUB(dnst_reset)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%charge = 0.0_wp
    if(this%nspin>1) call storage_reset(this%total)
    call storage_reset(this%data)

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
  subroutine dnst_get_info(this, size, nspin, fine, use)
    type(dnst_t),      intent(in)  :: this
    integer, optional, intent(out) :: size
    integer, optional, intent(out) :: nspin
    logical, optional, intent(out) :: fine
    logical, optional, intent(out) :: use

    PUSH_SUB(dnst_get_info)

    ASSERT(associated(this%config))
    if(present(nspin)) nspin = this%nspin
    call storage_get(this%data, size=size, fine=fine, alloc=use)

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
    type(dnst_t), target, intent(in) :: this
    type(json_object_t), pointer     :: that

    PUSH_SUB(dnst_get_config)

    nullify(that)
    if(associated(this%config)) that=>this%config

    POP_SUB(dnst_get_config)
  end subroutine dnst_get_config

  ! ---------------------------------------------------------
  subroutine dnst_get_simulation(this, that)
    type(dnst_t), target, intent(in) :: this
    type(simulation_t),  pointer     :: that

    PUSH_SUB(dnst_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(dnst_get_simulation)
  end subroutine dnst_get_simulation

  ! ---------------------------------------------------------
  subroutine dnst_get_storage(this, that, total)
    type(dnst_t),     target, intent(in) :: this
    type(storage_t), pointer             :: that
    logical,        optional, intent(in) :: total

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
    type(dnst_t),                 intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    integer,            optional, intent(in) :: spin
    logical,            optional, intent(in) :: total

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
    type(dnst_t),                   intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

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
      this%updt = that%updt
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
    this%static = .false.
    this%xtrnl = .false.
    this%updt = .false.
    this%nspin = 0
    SAFE_DEALLOCATE_A(this%charge)
    call storage_end(this%data)

    POP_SUB(dnst_end)
  end subroutine dnst_end

end module dnst_oct_m

!! Local Variables:
!! mode: f90
!! End:

