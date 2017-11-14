#include "global.h"

module fio_density_oct_m

  use base_density_oct_m
  use basis_oct_m
  use dnst_oct_m
  use global_oct_m
  use intrpl_oct_m
  use io_function_oct_m
  use json_oct_m
  use kinds_oct_m
  use mesh_oct_m
  use messages_oct_m
  use path_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use space_oct_m
  use storage_oct_m

  implicit none

  private

  public ::              &
    fio_density__init__, &
    fio_density__load__

  public ::               &
    fio_density_intrpl_t

  public ::                   &
    fio_density_intrpl_init,  &
    fio_density_intrpl_start, &
    fio_density_intrpl_eval,  &
    fio_density_intrpl_stop,  &
    fio_density_intrpl_get,   &
    fio_density_intrpl_copy,  &
    fio_density_intrpl_end

  type :: fio_density_intrpl_t
    private
    type(base_density_t), pointer :: self =>null()
    type(simulation_t),   pointer :: sim  =>null()
    type(json_object_t)           :: config
    type(base_density_t)          :: dnst
    type(intrpl_t)                :: intrp
  end type fio_density_intrpl_t

  interface fio_density_intrpl_get
    module procedure fio_density_intrpl_get_info
    module procedure fio_density_intrpl_get_type
  end interface fio_density_intrpl_get

contains

  ! ---------------------------------------------------------
  subroutine fio_density__init__(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(fio_density__init__)

    call base_density__init__(this, init)

    POP_SUB(fio_density__init__)

  contains
    
    subroutine init(this, config)
      type(dnst_t),        intent(out) :: this
      type(json_object_t), intent(in)  :: config
      
      real(kind=wp) :: chrb, chrt
      integer       :: ispn, nspin, ierr
      logical       :: ipol

      PUSH_SUB(fio_density__init__.init)

      call dnst_init(this, config)
      call dnst_get(this, nspin=nspin)
      call json_get(config, "invertpolarization", ipol, ierr)
      if(ierr/=JSON_OK) ipol = .false.
      if(ipol.and.(nspin>1))then
        do ispn = 1, nspin/2
          call dnst_get(this, chrb, spin=ispn)
          call dnst_get(this, chrt, spin=nspin-ispn+1)
          call dnst_set(this, chrt, spin=ispn)
          call dnst_set(this, chrb, spin=nspin-ispn+1)
        end do
      end if

      POP_SUB(fio_density__init__.init)
    end subroutine init
    
  end subroutine fio_density__init__

  ! ---------------------------------------------------------
  subroutine fio_density__read__(this, dir, file, ispin)
    type(base_density_t), intent(inout) :: this
    character(len=*),     intent(in)    :: dir
    character(len=*),     intent(in)    :: file
    integer,              intent(in)    :: ispin

    real(kind=wp), dimension(:,:), pointer :: dnst
    type(simulation_t),            pointer :: sim
    type(mesh_t),                  pointer :: mesh
    character(len=MAX_PATH_LEN)            :: fpth
    integer                                :: ierr

    PUSH_SUB(fio_density__read__)

    nullify(dnst, sim, mesh)
    call path_join(dir, file, fpth)
    call base_density_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, mesh, fine=.true.)
    ASSERT(associated(mesh))
    nullify(sim)
    call base_density_get(this, dnst)
    ASSERT(associated(dnst))
    call dio_function_input(fpth, mesh, dnst(:,ispin), ierr)
    nullify(dnst, mesh)
    if(ierr/=0)then
      call base_density_end(this)
      message(1) = "Could not read the input file: '"//trim(adjustl(file))//".obf'"
      message(2) = "in the directory: '"//trim(adjustl(dir))//"'"
      write(unit=message(3), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(3)
    end if

    POP_SUB(fio_density__read__)
  end subroutine fio_density__read__

  ! ---------------------------------------------------------
  subroutine fio_density__load__(this)
    type(base_density_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: ispn, nspin, ierr
    logical                      :: load, ipol

    PUSH_SUB(fio_density__load__)

    nullify(cnfg, list)
    call base_density_get(this, use=load)
    if(load)then
      call base_density_get(this, cnfg)
      ASSERT(associated(cnfg))
      call json_get(cnfg, "dir", dir, ierr)
      ASSERT(ierr==JSON_OK)
      call json_get(cnfg, "files", list, ierr)
      ASSERT(ierr==JSON_OK)
      call base_density_get(this, nspin=nspin)
      ASSERT(json_len(list)==nspin)
      call json_get(cnfg, "invertpolarization", ipol, ierr)
      if(ierr/=JSON_OK) ipol = .false.
      do ispn = 1, nspin
        call json_get(list, ispn, file, ierr)
        ASSERT(ierr==JSON_OK)
        if(ipol)then
          call fio_density__read__(this, trim(adjustl(dir)), trim(adjustl(file)), nspin-ispn+1)
        else
          call fio_density__read__(this, trim(adjustl(dir)), trim(adjustl(file)), ispn)
        end if
      end do
      nullify(cnfg, list)
      call base_density__update__(this)
    end if
    !call base_density_set(this, static=.true.)

    POP_SUB(fio_density__load__)
  end subroutine fio_density__load__

  ! ---------------------------------------------------------
  pure subroutine fio_density_adjust_spin_1_n(this, that)
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

  end subroutine fio_density_adjust_spin_1_n

  ! ---------------------------------------------------------
  pure subroutine fio_density_adjust_spin_2_n(this, that)
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

  end subroutine fio_density_adjust_spin_2_n

  ! ---------------------------------------------------------
  pure subroutine fio_density_adjust_spin(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    select case(size(this))
    case(1)
      call fio_density_adjust_spin_1_n(this(1), that)
    case(2)
      call fio_density_adjust_spin_2_n(this, that)
    case default
      this = -1.0_wp
    end select

  end subroutine fio_density_adjust_spin

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_init(this, that, type)
    type(fio_density_intrpl_t),   intent(out) :: this
    type(base_density_t), target, intent(in)  :: that
    integer,            optional, intent(in)  :: type

    type(json_object_t), pointer :: cnfg
    type(storage_t),     pointer :: data

    PUSH_SUB(fio_density_intrpl_init)

    nullify(cnfg, data)
    this%self => that
    nullify(this%sim)
    call base_density_get(this%self, cnfg)
    ASSERT(associated(cnfg))
    call json_copy(this%config, cnfg)
    nullify(cnfg)
    call json_set(this%config, "reduce", .false.)
    call json_set(this%config, "external", .false.)
    call json_set(this%config, "default", .true.)
    call base_density_get(that, data)
    ASSERT(associated(data))
    call intrpl_init(this%intrp, data, type=type)
    nullify(data)

    POP_SUB(fio_density_intrpl_init)
  end subroutine fio_density_intrpl_init

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_start(this, sim, nspin)
    type(fio_density_intrpl_t), intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    integer,                    intent(in)    :: nspin

    PUSH_SUB(fio_density_intrpl_start)

    ASSERT(associated(this%self))
    ASSERT(.not.associated(this%sim))
    ASSERT(nspin>0)
    ASSERT(nspin<3)
    this%sim => sim
    call json_set(this%config, "nspin", nspin)
    call base_density_init(this%dnst, this%config)
    call base_density_start(this%dnst, this%sim)
    call base_density_set(this%dnst, static=.true.)

    POP_SUB(fio_density_intrpl_start)
  end subroutine fio_density_intrpl_start

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl__charge__(this)
    type(fio_density_intrpl_t),  intent(inout) :: this

    real(kind=wp), dimension(:), allocatable :: ichr, ochr
    real(kind=wp)                            :: chrg
    integer                                  :: isp, ispn, ospn

    PUSH_SUB(fio_density_intrpl__charge__)

    if(associated(this%self).and.associated(this%sim))then
      call base_density_get(this%self, nspin=ispn)
      ASSERT(ispn>0)
      ASSERT(ispn<3)
      call base_density_get(this%dnst, nspin=ospn)
      ASSERT(ospn>0)
      ASSERT(ospn<3)
      if(ispn/=ospn)then
        SAFE_ALLOCATE(ichr(ispn))
        do isp = 1, ispn
          call base_density_get(this%self, ichr(isp), spin=isp)
        end do
        SAFE_ALLOCATE(ochr(ospn))
        call fio_density_adjust_spin(ochr, ichr)
        SAFE_DEALLOCATE_A(ichr)
        do isp = 1, ospn
          call base_density_set(this%dnst, ochr(isp), spin=isp)
        end do
        SAFE_DEALLOCATE_A(ochr)
      else
        do isp = 1, ispn
          call base_density_get(this%self, chrg, spin=isp)
          call base_density_set(this%dnst, chrg, spin=isp)
        end do
      end if
    end if

    POP_SUB(fio_density_intrpl__charge__)
  end subroutine fio_density_intrpl__charge__

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl__eval__(this, x, val)
    type(fio_density_intrpl_t),  intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: val

    integer :: ierr

    PUSH_SUB(fio_density_intrpl__eval__)

    ierr = INTRPL_NI
    if(associated(this%self)) call intrpl_eval(this%intrp, x, val, ierr)
    if(ierr/=INTRPL_OK) val = 0.0_wp

    POP_SUB(fio_density_intrpl__eval__)
  end subroutine fio_density_intrpl__eval__

  ! ---------------------------------------------------------
  subroutine  fio_density_intrpl_eval(this, that, config)
    type(fio_density_intrpl_t), target, intent(inout) :: this
    type(base_density_t),              pointer        :: that
    type(json_object_t),                intent(in)    :: config

    real(kind=wp), dimension(:,:),   pointer :: dnst
    real(kind=wp), dimension(:), allocatable :: x, irho
    type(space_t),                   pointer :: space
    type(mesh_t),                    pointer :: mesh
    type(basis_t)                            :: basis
    integer                                  :: indx, np, nspin

    PUSH_SUB(fio_density_intrpl_eval)

    nullify(that, dnst, space, mesh)
    call simulation_get(this%sim, space)
    ASSERT(associated(space))
    ASSERT(space%dim>0)
    SAFE_ALLOCATE(x(space%dim))
    call basis_init(basis, space, config)
    call simulation_get(this%sim, mesh)
    ASSERT(associated(mesh))
    call fio_density_intrpl_get(this, dim=nspin)
    ASSERT(nspin>0)
    ASSERT(nspin<3)
    SAFE_ALLOCATE(irho(nspin))
    that => this%dnst
    call base_density_set(this%dnst, static=.false.)
    call base_density_get(this%dnst, size=np)
    call base_density_get(this%dnst, dnst)
    ASSERT(associated(dnst))
    call base_density__reset__(this%dnst)
    call fio_density_intrpl__charge__(this)
    do indx = 1, np
      call basis_to_internal(basis, mesh%x(indx,1:space%dim), x)
      call fio_density_intrpl__eval__(this, x, irho)
      call fio_density_adjust_spin(dnst(indx,:), irho)
    end do
    call base_density__update__(this%dnst)
    call base_density_set(this%dnst, static=.true.)
    SAFE_DEALLOCATE_A(irho)
    SAFE_DEALLOCATE_A(x)
    nullify(dnst, space, mesh)
    call basis_end(basis)

    POP_SUB(fio_density_intrpl_eval)
  end subroutine fio_density_intrpl_eval

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_stop(this)
    type(fio_density_intrpl_t), intent(inout) :: this

    PUSH_SUB(fio_density_intrpl_stop)

    ASSERT(associated(this%self))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call base_density_end(this%dnst)

    POP_SUB(fio_density_intrpl_stop)
  end subroutine fio_density_intrpl_stop

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_get_info(this, type, dim, default)
    type(fio_density_intrpl_t), intent(in)  :: this
    integer,          optional, intent(out) :: type
    integer,          optional, intent(out) :: dim
    real(kind=wp),    optional, intent(out) :: default

    PUSH_SUB(fio_density_intrpl_get_info)

    call intrpl_get(this%intrp, type, dim, default)

    POP_SUB(fio_density_intrpl_get_info)
  end subroutine fio_density_intrpl_get_info

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_get_type(this, that)
    type(fio_density_intrpl_t), intent(in) :: this
    type(base_density_t),      pointer     :: that

    PUSH_SUB(fio_density_intrpl_get_type)

    nullify(that)
    if(associated(this%self)) that => this%self

    POP_SUB(fio_density_intrpl_get_type)
  end subroutine fio_density_intrpl_get_type

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_copy(this, that)
    type(fio_density_intrpl_t), intent(out) :: this
    type(fio_density_intrpl_t), intent(in)  :: that

    PUSH_SUB(fio_density_intrpl_copy)

    call fio_density_intrpl_end(this)
    if(associated(that%self))then
      this%self => that%self
      this%sim => that%sim
      call json_copy(this%config, that%config)
      if(associated(that%sim)) call base_density_copy(this%dnst, that%dnst)
      call intrpl_copy(this%intrp, that%intrp)
    end if

    POP_SUB(fio_density_intrpl_copy)
  end subroutine fio_density_intrpl_copy

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_end(this)
    type(fio_density_intrpl_t), intent(inout) :: this

    PUSH_SUB(fio_density_intrpl_end)

    if(associated(this%self))then
      if(associated(this%sim)) call base_density_end(this%dnst)
      nullify(this%sim)
      call json_end(this%config)
      call intrpl_end(this%intrp)
    end if
    nullify(this%self)

    POP_SUB(fio_density_intrpl_end)
  end subroutine fio_density_intrpl_end

end module fio_density_oct_m

!! Local Variables:
!! mode: f90
!! End:
