#include "global.h"

module frozen_density_m

  use base_density_m
  use basis_m
  use fio_density_m
  use global_m
  use json_m
  use kinds_m
  use mesh_m
  use messages_m
  use profiling_m
  use simulation_m
  use space_m

  implicit none

  private

  public ::                &
    frozen_density__acc__

contains

  ! ---------------------------------------------------------
  subroutine frozen_density__acc__charge(this, that)
    type(base_density_t), intent(inout) :: this !> frozen
    type(base_density_t), intent(in)    :: that !> fio

    real(kind=wp), dimension(:), allocatable :: ochrg, ichrg
    real(kind=wp)                            :: charge
    integer                                  :: ispin, nspin

    PUSH_SUB(frozen_density__acc__charge)

    call base_density_get(that, nspin=nspin)
    SAFE_ALLOCATE(ichrg(1:nspin))
    do ispin = 1, nspin
      call base_density_get(that, ichrg(ispin), ispin)
    end do
    call base_density_get(this, nspin=nspin)
    SAFE_ALLOCATE(ochrg(1:nspin))
    call frozen_density_adjust_spin(ochrg, ichrg)
    SAFE_DEALLOCATE_A(ichrg)
    do ispin = 1, nspin
      call base_density_get(this, charge, ispin)
      call base_density_set(this, (charge+ochrg(ispin)), ispin)
    end do
    SAFE_DEALLOCATE_A(ochrg)

    POP_SUB(frozen_density__acc__charge)
  end subroutine frozen_density__acc__charge

  ! ---------------------------------------------------------
  subroutine  frozen_density__acc__intrpl(this, intrpl, config)
    type(base_density_t),       intent(inout) :: this
    type(fio_density_intrpl_t), intent(in)    :: intrpl
    type(json_object_t),        intent(in)    :: config

    real(kind=wp), dimension(:,:),   pointer :: dnst
    real(kind=wp), dimension(:), allocatable :: x, orho, irho
    type(simulation_t),              pointer :: sim
    type(space_t),                   pointer :: space
    type(mesh_t),                    pointer :: mesh
    type(basis_t)                            :: basis
    integer                                  :: indx, jndx, np, nspin

    PUSH_SUB(frozen_density__acc__intrpl)

    nullify(dnst, sim, space, mesh)
    call base_density_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, space)
    ASSERT(associated(space))
    call basis_init(basis, space, config)
    call simulation_get(sim, mesh)
    ASSERT(associated(mesh))
    nullify(sim)
    call base_density_get(this, dnst)
    ASSERT(associated(dnst))
    SAFE_ALLOCATE(x(space%dim))
    call fio_density_get(intrpl, dim=nspin)
    ASSERT(nspin>0)
    ASSERT(nspin<3)
    SAFE_ALLOCATE(irho(nspin))
    call base_density_get(this, size=np, nspin=nspin)
    ASSERT(nspin>0)
    ASSERT(nspin<3)
    SAFE_ALLOCATE(orho(nspin))
    do indx = 1, np
      call basis_to_internal(basis, mesh%x(indx,1:space%dim), x)
      call fio_density_eval(intrpl, x, irho)
      call frozen_density_adjust_spin(orho, irho)
      do jndx = 1, nspin
        dnst(indx,jndx) = dnst(indx,jndx) + orho(jndx)
      end do
    end do
    SAFE_DEALLOCATE_A(irho)
    SAFE_DEALLOCATE_A(orho)
    SAFE_DEALLOCATE_A(x)
    nullify(dnst, space, mesh)
    call basis_end(basis)

    POP_SUB(frozen_density__acc__intrpl)
  end subroutine frozen_density__acc__intrpl

  ! ---------------------------------------------------------
  subroutine frozen_density__acc__(this, that, config)
    type(base_density_t), intent(inout) :: this !> frozen
    type(base_density_t), intent(in)    :: that !> fio
    type(json_object_t),  intent(in)    :: config

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    type(fio_density_intrpl_t)   :: intrp
    integer                      :: type, ierr
    logical                      :: acc

    PUSH_SUB(frozen_density__acc__)

    nullify(cnfg, list)
    call base_density_get(that, use=acc)
    if(acc)then
      call frozen_density__acc__charge(this, that)
      call json_get(config, "type", type, ierr)
      if(ierr==JSON_OK)then
        call fio_density_init(intrp, that, type)
      else
        call fio_density_init(intrp, that)
      end if
      call json_get(config, "positions", list, ierr)
      ASSERT(ierr==JSON_OK)
      ASSERT(json_len(list)>0)
      call json_init(iter, list)
      do
        nullify(cnfg)
        call json_next(iter, cnfg, ierr)
        if(ierr/=JSON_OK)exit
        call frozen_density__acc__intrpl(this, intrp, cnfg)
      end do
      call json_end(iter)
      nullify(cnfg, list)
      call fio_density_end(intrp)
      call base_density__update__(this)
    end if

    POP_SUB(frozen_density__acc__)

  end subroutine frozen_density__acc__

  ! ---------------------------------------------------------
  pure subroutine frozen_density_adjust_spin_1_n(this, that)
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

  end subroutine frozen_density_adjust_spin_1_n

  ! ---------------------------------------------------------
  pure subroutine frozen_density_adjust_spin_2_n(this, that)
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

  end subroutine frozen_density_adjust_spin_2_n

  ! ---------------------------------------------------------
  pure subroutine frozen_density_adjust_spin(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    select case(size(this))
    case(1)
      call frozen_density_adjust_spin_1_n(this(1), that)
    case(2)
      call frozen_density_adjust_spin_2_n(this, that)
    case default
      this = -1.0_wp
    end select

  end subroutine frozen_density_adjust_spin

end module frozen_density_m
 
!! Local Variables:
!! mode: f90
!! End:
