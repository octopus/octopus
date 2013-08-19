#include "global.h"

#define TEMPLATE_NAME frozen
#define SUBTEMPLATE_NAME fio
#define EXTERNAL
#define IONIC
#define TNADD
#include "thamiltonian.F90"
#include "tcalc.F90"
#include "twrap.F90"
#undef TNADD
#undef IONIC
#undef EXTERNAL
#undef SUBTEMPLATE_NAME
#undef TEMPLATE_NAME

module frozen_config_m

  use global_m
  use messages_m
  use profiling_m

  use datasets_m,    only: datasets_check
  use fio_config_m,  only: fio_config_parse
  use functional_m,  only: XC_NONE
  use json_m,        only: json_object_t, json_array_t, json_init, json_set, json_append
  use kinds_m,       only: wp
  use parser_m,      only: block_t, parse_block, parse_block_end, parse_block_n, parse_block_cols, parse_integer, parse_float
  use unit_system_m, only: units_inp
  use XC_F90(lib_m), only: XC_F90(family_from_id)

  implicit none

  private
  public ::              &
    frozen_config_parse

contains

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_list(this, ndim)
    type(json_array_t), intent(out) :: this
    integer,            intent(in)  :: ndim
    !
    type(json_object_t), pointer :: cnfg
    type(block_t)                :: block
    integer                      :: lin, nlin, ncls
    !
    !%Variable FrozenSystems
    !%Type block
    !%Section System
    !%Description
    !% Lists the directories and the optional parameters to use on the reading of the frozen systems.
    !% the parameters are: the type of interpolation to use, wether to use or not the frozen 
    !% potential and the translation and rotation to apply to the system.
    !%
    !% <tt>%FrozenSystems
    !% <br>&nbsp;&nbsp;'directory_1' | nearest | yes | x | y | z | theta_xy | theta_xz | theta_yz
    !% <br>&nbsp;&nbsp;'directory_2' | nearest |  no | x | y | z | theta_xy
    !% <br>&nbsp;&nbsp;'directory_2' | nearest | yes | x 
    !% <br>%</tt>
    !%
    !%End
    !
    nullify(cnfg)
    call json_init(this)
    if(parse_block(datasets_check('FrozenSystems'),block)==0) then
      nlin=parse_block_n(block)
      if(nlin>0)then
        do lin=1, nlin
          ncls=parse_block_cols(block, lin-1)
          ASSERT(ncls>0)
          SAFE_ALLOCATE(cnfg)
          call fio_config_parse(cnfg, block, ndim, lin, ncls)
          call json_append(this, cnfg)
          nullify(cnfg)
        end do
      end if
      call parse_block_end(block)
    end if
    return
  end subroutine frozen_config_parse_list

  ! ---------------------------------------------------------
  subroutine frozen_space_config_parse(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim
    !
    call json_init(this)
    call json_set(this, "dimensions", ndim)
    return
  end subroutine frozen_space_config_parse

  ! ---------------------------------------------------------
  subroutine frozen_density_config_parse(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: nspin
    !
    call json_init(this)
    call json_set(this, "SpinComponents", nspin)
    return
  end subroutine frozen_density_config_parse

  ! ---------------------------------------------------------
  subroutine frozen_states_config_parse(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call frozen_density_config_parse(cnfg, nspin)
    call json_set(this, "density", cnfg)
    nullify(cnfg)
    return
  end subroutine frozen_states_config_parse

  ! ---------------------------------------------------------
  subroutine frozen_system_config_parse(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call frozen_space_config_parse(cnfg, ndim)
    call json_set(this, "space", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call frozen_states_config_parse(cnfg, nspin)
    call json_set(this, "states", cnfg)
    nullify(cnfg)
    return
  end subroutine frozen_system_config_parse

  ! ---------------------------------------------------------
  subroutine frozen_simulation_config_parse(this)
    type(json_object_t), intent(out) :: this
    !
    type(json_object_t), pointer :: cnfg
    !
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call json_init(cnfg)
    call json_set(this, "grid", cnfg)
    nullify(cnfg)
    return
  end subroutine frozen_simulation_config_parse

!!$  ! ---------------------------------------------------------
!!$  subroutine functional_parse_config(this, nel)
!!$    type(json_object_t),     intent(out) :: this
!!$    real(kind=wp), optional, intent(in)  :: nel
!!$    !
!!$    real(kind=wp) :: rtmp
!!$    integer       :: itmp
!!$    !
!!$    call json_init(this)
!!$    if(present(nel))call json_set(this, "nel", nel)
!!$    call parse_float(datasets_check('Xalpha'), 1.0_wp, rtmp)
!!$    call json_set(this, "Xalpha", rtmp)
!!$    call parse_integer(datasets_check('Interaction1D'), INT_SOFT_COULOMB, itmp)
!!$    call json_set(this, "Interaction1D", itmp)
!!$    call parse_float(datasets_check("Interaction1DScreening"), 1.0_wp, rtmp)
!!$    call json_set(this, "Interaction1DScreening", rtmp)
!!$    call parse_integer(datasets_check("LB94_modified"), 0, itmp)
!!$    call json_set(this, "LB94_modified", itmp)
!!$    call parse_float(datasets_check('LB94_threshold'), 1.0e-6_wp, rtmp)
!!$    call json_set(this, "LB94_threshold", rtmp)
!!$    return
!!$  end subroutine functional_parse_config
  
  ! ---------------------------------------------------------
  subroutine frozen_tnadd_config_parse(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: nspin
    !
    real(kind=wp) :: factor, energy
    integer       :: id, family
    !
    !%Variable TnaddFunctional
    !%Type integer
    !%Section Hamiltonian
    !%Description
    !% Chooses the Kinetic Functional to use and passes any additional parameters.
    !%Variable TnaddFactor
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian
    !%Description
    !%
    !%End
    call json_init(this)
    call json_set(this, "SpinComponents", nspin)
    call parse_integer(datasets_check('TnaddFunctional'), XC_NONE, id)
    call json_set(this, "functional", id)
    if(id>XC_NONE)then
      family=XC_F90(family_from_id)(id)
      call json_set(this, "family", family)
      call parse_float(datasets_check('TnaddFactor'), 1.0_wp, factor)
      if(abs(factor)<1.0e-16_wp)then
        message(1)="The 'TnaddFactor' value specified may be too small."
        call messages_warning(1)
      end if
      call json_set(this, "factor", factor)
    end if
    return
  end subroutine frozen_tnadd_config_parse

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_config_parse(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    !
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call json_init(cnfg)
    call json_set(this, "external", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call json_init(cnfg)
    call json_set(this, "ionic", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call frozen_tnadd_config_parse(cnfg, nspin)
    call json_set(this, "tnadd", cnfg)
    nullify(cnfg)
    return
  end subroutine frozen_hamiltonian_config_parse

  ! ---------------------------------------------------------
  subroutine frozen_calc_config_parse(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call frozen_system_config_parse(cnfg, ndim, nspin)
    call json_set(this, "system", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call frozen_simulation_config_parse(cnfg)
    call json_set(this, "simulation", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call frozen_hamiltonian_config_parse(cnfg, nspin)
    call json_set(this, "hamiltonian", cnfg)
    nullify(cnfg)
    return
  end subroutine frozen_calc_config_parse

  ! ---------------------------------------------------------
  subroutine frozen_config_parse(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    !
    nullify(cnfg)
    call json_init(this)
    call json_set(this, "temporary", .false.)
    call json_set(this, "pass_simulation", .true.)
    SAFE_ALLOCATE(cnfg)
    call frozen_calc_config_parse(cnfg, ndim, nspin)
    call json_set(this, "calculation", cnfg)
    nullify(cnfg, list)
    SAFE_ALLOCATE(list)
    call frozen_config_parse_list(list, ndim)
    call json_set(this, "list", list)
    nullify(list)
    return
  end subroutine frozen_config_parse

end module frozen_config_m

!! Local Variables:
!! mode: f90
!! End:
