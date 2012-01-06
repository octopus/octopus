#include "global.h"

module fd_m
  use datasets_m, only: datasets_check
  use fd_basis_m, only: fd_basis_t, fd_basis_init, fd_basis_copy, fd_basis_end
  use fd_epot_m, only: fd_epot_t, fd_epot_init, fd_epot_vpsl, fd_epot_copy, fd_epot_end
  use fd_geometry_m, only: fd_atoms_add, fd_atoms_get, fd_atoms_get_classical
  use fd_states_m, only: fd_states_t, &
    fd_states_init, fd_states_density, fd_states_copy, fd_states_end
  use geometry_m, only: geometry_t, atom_classical_t, &
    geometry_init_from_dump, geometry_copy, geometry_end
  use global_m
  use io_m, only: io_open, io_close
  use loct_m, only: loct_pointer_copy
  use mesh_m, only: mesh_t
  use messages_m
  use parser_m, only: block_t, parse_block, parse_block_n, parse_block_cols, & 
    parse_block_string,parse_block_integer,  parse_block_float, parse_block_end
  use profiling_m
  use simul_box_m, only: simul_box_in_box
  use space_m, only: space_t, space_init_from_dump, space_copy, space_end
  use unit_m, only: units_to_atomic
  use unit_system_m, only: units_inp

  implicit none

  private
  public ::             &
    fd_t,               &
    fd_init,            &
    fd_used,            &
    fd_densities_init,  &
    fd_potentials_init, &
    fd_geometry_update, &
    fd_copy,            &
    fd_end

  type fd_sys_t
    type(space_t)     :: space
    type(fd_basis_t)  :: basis
    type(geometry_t)  :: geo
    type(fd_states_t) :: st
    type(fd_epot_t)   :: ep
  end type fd_sys_t

  type fd_t
    private
    integer                               :: nst=0
    type(fd_sys_t), pointer, dimension(:) :: sys=>null()
  end type fd_t

contains
  
  subroutine fd_init(this)
    type(fd_t), intent(inout) :: this
    !
    FLOAT, dimension(2*MAX_DIM) :: vtmp
    character(len=256)          :: dname
    type(block_t)               :: blk
    integer                     :: iunit, whpot, intrp
    integer                     :: gb, lin, col, ncls
    !
    !%Variable FrozenSystems
    !%Type block
    !%Section States
    !%Description
    !%End
    if(parse_block(datasets_check('FrozenSystems'), blk)==0) then
      this%nst=parse_block_n(blk)
      if(this%nst>0)then
        SAFE_ALLOCATE(this%sys(this%nst))
        do lin=1, this%nst
          ncls=parse_block_cols(blk, lin-1)
          ASSERT(ncls>2)
          call parse_block_string(blk, lin-1, 0, dname)
          print *,"*** fd_init dname: ", dname
          call parse_block_integer(blk, lin-1, 1, whpot)
          print *,"*** fd_init whpot: ", whpot
          call parse_block_integer(blk, lin-1, 2, intrp)
          print *,"*** fd_init intrp: ", intrp
          vtmp=0.0
          do col=3, ncls-1
            call parse_block_float(blk, lin-1, col, vtmp(col-2))
            print *,"*** fd_init vtmp(", col-2, "): ", vtmp(col-2)
          end do
          vtmp(:MAX_DIM)=units_to_atomic(units_inp%length, vtmp(:MAX_DIM))
          print *,"*** fd_init vtmp: ", vtmp
          iunit=io_open(trim(dname)//'/data.dmp', &
            action='read', is_tmp=.true., form='unformatted')
          read(iunit) gb
          ASSERT(gb==GUARD_BITS)
          read(iunit) gb
          ASSERT(gb==GUARD_BITS)
          call space_init_from_dump(this%sys(lin)%space, iunit)
          call geometry_init_from_dump(this%sys(lin)%geo, this%sys(lin)%space, iunit)
          call fd_basis_init(this%sys(lin)%basis, vtmp(:MAX_DIM), vtmp(MAX_DIM+1:), iunit)
          call fd_states_init(this%sys(lin)%st, &
            this%sys(lin)%basis, this%sys(lin)%geo, intrp, iunit)
          read(iunit) gb
          ASSERT(gb==GUARD_BITS)
          call fd_epot_init(this%sys(lin)%ep, &
            this%sys(lin)%basis, this%sys(lin)%geo, whpot, intrp, iunit)
          read(iunit) gb
          ASSERT(gb==GUARD_BITS)
          call io_close(iunit)
        end do
        call parse_block_end(blk)
      end if
    end if
    return
  end subroutine fd_init

  elemental function fd_used(this) result(usd)
    type(fd_t), intent(in) :: this
    !
    logical :: usd
    !
    usd=(this%nst>0)
    return
  end function fd_used

  subroutine fd_densities_init(this, mesh, geo, density)
    type(fd_t),            intent(in)    :: this
    type(mesh_t),          intent(in)    :: mesh
    type(geometry_t),      intent(in)    :: geo
    FLOAT, dimension(:,:), intent(inout) :: density
    !
    character(len=9) :: fname
    FLOAT, dimension(size(density,dim=2)) :: rho
    integer                               :: i, j, iunit
    !
    density=M_ZERO
    ASSERT(size(density,dim=1)==mesh%np_part_global)
    do i=1, this%nst
      do j=1, mesh%np_part_global
        if(simul_box_in_box(mesh%sb, geo, mesh%x(j,:)))then
          call fd_states_density(this%sys(i)%st, mesh%x(j,:), rho)
          density(j,:)=density(j,:)+rho
        end if
      end do
      write(fname,fmt="(a4,i1.1,a4)") "rho_", i, ".plt"
      iunit=io_open(fname, action="write", is_tmp=.true.)
      do j = 1, mesh%np_part_global
        if((mesh%x(j,2)==0.0).and.(mesh%x(j,3)==0.0))then
          write(iunit,*) mesh%x(j,1), density(j,:)
        end if
      end do
      call io_close(iunit)
    end do
    return
  end subroutine fd_densities_init

  subroutine fd_potentials_init(this, mesh, geo, vpsl)
    type(fd_t),          intent(in)    :: this
    type(mesh_t),        intent(in)    :: mesh
    type(geometry_t),    intent(in)    :: geo
    FLOAT, dimension(:), intent(inout) :: vpsl
    !
    FLOAT   :: pot
    integer :: i, j
    !
    vpsl=M_ZERO
    ASSERT(size(vpsl)==mesh%np_global)
    do i=1, this%nst
      do j=1, mesh%np_global
        if(simul_box_in_box(mesh%sb, geo, mesh%x(j,:)))then
          call fd_epot_vpsl(this%sys(i)%ep, mesh%x(j,:), pot)
          vpsl(j)=vpsl(j)+pot
        end if
      end do
    end do
    return
  end subroutine fd_potentials_init

  subroutine fd_geometry_update(this, geo)
    type(fd_t),       intent(in)    :: this
    type(geometry_t), intent(inout) :: geo
    !
    type(atom_classical_t), pointer, dimension(:) :: atms
    integer                                       :: i, natms
    !
    do i=1, this%nst
      nullify(atms)
      natms=this%sys(i)%geo%natoms
      call fd_atoms_get(this%sys(i)%geo, this%sys(i)%basis, atms)
      call fd_atoms_add(geo%nfatoms, geo%fatom, natms, atms)
      SAFE_DEALLOCATE_P(atms)
    end do
    do i=1, this%nst
      nullify(atms)
      natms=this%sys(i)%geo%ncatoms
      call fd_atoms_get_classical(this%sys(i)%geo, this%sys(i)%basis, atms)
      call fd_atoms_add(geo%ncatoms, geo%catom, natms, atms)
      SAFE_DEALLOCATE_P(atms)
    end do
    return
  end subroutine fd_geometry_update

  subroutine fd_copy(this_out, this_in)
    type(fd_t), intent(inout) :: this_out
    type(fd_t), intent(in)    :: this_in
    !
    integer :: i
    !
    this_out%nst=this_in%nst
    SAFE_ALLOCATE(this_out%sys(this_in%nst))
    do i = 1, this_in%nst
      call space_copy(this_out%sys(i)%space, this_in%sys(i)%space)
      call fd_basis_copy(this_out%sys(i)%basis, this_in%sys(i)%basis)
      call geometry_copy(this_out%sys(i)%geo, this_in%sys(i)%geo)
      call fd_states_copy(this_out%sys(i)%st, this_in%sys(i)%st)
      call fd_epot_copy(this_out%sys(i)%ep, this_in%sys(i)%ep)
    end do
    return
  end subroutine fd_copy

  subroutine fd_end(this)
    type(fd_t), intent(inout) :: this
    !
    integer :: i
    !
    do i = 1, this%nst
      call space_end(this%sys(i)%space)
      call fd_basis_end(this%sys(i)%basis)
      call geometry_end(this%sys(i)%geo)
      call fd_states_end(this%sys(i)%st)
      call fd_epot_end(this%sys(i)%ep)
    end do
    SAFE_DEALLOCATE_P(this%sys)
    this%nst=0
    return
  end subroutine fd_end

end module fd_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
