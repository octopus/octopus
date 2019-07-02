!! Copyright (C) 2019 H. Huebener, X. Andrade
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module scdm_oct_m
  use blacs_proc_grid_oct_m
  use blacs_oct_m
  use comm_oct_m
  use cube_oct_m
  use derivatives_oct_m
  use fft_oct_m
  use global_oct_m
  use index_oct_m
  use lalg_adv_oct_m
  use mesh_oct_m
  use mesh_cube_map_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use poisson_oct_m
  use poisson_fft_oct_m
  use profiling_oct_m
  use scalapack_oct_m
  use simul_box_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_parallel_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::               &
       scdm_t,            &
       scdm_init,         &
       dscdm_localize,    &
       zscdm_localize,    &
       scdm_rotate_states
  
  type scdm_t
    private
    type(states_elec_t)      :: st          !< localized orthogonal states
    FLOAT, pointer,   public :: center(:,:) !< coordinates of centers of states (in same units as mesh%x)
    integer,          public :: box_size    !< number of mesh points in the dimension of local box around scdm states 
                                            !! NOTE: this could be dynamic and state dependent
    integer,          public :: full_box    !< = (2*box_size+1)**3, i.e. number of points in box
    type(mesh_t)             :: boxmesh     !< mesh describing the small box
    type(cube_t)             :: boxcube     !< cube of the small box (used for fft in poisson solver
                                            !! has doubled size for truncation)
    integer, pointer, public :: box(:,:,:,:)  !< indices of global points that are contained in the local box for each state
    
    integer                  :: full_cube_n(3) !< dimension of cube of fullsimulation cell
    
    FLOAT, pointer,   public :: dpsi(:,:)   !< scdm states in their local box
    CMPLX, pointer,   public :: zpsi(:,:)   ! ^
    type(poisson_t),  public :: poisson     !< solver used to compute exchange with localized scdm states
    type(poisson_fft_t)      :: poisson_fft !< used for above poisson solver

    logical                  :: verbose     !< write info about SCDM procedure
    logical,          public :: psi_scdm    !< Hamiltonian is applied to an SCDM state

    integer                  :: nst         !< total number of states, copy os st%nst
    
    ! parallelization of scdm states
    type(mpi_grp_t)          :: st_grp      !< MPI group for states parallelization, inherited from st
    type(mpi_grp_t)          :: dom_grp     !< MPI group for domain parallelization, inherited from mesh
    type(mpi_grp_t),  public :: st_exx_grp  !< MPI group for state parallelization in the exchange operator
                                            !! this is a copy of the domain group, i.e. the domain group is
                                             !! used for states parallelization in the exchange operator
    integer,          public :: st_exx_start!< index of state distribution in the exchange operator
    integer,          public :: st_exx_end  !.
    logical                  :: root        !< this is a redundat flag equal to mpi_world%rank==0
#ifdef HAVE_SCALAPACK 
    type(blacs_proc_grid_t)  :: proc_grid  !< blacs context for RRQR on transpose states with scalapack
#endif
  end type scdm_t

  logical,public    :: scdm_is_init=.false.  ! is initialized
  logical,public    :: scdm_is_local=.false.  ! is localized

  type(profile_t), save :: prof_scdm, prof_scdm_QR, prof_scdm_matmul1, prof_scdm_matmul3
  
contains

!> this initializes the states and solver to compute exact exchange using the method described in
!! A. Damle, L. Lin, L. Ying: Compressed representation of Kohn-Sham orbitals via 
!!                            selected columns of the density matrix
!! http://arxiv.org/abs/1408.4926 (accepted in JCTC as of 17th March 2015)
subroutine scdm_init(st, namespace, der, fullcube, scdm, operate_on_scdm)
  type(states_elec_t), intent(in)  :: st !< this contains the KS set (for now from hm%hf_st which is confusing)
  type(namespace_t),   intent(in)  :: namespace
  type(derivatives_t)              :: der
  type(cube_t)                     :: fullcube !< cube of the full cell
  type(scdm_t)                     :: scdm
  logical,                optional :: operate_on_scdm  !< apply exchange to SCDM states by performing a basis rotation on the st object

  integer :: ii, jj, kk, ip
  logical :: operate_on_scdm_
  
  integer,  allocatable:: istart(:)
  integer,  allocatable:: iend(:)
  integer,  allocatable:: ilsize(:)
  integer :: box(3)
  FLOAT :: dummy
  FLOAT :: rcut ! orbital cutoff radius (box size) NOTE: this could be dynamic and state dependent
  
  PUSH_SUB(scdm_init)
  ! check if already initialized
  if (scdm_is_init) then
    POP_SUB(scdm_init)
    return
  end if
  
  if (st%d%nik > 1) call messages_not_implemented("SCDM with k-point sampling")
  if (der%mesh%sb%periodic_dim > 0 .and. der%mesh%sb%periodic_dim /= 3) &
       call messages_not_implemented("SCDM with mixed-periodicity")  
  
  ! determine whether we are applying the scdm exchange operator to scdm states
  ! NOTE: this should be always the case, but for now only in td
  ! set default
  if(present(operate_on_scdm)) then
    operate_on_scdm_ = operate_on_scdm
  else
    operate_on_scdm_ = .false.
  end if
  
  scdm%psi_scdm = operate_on_scdm_
  
  ! set mpi groups
  scdm%st_grp     = st%mpi_grp
  scdm%dom_grp    = der%mesh%mpi_grp
  scdm%st_exx_grp = der%mesh%mpi_grp ! this is used only when the exchange operator is applied

  scdm%root = (mpi_world%rank ==0)
  
  scdm%nst   = st%nst
  
  ! initialize state object for the SCDM states by copying
  call states_elec_copy(scdm%st,st)
  
  !%Variable SCDM_verbose
  !%Type logical
  !%Default no
  !%Section Hamiltonian
  !%Description
  !% Output detailed information on SCDM procedure.
  !%End
  call parse_variable(namespace, 'SCDM_verbose', .false., scdm%verbose)
  
  scdm%full_cube_n = fullcube%rs_n_global
  
  ! allocate centers
  SAFE_ALLOCATE(scdm%center(1:3,1:scdm%st%nst))
  
  ! make a cube around the center points
  ! with side length NOTE: this should be dynamic
  
  !%Variable SCDMCutoffRadius
  !%Type float
  !%Default 3. Ang
  !%Section Hamiltonian
  !%Description
  !% Controls the size of the box on which the SCDM states are defined (box size = 2*radius).
  !%End  
  call parse_variable(namespace, 'SCDMCutoffRadius', 3._8, rcut, units_inp%length)
  if (scdm%root.and.scdm%verbose) call messages_print_var_value(stdout, 'SCDM cutoff', rcut)
  ! box_size is half the size of the  box
  scdm%box_size = 0
  do ii = 1, 3
    scdm%box_size = max(scdm%box_size,ceiling(rcut/der%mesh%spacing(ii)))
  end do
  
  if (scdm%root .and. scdm%verbose) then
    call messages_print_var_value(stdout,'SCDM box_size', scdm%box_size)
    call messages_print_var_value(stdout,'SCDM box_size[Ang]', scdm%box_size*der%mesh%spacing(1)*0.529177249)
  end if
  scdm%full_box = (2*scdm%box_size+1)**3
  !check if scdm is not bigger than fft-grid of full simualtion cell  
  if (scdm%full_box > der%mesh%np_global) then
    message(1) = 'SCDM box larger than mesh, no point in using it'
    call messages_fatal(1,only_root_writes = .true.)
  end if
  dummy = 2*(2*scdm%box_size+1)*der%mesh%spacing(1)*0.529177249
  if (scdm%root .and. scdm%verbose) call messages_print_var_value(stdout, 'SCDM fullbox[Ang]', dummy)
  SAFE_ALLOCATE(scdm%box(1:scdm%box_size*2+1,1:scdm%box_size*2+1,1:scdm%box_size*2+1,1:scdm%st%nst))
  
  ! the localzied states defined in the box are distributed over state index for the exchange operator
  ! but using group of domain parallelization, here named st_exx_grp
  SAFE_ALLOCATE(istart(1:scdm%st_exx_grp%size))
  SAFE_ALLOCATE(iend(1:scdm%st_exx_grp%size))
  SAFE_ALLOCATE(ilsize(1:scdm%st_exx_grp%size))
  
  call multicomm_divide_range(st%nst, scdm%st_exx_grp%size, istart, iend, lsize=ilsize)
  scdm%st_exx_start = istart(scdm%st_exx_grp%rank+1)
  scdm%st_exx_end = iend(scdm%st_exx_grp%rank+1)
  
  ! allocate local boxes for the SCDM states
  if (.not.states_are_real(st)) then
    ! localized SCDM states 
    SAFE_ALLOCATE(scdm%zpsi(1:scdm%full_box,1:scdm%nst)) ! this can be distributed in memory 
  else ! real
    SAFE_ALLOCATE(scdm%dpsi(1:scdm%full_box,1:scdm%nst)) ! this can be distributed in memory  
  end if
  
  ! create a mesh object for the small box (for now each scdm state is in the same box, should be dynamic)
  ! only initialize values needed in the following (e.g. by poisson_fft_init)
  scdm%boxmesh%spacing(:) = minval(der%mesh%spacing(:))
  SAFE_ALLOCATE(scdm%boxmesh%sb)
  scdm%boxmesh%sb%periodic_dim = 0
  scdm%boxmesh%sb%dim = 3
  scdm%boxmesh%sb%klattice_primitive(1:3,1:3) = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
  scdm%boxmesh%sb%rlattice_primitive(1:3,1:3) = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
  
  !set mesh points 
  scdm%boxmesh%np = scdm%full_box
  scdm%boxmesh%np_global = scdm%boxmesh%np
  scdm%boxmesh%np_part = scdm%boxmesh%np_global
  scdm%boxmesh%np_part_global = scdm%boxmesh%np_global
  ! set index type of mesh
  scdm%boxmesh%idx%is_hypercube = .false.
  ! mesh has to be centered around zero with left overhang otherwise mesh_cub_map does not seem to work
  scdm%boxmesh%idx%nr(1,:) = -(scdm%box_size)
  scdm%boxmesh%idx%nr(2,:) =  (scdm%box_size) 
  
  scdm%boxmesh%idx%dim = 3
  scdm%boxmesh%idx%ll(:) = scdm%boxmesh%idx%nr(2,:) - scdm%boxmesh%idx%nr(1,:) + 1
  scdm%boxmesh%idx%enlarge(:) = 0
  SAFE_ALLOCATE(scdm%boxmesh%idx%lxyz(1:scdm%boxmesh%np,1:scdm%boxmesh%idx%dim))
  ! need to copy indices because otherwise line gets too long (precompiler?)
  ii = -(scdm%box_size)
  jj = (scdm%box_size)
  SAFE_ALLOCATE(scdm%boxmesh%idx%lxyz_inv(ii:jj,ii:jj,ii:jj))

  ip = 0
  do ii = scdm%boxmesh%idx%nr(1,1),scdm%boxmesh%idx%nr(2,1)
    do jj = scdm%boxmesh%idx%nr(1,2),scdm%boxmesh%idx%nr(2,2)
      do kk = scdm%boxmesh%idx%nr(1,3),scdm%boxmesh%idx%nr(2,3)
        ip = ip +1
        scdm%boxmesh%idx%lxyz(ip,1) = ii
        scdm%boxmesh%idx%lxyz(ip,2) = jj
        scdm%boxmesh%idx%lxyz(ip,3) = kk
        scdm%boxmesh%idx%lxyz_inv(ii,jj,kk) = ip
      end do
    end do
  end do
  
  scdm%boxmesh%parallel_in_domains = .false.
  
  call mesh_cube_map_init(scdm%boxmesh%cube_map, scdm%boxmesh%idx, scdm%boxmesh%np_global)
  
  ! create a cube object for the small box, with double size for coulomb truncation
  box(1:3) = scdm%boxmesh%idx%ll(1:3)*2
  call cube_init(scdm%boxcube, box, scdm%boxmesh%sb,fft_type=FFT_REAL, fft_library=FFTLIB_FFTW)
  
  ! set up poisson solver used for the exchange operator with scdm states
  ! this replicates poisson_kernel_init()
  scdm%poisson%poisson_soft_coulomb_param = M_ZERO
  if (der%mesh%sb%periodic_dim.eq.3) then
    call poisson_fft_init(scdm%poisson_fft, namespace, scdm%boxmesh, scdm%boxcube, &
         kernel=POISSON_FFT_KERNEL_HOCKNEY,fullcube=fullcube)
  else !non periodic case
    call poisson_fft_init(scdm%poisson_fft, namespace, scdm%boxmesh, scdm%boxcube, kernel=POISSON_FFT_KERNEL_SPH)
  end if
  
  ! create poisson object
  SAFE_ALLOCATE(scdm%poisson%der)
  SAFE_ALLOCATE(scdm%poisson%der%mesh)
  scdm%poisson%der%mesh = scdm%boxmesh
  scdm%poisson%der%mesh%vp%npart = 1
  scdm%poisson%method = POISSON_FFT
  scdm%poisson%kernel = POISSON_FFT_KERNEL_SPH
  scdm%poisson%cube = scdm%boxcube
  scdm%poisson%fft_solver = scdm%poisson_fft

#ifdef HAVE_SCALAPACK
  if(st%scalapack_compatible) then
    ! create a blacs context with the transpose row and col numbers
    scdm%proc_grid%npcol = st%dom_st_proc_grid%nprow
    scdm%proc_grid%nprow = st%dom_st_proc_grid%npcol
    CALL blacs_get( -1, 0, scdm%proc_grid%context )
    CALL blacs_gridinit(scdm%proc_grid%context, 'Row-major', scdm%proc_grid%nprow, scdm%proc_grid%npcol )
    CALL blacs_gridinfo(scdm%proc_grid%context,scdm%proc_grid%nprow,scdm%proc_grid%npcol,scdm%proc_grid%myrow,scdm%proc_grid%mycol)
 end if
#endif
  ! set flag to do this only once
  scdm_is_init = .true.
  
  call messages_write('done SCDM init')
  
  POP_SUB(scdm_init)
end subroutine scdm_init

!> wrapper routine to rotate  KS states into their SCDM representation
subroutine scdm_rotate_states(st,mesh,scdm)
  type(states_elec_t), intent(inout)  :: st
  type(mesh_t),        intent(in)     :: mesh
  type(scdm_t),        intent(inout)  :: scdm
  
  PUSH_SUB(scdm_rotate_states)
  
  if (.not.states_are_real(st)) then
    call zscdm_rotate_states(st,mesh,scdm)
  else
    call dscdm_rotate_states(st,mesh,scdm)
  end if
  
  POP_SUB(scdm_rotate_states)
  
end subroutine scdm_rotate_states
  
subroutine check_box_in_index(idx,center,size,out)
  type(index_t),    intent(in)  :: idx
  integer, intent(in)           :: center(:)
  integer, intent(in)           :: size
  logical, intent(out)          :: out(3)
  
  ! internal
  integer :: ix(3), corner(3,8), i1, idim
  
  PUSH_SUB(check_box_in_index)
  out(1:3) = .false.
  
  ! make the sign pattern for corners
  corner(:,1) = (/1,1,1/)
  corner(:,2) = (/1,1,-1/)
  corner(:,3) = (/1,-1,1/)
  corner(:,4) = (/-1,1,1/)
  corner(:,5) = (/1,-1,-1/)
  corner(:,6) = (/-1,1,-1/)
  corner(:,7) = (/-1,-1,1/)
  corner(:,8) = (/-1,-1,-1/)
  
  do idim=1,3
    do i1=1,8
      ix(:)=center(:) + size*corner(:,i1)
      if (ix(idim).lt.idx%nr(1,idim).or.ix(idim).gt.idx%nr(2,idim)) then
        out(idim) = .true. 
        exit
      end if
    end do
  end do
  POP_SUB(check_box_in_index)
end subroutine check_box_in_index

#include "undef.F90"
#include "real.F90"
#include "scdm_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "scdm_inc.F90"

end module scdm_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
