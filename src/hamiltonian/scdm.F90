!! Copyright (C) 2015 H. Huebener
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
!! $Id$

#include "global.h"

module scdm_m
  use batch_m
  use batch_ops_m
  use blas_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use cmplxscl_m
  use cube_m
  use cube_function_m
  use derivatives_m
  use fft_m
  use nfft_m
  use geometry_m
  use global_m
  use grid_m
  use hardware_m
  use index_m
  use io_m
  use io_function_m
  use kpoints_m
  use lalg_basic_m
  use math_m
  use mesh_m
  use mesh_cube_map_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use multicomm_m
  use opencl_m
  use opencl_m
  use par_vec_m
  use parser_m
  use poisson_m
  use poisson_fft_m
  use profiling_m
  use simul_box_m
  use smear_m
  use states_m
  use states_calc_m
  use states_dim_m
  use types_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use xc_m
  use XC_F90(lib_m)

  implicit none

  private
  public ::               &
       scdm_t,            &
       scdm_init,         &
       dscdm_localize,    &
       zscdm_localize,    &
       scdm_rotate_states
  
  type scdm_t
    type(states_t)   :: st          !< localized orthogonal states
    type(poisson_t)  :: poisson1    !< solver (not used, only for testing)
    type(cube_t)     :: cube        !< mesh cube for fft
    FLOAT, pointer   :: center(:,:) !< coordinates of centers of states (in same units as mesh%x)
    FLOAT            :: rcut        !< orbital cutoff radius (box size) NOTE: this could be dynamic and state dependent
    integer          :: box_size    !< number of mesh points in the dimension of local box around scdm states 
                                    !! NOTE: this could be dynamic and state dependent
    integer          :: full_box    !< = (2*box_size+1)**3, i.e. number of points in box
    type(mesh_t)     :: boxmesh     !< mesh describing the small box
    type(cube_t)     :: boxcube     !< cube of the small box (used for fft in poisson solver
                                    !! has doubled size for truncation)
    integer, pointer :: box(:,:,:,:)  !< indices of global points that are contained in the local box for each state
    
    integer          :: full_cube_n(3) !< dimension of cube of fullsimulation cell
    
    FLOAT, pointer   :: dpsi(:,:)   !< scdm states in their local box
    CMPLX, pointer   :: zpsi(:,:)   ! ^
    type(poisson_t)  :: poisson     !< solver used to compute exchange with localized scdm states
    type(poisson_fft_t) :: poisson_fft !< used for above poisson solver
    type(cmplxscl_t)    :: cmplxscl
    logical, pointer :: periodic(:) !< tracks which scdm states are split by the periodic boundary conditions

    logical          :: re_ortho_normalize=.false. !< orthonormalize the scdm states
    logical          :: verbose     !< write info about SCDM procedure
    logical          :: psi_scdm    !< Hamiltonian is applied to an SCDM state
    
    ! parallelization of scdm states
    logical          :: root        !< this is a redundat flag equal to mesh%vp%rank==0
    integer          :: nst
    integer          :: st_start    !< the distributed index
    integer          :: st_end      ! .
    integer          :: lnst        ! .

  end type scdm_t

  logical,public    :: scdm_is_init=.false.  ! is initialized
  logical,public    :: scdm_is_local=.false.  ! is localized
  
  !> debug stuff
  type(geometry_t), public   :: scdm_geo

contains

  !> this initializes the states and solver to compute exact exchange using the method described in
  !! A. Damle, L. Lin, L. Ying: Compressed representation of Kohn-Sham orbitals via 
  !!                            selected columns of the density matrix
  !! http://arxiv.org/abs/1408.4926 (accepted in JCTC as of 17th March 2015)
  subroutine scdm_init(st,der,fullcube,scdm)
    
    type(states_t), intent(in)  :: st !< this contains the KS set (for now from hm%hf_st which is confusing)
    type(derivatives_t) :: der
    type(cube_t) :: fullcube !< cube of the full cell
    type(scdm_t) :: scdm
    
    type(cmplxscl_t) :: cmplxscl
    integer :: ii, jj, kk, ip, rank
    integer :: inp_calc_mode
    !debug
    integer :: temp(3)
    
    integer,  allocatable:: istart(:)
    integer,  allocatable:: iend(:)
    integer,  allocatable:: ilsize(:)
    integer :: box(3)
    FLOAT :: dummy, enlarge(3)
    
    PUSH_SUB(scdm_init)
    ! check if already initialized
    if (scdm_is_init) then
      POP_SUB(scdm_init)
      return
    end if
    
    if (st%lnst /= st%nst) call messages_not_implemented("SCDM with state parallelization")
    if (st%d%nik > 1) call messages_not_implemented("SCDM with k-point sampling")
    if (der%mesh%sb%periodic_dim > 0 .and. der%mesh%sb%periodic_dim /= 3) &
         call messages_not_implemented("SCDM with mixed-periodicity")  

    ! determine whether we can apply scdm exchange operator to scdm states
    ! NOTE: this should be always the case, but for now only in td
    call parse_variable('CalculationMode', 0, inp_calc_mode)
    if(inp_calc_mode == 3) then
      scdm%psi_scdm = .true.
    else
      scdm%psi_scdm = .false.
    end if
    
    
#ifdef HAVE_MPI
    call MPI_Comm_Rank( der%mesh%mpi_grp%comm, rank, mpi_err)
#endif
    scdm%root = (rank ==0)

    ! inherit some indices from st
    scdm%st%d%dim = st%d%dim
    scdm%st%nst   = st%nst
    scdm%nst   = st%nst
    scdm%st%d%nik = st%d%nik
    scdm%st%d     = st%d
    scdm%cmplxscl = st%cmplxscl
    scdm%st%st_start = st%st_start
    scdm%st%st_end = st%st_end
    scdm%st%lnst = st%lnst

    !%Variable SCDM_reorthonormalize
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% If set to yes, and <tt>scdm_EXX = yes</tt>, the SCDM states are orthonormalized 
    !% on the domain defined by <tt>SCDMCutoffRadius<tt> (as opposed to the full simulation cell).
    !%End 
    call parse_variable('SCDM_reorthonormalize', .false., scdm%re_ortho_normalize)
    if (scdm%re_ortho_normalize) scdm%st%d%orth_method = ORTH_CHOLESKY_SERIAL

    !%Variable SCDM_verbose
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% Output detailed information on SCDM procedure.
    !%End
    call parse_variable('SCDM_verbose', .false., scdm%verbose)

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
    call parse_variable('SCDMCutoffRadius', 3._8, scdm%rcut, units_inp%length)
    if (scdm%root.and.scdm%verbose) call messages_print_var_value(stdout, 'SCDM cutoff', scdm%rcut)
    ! box_size is half the size of the  box
    scdm%box_size = 0
    do ii = 1, 3
      scdm%box_size = max(scdm%box_size,ceiling(scdm%rcut/der%mesh%spacing(ii)))
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
    SAFE_ALLOCATE(scdm%box(1:scdm%box_size*2+1,1:scdm%box_size*2+1,1:scdm%box_size*2+1,scdm%st%nst))

    ! the localzied states defined in the box are distributed over state index
    SAFE_ALLOCATE(istart(1:der%mesh%mpi_grp%size))
    SAFE_ALLOCATE(iend(1:der%mesh%mpi_grp%size))
    SAFE_ALLOCATE(ilsize(1:der%mesh%mpi_grp%size))

    call multicomm_divide_range(st%nst, der%mesh%mpi_grp%size, istart, iend, lsize=ilsize)
    scdm%st_start = istart(der%mesh%vp%rank+1)
    scdm%st_end = iend(der%mesh%vp%rank+1)
    scdm%lnst = ilsize(der%mesh%vp%rank+1)

    ! allocate local chunk of states
    ! localized SCDM states in full box (not really needed, but convenient for distribution)
    ! root process holds all states NOTE: this is not great... clearly ... but will go away when SCDM procedure is parallel
    if (.not.states_are_real(st)) then
      if (scdm%root) then
        SAFE_ALLOCATE(scdm%st%zdontusepsi(1:der%mesh%np_global,1:scdm%st%d%dim,1:scdm%st%nst,1:scdm%st%d%nik))
      else
        SAFE_ALLOCATE(scdm%st%zdontusepsi(1:der%mesh%np_global,1:scdm%st%d%dim,1:scdm%lnst,1:scdm%st%d%nik))
      end if
      ! localized SCDM states 
      SAFE_ALLOCATE(scdm%zpsi(1:scdm%full_box,1:scdm%lnst))
    else ! real
      if (scdm%root) then
        SAFE_ALLOCATE(scdm%st%ddontusepsi(1:der%mesh%np_global, 1:scdm%st%d%dim, 1:scdm%st%nst, 1:scdm%st%d%nik))
      else
        SAFE_ALLOCATE(scdm%st%ddontusepsi(1:der%mesh%np_global, 1:scdm%st%d%dim, 1:scdm%lnst, 1:scdm%st%d%nik))
      end if
      ! localized SCDM states
      SAFE_ALLOCATE(scdm%dpsi(1:scdm%full_box,1:scdm%lnst))
    end if
    
    SAFE_ALLOCATE(scdm%periodic(1:scdm%lnst))
    
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
!?
    scdm%boxmesh%idx%enlarge(:) = 0
    SAFE_ALLOCATE(scdm%boxmesh%idx%lxyz(1:scdm%boxmesh%np,1:scdm%boxmesh%idx%dim))
    ! need to copy indices because otherwise line gets too long (precompiler?)
    ii = -(scdm%box_size)
    jj = (scdm%box_size)
    SAFE_ALLOCATE(scdm%boxmesh%idx%lxyz_inv(ii:jj,ii:jj,ii:jj))
    !
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
    
    call mesh_cube_map_init(scdm%boxmesh%cube_map, scdm%boxmesh%idx, scdm%boxmesh%np_global)

    ! create a cube object for the small box, with double size for coulomb truncation
    box(1:3) = scdm%boxmesh%idx%ll(1:3)*2
    call cube_init(scdm%boxcube, box, scdm%boxmesh%sb,fft_type=FFT_REAL, fft_library=FFTLIB_FFTW)
    
    ! Joseba recommends including this
    !if (der%mesh%parallel_in_domains .and. this%cube%parallel_in_domains) then
    !    call mesh_cube_parallel_map_init(this%mesh_cube_map, der%mesh, this%cube)
    !end if
    
    ! set up poisson solver used for the exchange operator with scdm states
    ! this replictaes poisson_kernel_init()
    scdm%poisson%poisson_soft_coulomb_param = M_ZERO
    if (der%mesh%sb%periodic_dim.eq.3) then
      call poisson_fft_init(scdm%poisson_fft, scdm%boxmesh, scdm%boxcube, &
                            kernel=POISSON_FFT_KERNEL_HOCKNEY,fullcube=fullcube)
    else !non periodic case
      call poisson_fft_init(scdm%poisson_fft, scdm%boxmesh, scdm%boxcube, kernel=POISSON_FFT_KERNEL_SPH)
    end if

    ! create poisson object
    SAFE_ALLOCATE(scdm%poisson%der)
    SAFE_ALLOCATE(scdm%poisson%der%mesh)
    scdm%poisson%der%mesh = scdm%boxmesh
    scdm%poisson%method = POISSON_FFT
    scdm%poisson%kernel = POISSON_FFT_KERNEL_SPH
    scdm%poisson%cube = scdm%boxcube
    scdm%poisson%fft_solver = scdm%poisson_fft

    ! set flag to do this only once
    scdm_is_init = .true.

    call messages_write('done SCDM init')

    POP_SUB(scdm_init)
  end subroutine scdm_init

  !> wrapper routine to rotate  KS states into their SCDM representation
  subroutine scdm_rotate_states(st,mesh,scdm)
    type(states_t), intent(inout)  :: st
    type(mesh_t), intent(in)       :: mesh
    type(scdm_t), intent(inout)    :: scdm

    PUSH_SUB(scdm_rotate_states)

    if (.not.states_are_real(st)) then
      call zscdm_rotate_states(st,mesh,scdm)
    else
      call dscdm_rotate_states(st,mesh,scdm)
    end if

    POP_SUB(scdm_rotate_states)

  end subroutine scdm_rotate_states

  
  
  !> wrapper routine for real rank-revealing QR decompisition
  !! of the n*np matrix kst, returning the pivot vector jpvt
  subroutine dRRQR(nn, np, kst, jpvt)
    integer, intent(in)  :: nn
    integer, intent(in)  :: np
    FLOAT, intent(inout) :: kst(:,:)
    integer, intent(out) :: jpvt(np)

    integer            :: lwork, info
    FLOAT              :: tau(nn)
    FLOAT, allocatable :: work(:)

    PUSH_SUB(dRRQR)
    ! dummy call to obtain dimension of work
    SAFE_ALLOCATE(work(1:1))
    call DGEQP3(nn, np, kst, nn, jpvt, tau, work, -1, info )
    if (info /= 0) then
      write(message(1),'(A28,I2)') 'Illegal argument in DGEQP3: ', info
      call messages_fatal(1)
    end if
    ! Note: scalapack routine is called P?GEQPF()

    lwork = work(1)
    SAFE_DEALLOCATE_A(work)
    SAFE_ALLOCATE(work(1:lwork))

    jpvt(:) = 0
    tau(:) = 0.
    ! actual call
    call DGEQP3(nn, np, kst, nn, jpvt, tau, work, lwork, info)
    if (info /= 0) then
      write(message(1),'(A28,I2)') 'Illegal argument in DGEQP3: ', info
      call messages_fatal(1)
    end if

    POP_SUB(dRRQR)
  end subroutine dRRQR

  !> wrapper routine for complex rank-revealing QR decomposition
  !! of the n*np matrix kst, returning the pivot vector jpvt
  subroutine zRRQR(nn, np, kst, jpvt)
    integer, intent(in)  :: nn
    integer, intent(in)  :: np
    CMPLX, intent(inout) :: kst(:,:)
    integer, intent(out) :: jpvt(np)

    integer            :: lwork,info
    CMPLX              :: tau(nn)
    CMPLX, allocatable :: work(:)
    FLOAT              :: rwork(2*np)

    PUSH_SUB(zRRQR)
    ! dummy call to obtain dimension of work
    SAFE_ALLOCATE(work(1:1))
    call ZGEQP3(nn, np, kst, nn, jpvt, tau, work, -1, rwork, info)
    if (info /= 0) then
      write(message(1),'(A28,I2)') 'Illegal argument in ZGEQP3: ', info
      call messages_fatal(1)
    end if
    ! Note: scalapack routine is called P?GEQPF()

    lwork = work(1)
    SAFE_DEALLOCATE_A(work)
    SAFE_ALLOCATE(work(1:lwork))

    jpvt(:) = 0
    tau(:) = 0.
    ! actual call
    call ZGEQP3(nn, np, kst, nn, jpvt, tau, work, lwork, rwork, info)
    if (info /= 0)then
      write(message(1),'(A28,I2)') 'Illegal argument in ZGEQP3: ', info
      call messages_fatal(1)
    end if

    POP_SUB(zRRQR)
  end subroutine zRRQR

  !> check if there are points outside index range of idx by
  !! checking the corners only. This is intended for rectangular cells
  !! should be generalized to arbitrary shapes (but then need to check the faces)
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

end module scdm_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
