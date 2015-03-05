!> @file
!! Wrapper for the MPI call (this file is preprocessed.)
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!!JOSEBA#if defined HAVE_CONFIG_H
!!JOSEBA#include <config.inc>
!!JOSEBA#endif
#include "config.h" 
!! JOSEBA


!> Module defining the routines which wrap the MPI calls
module wrapper_MPI
  ! TO BE REMOVED with f_malloc
  !use memory_profiling!, only: ndebug
  ! TO BE REMOVED with f_malloc
  use time_profiling, only: TIMING_UNINITIALIZED

  implicit none

  ! MPI handling
#ifdef HAVE_MPI2
  logical, parameter :: have_mpi2 = .true.  !< Flag to use in the code to switch between MPI1 and MPI2
#else
  integer :: MPI_IN_PLACE = 0               !< Fake MPI_IN_PLACE variable to allow compilation in sumrho.
  logical, parameter :: have_mpi2 = .false. !< Flag to use in the code to switch between MPI1 and MPI2
#endif

  !> MPI definitions and datatypes for density and wavefunctions
  #include "mpif.h"

  logical :: mpi_thread_funneled_is_supported=.false. !< Control the OMP_NESTED based overlap, checked by bigdft_mpi_init below

  !timing categories for MPI wrapper
  integer, parameter :: smallsize=5 !< limit for a communication with small size
  character(len=*), parameter, public :: tgrp_mpi_name='Communications'
  !timing categories
  integer, public, save :: TCAT_ALLRED_SMALL=TIMING_UNINITIALIZED
  integer, public, save :: TCAT_ALLRED_LARGE=TIMING_UNINITIALIZED
  integer, public, save :: TCAT_ALLGATHERV  =TIMING_UNINITIALIZED
  integer, public, save :: TCAT_GATHER      =TIMING_UNINITIALIZED
  
  !error codes
  integer, public, save :: ERR_MPI_WRAPPERS

  !> Interface for MPITYPE routine
  interface mpitype
     module procedure mpitype_i,mpitype_d,mpitype_r,mpitype_l,mpitype_c,mpitype_li
     module procedure mpitype_i1,mpitype_i2
     module procedure mpitype_d1,mpitype_d2
     module procedure mpitype_c1
  end interface mpitype

  interface mpimaxdiff
     module procedure mpimaxdiff_i0,mpimaxdiff_d0
     module procedure mpimaxdiff_i1,mpimaxdiff_i2
     module procedure mpimaxdiff_d1,mpimaxdiff_d2
  end interface mpimaxdiff

  !> Interface for MPI_ALLREDUCE routine, to be updated little by little
  interface mpiallred
     module procedure mpiallred_int,mpiallred_real, &
          & mpiallred_double,&!,mpiallred_double_1,mpiallred_double_2,&
          & mpiallred_log
     module procedure mpiallred_d1,mpiallred_d2
  end interface mpiallred

  interface mpigather
     module procedure mpigather_d1d1,mpigather_d2d1,mpigather_d1d2,mpigather_d2
     module procedure mpigather_i0i2,mpigather_d0d2,mpigather_i1i2,mpigather_i2
  end interface mpigather

  interface mpibcast
     module procedure mpibcast_i0,mpibcast_li0,mpibcast_d0
     module procedure mpibcast_c1,mpibcast_d1,mpibcast_d2,mpibcast_i1
  end interface mpibcast

#ifdef HAVE_MPI2
  interface mpi_get_to_allgatherv
     module procedure mpi_get_to_allgatherv_double
  end interface mpi_get_to_allgatherv

  interface mpiget
    module procedure mpiget_d0
  end interface mpiget

  interface mpiwindow
    module procedure mpiwindow_d0
  end interface mpiwindow
#endif

  interface mpitypesize
    module procedure mpitypesize_d0, mpitypesize_d1
  end interface mpitypesize

  !> Interface for MPI_ALLGATHERV routine
  interface mpiallgatherv
     module procedure mpiallgatherv_double
  end interface mpiallgatherv
  
  interface mpiiallred
      module procedure mpiiallred_double
  end interface mpiiallred

  interface mpiialltoallv
      module procedure mpiialltoallv_double
  end interface mpiialltoallv

  !> Global MPI communicator which contains all information related to the MPI process
  type, public :: mpi_environment
     integer :: mpi_comm !< MPI communicator
     integer :: iproc    !< Process Id
                         !! @ingroup RESERVED
     integer :: nproc    !< Number of MPI processes (in the given communicator)
                         !! @ingroup RESERVED
     integer :: igroup   !< MPI Group Id
     integer :: ngroup   !< Number of MPI groups
  end type mpi_environment

  public :: mpi_environment_null
  public :: mpi_environment_free
  public :: mpi_environment_set
  public :: mpi_environment_set1 !to be removed

  !>fake type to enhance documentation
  type, private :: doc
     !>number of entries in buffer (integer). Useful for buffer passed by reference 
     integer :: count
     !> rank of mpitask executing the operation (default value is root=0)
     integer :: root
     !> communicator of the communication
     integer :: comm
  end type doc

contains

  pure function mpi_environment_null() result(mpi)
    implicit none
    type(mpi_environment) :: mpi
    mpi%mpi_comm=MPI_COMM_NULL !better to put an invalid comm?
    mpi%igroup=-1
    mpi%ngroup=-1
    mpi%iproc=-1
    mpi%nproc=-1
  end function mpi_environment_null

  subroutine mpi_environment_free(mpi_env)
    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    implicit none
    type(mpi_environment), intent(inout) :: mpi_env
    !local variables
    integer :: ierr

    if (mpi_env%mpi_comm /= MPI_COMM_WORLD .and. &
         mpi_env%mpi_comm /= MPI_COMM_NULL) then
       call MPI_COMM_FREE(mpi_env%mpi_comm,ierr)
       if (ierr /=0) then
          call f_err_throw('Problem in MPI_COMM_FREE, ierr:'//&
               yaml_toa(ierr),err_name='BIGDFT_MPI_ERROR')
          return
       end if
    end if
    mpi_env=mpi_environment_null()
  end subroutine mpi_environment_free


  !> Set the MPI environment (i.e. taskgroup or MPI communicator)
  subroutine mpi_environment_set(mpi_env,iproc,nproc,mpi_comm,groupsize)
    use dynamic_memory
    use yaml_output
    implicit none
    integer, intent(in) :: iproc     !<  Proc id
    integer, intent(in) :: nproc     !<  Total number of MPI processes
    integer, intent(in) :: mpi_comm  !<  Global MPI_communicator
    integer, intent(in) :: groupsize !<  Number of MPI processes by (task)group
                                     !!  if 0 one taskgroup (MPI_COMM_WORLD)   
    type(mpi_environment), intent(out) :: mpi_env  !< MPI environment (out)
    !local variables
    integer :: j
    integer, dimension(:), allocatable :: group_list

    call f_routine(id='mpi_environment_set')
    mpi_env=mpi_environment_null()
!!$    mpi_env%igroup=0
!!$    mpi_env%ngroup=1
!!$    mpi_env%iproc=iproc
!!$    mpi_env%nproc=nproc
    mpi_env%mpi_comm=mpi_comm
    mpi_env%igroup=iproc/groupsize
    mpi_env%ngroup=nproc/groupsize
    mpi_env%iproc=mod(iproc,groupsize)
    mpi_env%nproc=groupsize
    if (groupsize /= nproc) then
       !define the strategy for the taskgroups
       group_list=f_malloc(groupsize,id='group_list')
       !iproc in the same group are close to each other
       do j=0,groupsize-1
          group_list(j+1)=mpi_env%igroup*groupsize+j
       enddo
       call create_group_comm(mpi_comm,groupsize,group_list,mpi_env%mpi_comm)
       if (iproc == 0) then
          call yaml_map('Total No. of Taskgroups created',nproc/mpi_env%nproc)
       end if
       call f_free(group_list)
    end if

    call f_release_routine()
  end subroutine mpi_environment_set


!!! PSolver n1-n2 plane mpi partitioning !!! 
  !> This is exactly like mpi_environment_set but it always creates groups
  !! the routine above should be modified accordingly
!!$  subroutine mpi_environment_set2(mpi_env,iproc,nproc,mpi_comm,groupsize)
!!$    use yaml_output
!!$    implicit none
!!$    integer, intent(in) :: iproc,nproc,mpi_comm,groupsize
!!$    type(mpi_environment), intent(out) :: mpi_env
!!$    !local variables
!!$    integer :: j
!!$    integer, dimension(:), allocatable :: group_list
!!$
!!$    call f_routine(id='mpi_environment_set2')
!!$    mpi_env=mpi_environment_null()
!!$
!!$    mpi_env%mpi_comm=mpi_comm
!!$
!!$    mpi_env%igroup=iproc/groupsize
!!$    mpi_env%ngroup=nproc/groupsize
!!$    mpi_env%iproc=mod(iproc,groupsize)
!!$    mpi_env%nproc=groupsize
!!$
!!$    !define the strategy for the taskgroups
!!$    group_list=f_malloc(groupsize,id='group_list')
!!$    !iproc in the same group are close to each other
!!$    do j=0,groupsize-1
!!$       group_list(j+1)=mpi_env%igroup*groupsize+j
!!$    enddo
!!$
!!$    call create_group_comm(mpi_comm,nproc,mpi_env%igroup,mpi_env%nproc,group_list,mpi_env%mpi_comm)
!!$!    if (iproc == 0) then
!!$!       call yaml_map('Total No. of Taskgroups created',nproc/mpi_env%nproc)
!!$!    end if
!!$    call f_free(group_list)
!!$    call f_release_routine()
!!$  end subroutine mpi_environment_set2


  !> This is a different procedure to assign the iproc according to the groups.
  subroutine mpi_environment_set1(mpi_env,iproc,mpi_comm,groupsize,ngroup)
    use yaml_output
    use dynamic_memory
    implicit none
    integer, intent(in) :: iproc,mpi_comm,groupsize,ngroup
    type(mpi_environment), intent(out) :: mpi_env
    !local variables
    integer :: j
    integer, dimension(:), allocatable :: group_list

    call f_routine(id='mpi_environment_set1')

    mpi_env=mpi_environment_null()

    mpi_env%igroup=-1

    mpi_env%ngroup=ngroup
    if (iproc < groupsize*ngroup) mpi_env%igroup=mod(iproc,ngroup)
    mpi_env%iproc=iproc/ngroup
    mpi_env%nproc=groupsize
    mpi_env%mpi_comm=mpi_comm

    !define the strategy for the taskgroups
    group_list=f_malloc(groupsize,id='group_list')
    !round-robin strategy
    if (mpi_env%igroup >0) then
       do j=0,groupsize-1
          group_list(j+1)=mpi_env%igroup+j*mpi_env%ngroup
       enddo
    else
       !these processes have MPI_COMM_NULL
       group_list=-1
       mpi_env%mpi_comm=MPI_COMM_NULL
    end if

    !call create_group_comm1(mpi_comm,nproc,mpi_env%igroup,ngroup,mpi_env%nproc,mpi_env%mpi_comm)
    call create_group_comm(mpi_comm,mpi_env%nproc,group_list,mpi_env%mpi_comm)
!    if (iproc == 0) then
!       call yaml_map('Total No. of Taskgroups created',ngroup)
!    end if
    call f_free(group_list)
    call f_release_routine()
  end subroutine mpi_environment_set1


  !> Create communicators associated to the groups of size group_size
  subroutine create_group_comm(base_comm,group_size,group_list,group_comm)
    implicit none
    integer, intent(in) :: base_comm,group_size
    integer, dimension(group_size), intent(in) :: group_list !< list of id of the group identified by group_id in units of base_comm
    integer, intent(out) :: group_comm
    !local variables
    integer :: grp,ierr,base_grp

    !take the base group
    call MPI_COMM_GROUP(base_comm,base_grp,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'group creation')
       return
    end if
    !create the groups with the list
    call MPI_GROUP_INCL(base_grp,group_size,group_list,grp,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'group inclusion')
       return
    end if
    !free base group
    call MPI_GROUP_FREE(base_grp,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'base_group free')
       return
    end if
    !create the communicator (the communicator can be also null)
    call MPI_COMM_CREATE(base_comm,grp,group_comm,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'communicator creator')
       return
    end if
    !free temporary group
    call MPI_GROUP_FREE(grp,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'new_group free')
       return
    end if

    contains
      
      subroutine check_ierr(ierr,message)
        use yaml_output, only: yaml_toa
        use dictionaries, only: f_err_throw
        implicit none
        integer, intent(in) :: ierr
        character(len=*), intent(in) :: message
        if (ierr /= 0) then
           call f_err_throw('Problem in '//trim(message)//&
                ', ierr:'//yaml_toa(ierr),err_name='BIGDFT_MPI_ERROR')
        end if
      end subroutine check_ierr

  end subroutine create_group_comm


  !!! PSolver n1-n2 plane mpi partitioning !!! 
  !> This routine is like create_group_comm with a different group_list
  subroutine create_group_comm1(base_comm,group_id,ngroup,group_size,group_comm)
    use dynamic_memory
    use yaml_output
    implicit none
    integer, intent(in) :: base_comm,group_size,group_id,ngroup
    integer, intent(out) :: group_comm
    !local variables
    character(len=*), parameter :: subname='create_group_comm'
    integer :: grp,ierr,i,j,base_grp,temp_comm!,i_stat,i_all
    integer, dimension(:), allocatable :: group_list

  ! allocate(group_list(group_size+ndebug),stat=i_stat)
    group_list = f_malloc(group_size,id='group_list')

    !take the base group
    call MPI_COMM_GROUP(base_comm,base_grp,ierr)
    if (ierr /=0) then
       call yaml_warning('Problem in group creation, ierr:'//yaml_toa(ierr))
       call MPI_ABORT(base_comm,1,ierr)
    end if
    do i=0,ngroup-1
       !define the new groups and thread_id
       do j=0,group_size-1
          group_list(j+1)=i+j*ngroup
       enddo
       call MPI_GROUP_INCL(base_grp,group_size,group_list,grp,ierr)
       if (ierr /=0) then
          call yaml_warning('Problem in group inclusion, ierr:'//yaml_toa(ierr))
          call MPI_ABORT(base_comm,1,ierr)
       end if
       call MPI_COMM_CREATE(base_comm,grp,temp_comm,ierr)
       if (ierr /=0) then
          call yaml_warning('Problem in communicator creator, ierr:'//yaml_toa(ierr))
          call MPI_ABORT(base_comm,1,ierr)
       end if
       !print *,'i,group_id,temp_comm',i,group_id,temp_comm
       if (i.eq. group_id) group_comm=temp_comm
    enddo

  !i_all=-product(shape(group_list ))*kind(group_list )
  ! deallocate(group_list,stat=i_stat)
    call f_free(group_list)
  end subroutine create_group_comm1


  !> Create a communicator between proc of same rank between the taskgroups.
  subroutine create_rank_comm(group_comm, rank_comm)
    use dynamic_memory
    use yaml_output
    implicit none
    integer, intent(in) :: group_comm
    integer, intent(out) :: rank_comm
    !local variables
    character(len=*), parameter :: subname='create_group_master'
    integer :: iproc_group, nproc, nproc_group, ngroups
    integer :: ierr, i, j
    integer, dimension(:), allocatable :: lrank, ids

    call MPI_COMM_RANK(group_comm, iproc_group, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
    call MPI_COMM_SIZE(group_comm, nproc_group, ierr)
    ngroups = nproc / nproc_group

    ! Put in lrank the group rank of each process, indexed by global iproc.
!   allocate(lrank(nproc+ndebug), stat = i_stat)
    lrank = f_malloc(nproc,id='lrank')
    call mpi_allgather(iproc_group, 1, MPI_INTEGER, lrank, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Put in ids, the global iproc of each process that share the same group iproc.
!   allocate(ids(ngroups+ndebug), stat = i_stat)
    ids = f_malloc(ngroups,id='ids')
    j = 1
    do i = 1, nproc
       if (lrank(i) == iproc_group) then
          ids(j) = i - 1
          j = j + 1
       end if
    end do
!  i_all=-product(shape(lrank ))*kind(lrank )
!   deallocate(lrank,stat=i_stat)
    call f_free(lrank)

!!$    call mpi_comm_rank(MPI_COMM_WORLD, iproc_group, ierr)
!!$    write(*,*) iproc_group, "->", ids
    
    ! Create a new comminucator for the list of ids.
    call create_group_comm(MPI_COMM_WORLD, ngroups, ids, rank_comm)
!  i_all=-product(shape(ids ))*kind(ids )
!   deallocate(ids,stat=i_stat)
    call f_free(ids)
  END SUBROUTINE create_rank_comm


  subroutine wmpi_init_thread(ierr)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(out) :: ierr
#ifdef HAVE_MPI_INIT_THREAD
    integer :: provided
    call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
    if (ierr /= MPI_SUCCESS) then
       write(*,*)'BigDFT_mpi_INIT: Error in MPI_INIT_THREAD',ierr
    else if (provided < MPI_THREAD_FUNNELED) then
       !write(*,*)'WARNING: MPI_THREAD_FUNNELED not supported!',provided,ierr
       !call MPI_INIT(ierr)
    else
       mpi_thread_funneled_is_supported=.true.
    endif
#else
    call MPI_INIT(ierr)      
    if (ierr /= MPI_SUCCESS) then
       call f_err_throw('An error in calling to MPI_INIT (THREAD) occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
#endif
  end subroutine wmpi_init_thread


  !> Finalization of the mpi
  subroutine mpifinalize()
    use dictionaries, only: f_err_throw
    implicit none
    !local variables
    integer :: ierr

    call MPI_FINALIZE(ierr)
    if (ierr /= MPI_SUCCESS) then
       call f_err_throw('An error in calling to MPI_INIT_THREAD occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpifinalize


  !> Initialize timings and also mpi errors
  subroutine mpi_initialize_timing_categories()
    use time_profiling, only: f_timing_category_group,f_timing_category
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none

    call f_timing_category_group(tgrp_mpi_name,&
         'Operations between MPI tasks')

    call f_timing_category('Allreduce, Small Size',tgrp_mpi_name,&
         'Allreduce operations for less than'//&
         trim(yaml_toa(smallsize))//' elements',&
         TCAT_ALLRED_SMALL)
    call f_timing_category('Allreduce, Large Size',tgrp_mpi_name,&
         'Allreduce operations for more than'//&
         trim(yaml_toa(smallsize))//' elements',&
         TCAT_ALLRED_LARGE)
    call f_timing_category('Allgatherv',tgrp_mpi_name,&
         'Variable allgather operations',&
         TCAT_ALLGATHERV)
    call f_timing_category('Gather',tgrp_mpi_name,&
         'Gather operations, in general moderate size arrays',&
         TCAT_GATHER)

    call f_err_define(err_name='ERR_MPI_WRAPPERS',err_msg='Error of MPI library',&
         err_id=ERR_MPI_WRAPPERS,&
         err_action='Some MPI library returned an error code, inspect runtime behaviour')

  end subroutine mpi_initialize_timing_categories

  pure function mpitype_i(data) result(mt)
    implicit none
    integer, intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER
  end function mpitype_i
  pure function mpitype_i1(data) result(mt)
    implicit none
    integer, dimension(:), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER
  end function mpitype_i1
  pure function mpitype_i2(data) result(mt)
    implicit none
    integer, dimension(:,:), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER
  end function mpitype_i2

  pure function mpitype_li(data) result(mt)
    implicit none
    integer(kind=8), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER8
  end function mpitype_li


  pure function mpitype_r(data) result(mt)
    implicit none
    real, intent(in) :: data
    integer :: mt
    mt=MPI_REAL
  end function mpitype_r
  pure function mpitype_d(data) result(mt)
    implicit none
    double precision, intent(in) :: data
    integer :: mt
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d
  pure function mpitype_d1(data) result(mt)
    implicit none
    double precision, dimension(:), intent(in) :: data
    integer :: mt
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d1
  pure function mpitype_d2(data) result(mt)
    implicit none
    double precision, dimension(:,:), intent(in) :: data
    integer :: mt
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d2
  pure function mpitype_l(data) result(mt)
    implicit none
    logical, intent(in) :: data
    integer :: mt
    mt=MPI_LOGICAL
  end function mpitype_l
  pure function mpitype_c(data) result(mt)
    implicit none
    character, intent(in) :: data
    integer :: mt
    mt=MPI_CHARACTER
  end function mpitype_c
  pure function mpitype_c1(data) result(mt)
    implicit none
    character, dimension(:), intent(in) :: data
    integer :: mt
    mt=MPI_CHARACTER
  end function mpitype_c1

  !> Function giving the mpi rank id for a given communicator
  function mpirank(comm)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(in) :: comm
    integer :: mpirank
    !local variables
    integer :: iproc,ierr

    call MPI_COMM_RANK(comm, iproc, ierr)
    if (ierr /=0) then
       iproc=-1
       mpirank=iproc
       call f_err_throw('An error in calling to MPI_COMM_RANK occurred',&
            err_id=ERR_MPI_WRAPPERS)
    end if
    mpirank=iproc

  end function mpirank

  !> Returns the number of mpi_tasks associated to a given communicator
  function mpisize(comm)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(in) :: comm
    integer :: mpisize
    !local variables
    integer :: nproc,ierr

    !verify the size of the receive buffer
    call MPI_COMM_SIZE(comm,nproc,ierr)
    if (ierr /=0) then
       nproc=0
       mpisize=nproc
       call f_err_throw('An error in calling to MPI_COMM_SIZE occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
    mpisize=nproc

  end function mpisize

  !> Performs the barrier of a given communicator, if present
  subroutine mpibarrier(comm)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(in), optional :: comm !< the communicator
    !local variables
    integer :: mpi_comm,ierr
    
    if (present(comm)) then
       mpi_comm=comm
    else
       mpi_comm=MPI_COMM_WORLD
    end if
    !call the barrier
    call MPI_BARRIER(mpi_comm,ierr)
    if (ierr /=0) then
       call f_err_throw('An error in calling to MPI_BARRIER occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpibarrier

  !> Gather the results of a given array into the root proc
  subroutine mpigather_d1d1(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    double precision, dimension(:), intent(in) :: sendbuf
    double precision, dimension(:), intent(inout) :: recvbuf
    include 'gather-inc.f90'   
  end subroutine mpigather_d1d1

  subroutine mpigather_d1d2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    double precision, dimension(:), intent(in) :: sendbuf
    double precision, dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'   
  end subroutine mpigather_d1d2

  subroutine mpigather_i1i2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    integer, dimension(:), intent(in) :: sendbuf
    integer, dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'   
  end subroutine mpigather_i1i2

  subroutine mpigather_i2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    integer, dimension(:,:), intent(in) :: sendbuf
    integer, dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'   
  end subroutine mpigather_i2


  subroutine mpigather_d2d1(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    double precision, dimension(:,:), intent(in) :: sendbuf
    double precision, dimension(:), intent(inout) :: recvbuf
    include 'gather-inc.f90'   
  end subroutine mpigather_d2d1

  subroutine mpigather_d2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    double precision, dimension(:,:), intent(in) :: sendbuf
    double precision, dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'   
  end subroutine mpigather_d2

  !> Gather the results of a given array into the root proc, version 
  !! working with adresses
  subroutine mpigather_i0i2(sendbuf,sendcount,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    integer, intent(inout) :: sendbuf
    integer, intent(in) :: sendcount
    integer, dimension(:,:), intent(inout) :: recvbuf
    !---like gather-inc
    integer, intent(in), optional :: root !< 0 if absent
    integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
    !local variables
    integer :: iroot,mpi_comm,ntot,ntotrecv,ntasks,ierr

    ntot=sendcount
    ntotrecv=size(recvbuf)

    include 'gather-inner-inc.f90'
    !-end gather-inc
  end subroutine mpigather_i0i2

  subroutine mpigather_d0d2(sendbuf,sendcount,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    double precision, intent(inout) :: sendbuf
    integer, intent(in) :: sendcount
    double precision, dimension(:,:), intent(inout) :: recvbuf
    !---like gather-inc
    integer, intent(in), optional :: root !< 0 if absent
    integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
    !local variables
    integer :: iroot,mpi_comm,ntot,ntotrecv,ntasks,ierr

    ntot=sendcount
    ntotrecv=size(recvbuf)

    include 'gather-inner-inc.f90'
    !-end gather-inc
  end subroutine mpigather_d0d2
  


  !> Interface for MPI_ALLGATHERV operations
  subroutine mpiallgatherv_double(buffer,counts,displs,me,mpi_comm,ierr)
    use dynamic_memory
    implicit none
    integer, dimension(0:), intent(in) :: counts
    integer, dimension(:), intent(in) :: displs
    integer, intent(in) :: mpi_comm, me
    real(kind=8), intent(inout) :: buffer
    integer, intent(out) :: ierr
#ifdef HAVE_MPI2
    call f_timer_interrupt(TCAT_ALLGATHERV)
    !case with MPI_IN_PLACE
    call MPI_ALLGATHERV(MPI_IN_PLACE,counts(me),mpitype(buffer),&
         buffer,counts,displs,mpitype(buffer),mpi_comm,ierr)
    call f_timer_resume()
#else
    !local variables
    real(kind=8), dimension(:), allocatable :: copybuf

    !Here we have a performance penalty by copying all buffer, instead of
    !just the send part, but I don`t see how to get buffer(displs(me))
    copybuf = f_malloc(sum(counts),id='copybuf')

    call dcopy(sum(counts),buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call f_timer_interrupt(TCAT_ALLGATHERV)
    call MPI_ALLGATHERV(copybuf(1+displs(me+1)),counts(me),mpitype(buffer),&
         buffer,counts,displs,mpitype(buffer),mpi_comm,ierr)
    call f_timer_resume()
    call f_free(copybuf)
#endif

    if (ierr /=0) stop 'MPIALLGATHERV_DBL'
  end subroutine mpiallgatherv_double

  !> Interface for MPI_ALLREDUCE operations
  subroutine mpiallred_int(sendbuf,count,op,comm,recvbuf)
    use dictionaries, only: f_err_throw,f_err_define
    use dynamic_memory
    implicit none
    integer, intent(inout) :: sendbuf
    integer, intent(inout), optional :: recvbuf
    integer, dimension(:), allocatable :: copybuf
    include 'allreduce-inc.f90'
  end subroutine mpiallred_int

  !> Interface for MPI_ALLREDUCE operations
  subroutine mpiallred_real(sendbuf,count,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    real, intent(inout) :: sendbuf
    real, intent(inout), optional :: recvbuf
    real, dimension(:), allocatable :: copybuf
    include 'allreduce-inc.f90'
  end subroutine mpiallred_real

  subroutine mpiallred_double(sendbuf,count,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    double precision, intent(inout) :: sendbuf
    double precision, intent(inout), optional :: recvbuf
    double precision, dimension(:), allocatable :: copybuf
    include 'allreduce-inc.f90'
  end subroutine mpiallred_double

  subroutine mpiallred_log(sendbuf,count,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    logical, intent(inout) :: sendbuf
    logical, intent(inout), optional :: recvbuf
    logical, dimension(:), allocatable :: copybuf
    include 'allreduce-inc.f90'
  end subroutine mpiallred_log

  subroutine mpiallred_d1(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    double precision, dimension(:), intent(inout) :: sendbuf
    double precision, dimension(:), intent(inout), optional :: recvbuf
    double precision, dimension(:), allocatable :: copybuf  
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_d1

  subroutine mpiallred_d2(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    double precision, dimension(:,:), intent(inout) :: sendbuf
    double precision, dimension(:,:), intent(inout), optional :: recvbuf
    double precision, dimension(:,:), allocatable :: copybuf  
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_d2

  recursive subroutine mpibcast_i0(buffer,count,root,comm,check)
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    implicit none
    integer, intent(inout) ::  buffer 
    include 'bcast-decl-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_i0

  subroutine mpibcast_li0(buffer,count,root,comm,check)
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    implicit none
    integer(kind=8), intent(inout) ::  buffer      
    include 'bcast-decl-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_li0

  recursive subroutine mpibcast_d0(buffer,count,root,comm,check)
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    implicit none
    double precision, intent(inout) ::  buffer 
    include 'bcast-decl-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d0

  subroutine mpibcast_c1(buffer,root,comm,check)
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    implicit none
    character, dimension(:), intent(inout) ::  buffer      
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_c1

  subroutine mpibcast_i1(buffer,root,comm,check)
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    implicit none
    integer, dimension(:), intent(inout) ::  buffer      
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_i1

  subroutine mpibcast_d1(buffer,root,comm,check)
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    implicit none
    double precision, dimension(:), intent(inout) ::  buffer      
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d1

  subroutine mpibcast_d2(buffer,root,comm,check)
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    implicit none
    double precision, dimension(:,:), intent(inout) ::  buffer      
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d2


  !> Detect the maximum difference between arrays all over a given communicator
  function mpimaxdiff_i0(n,array,root,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    integer, intent(in) :: n !<number of elements to be controlled
    integer, intent(inout) :: array !< starting point of the array
    integer, dimension(:,:), allocatable :: array_glob
    integer :: maxdiff 
    include 'maxdiff-decl-inc.f90'

    ndims = n
    maxdiff=0

    include 'maxdiff-inc.f90'
  end function mpimaxdiff_i0

  function mpimaxdiff_d0(n,array,root,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    integer, intent(in) :: n !<number of elements to be controlled
    double precision, intent(inout) :: array !< starting point of the array
    double precision, dimension(:,:), allocatable :: array_glob
    double precision :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = n
    maxdiff=0.d0

    include 'maxdiff-inc.f90'
  end function mpimaxdiff_d0

  function mpimaxdiff_d1(array,root,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    double precision, dimension(:), intent(in) :: array 
    double precision, dimension(:,:), allocatable :: array_glob
    double precision :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0.d0
    
    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_d1

  function mpimaxdiff_i1(array,root,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    integer, dimension(:), intent(in) :: array 
    integer, dimension(:,:), allocatable :: array_glob
    integer :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0
    
    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_i1

  function mpimaxdiff_i2(array,root,comm,bcast) result(maxdiff)

    use dynamic_memory
    implicit none
    !> array to be checked
    integer, dimension(:,:), intent(in) :: array 
    integer, dimension(:,:), allocatable :: array_glob
    integer :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0
    
    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_i2


  function mpimaxdiff_d2(array,root,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    double precision, dimension(:,:), intent(in) :: array 
    double precision, dimension(:,:), allocatable :: array_glob
    double precision :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0.d0
    
   include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_d2

  !!function mpitypesize_d(foo) result(sizeof)
  !!  use dictionaries, only: f_err_throw,f_err_define
  !!  implicit none
  !!  double precision, intent(in) :: foo
  !!  integer :: sizeof, ierr

  !!  call mpi_type_size(mpi_double_precision, sizeof, ierr)
  !!  if (ierr/=0) then
  !!      call f_err_throw('Error in mpi_type_size',&
  !!           err_id=ERR_MPI_WRAPPERS)
  !!  end if
  !!end function mpitypesize_d

  function mpitypesize_d0(foo) result(sizeof)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    double precision, intent(in) :: foo
    integer :: sizeof, ierr
    
    call mpi_type_size(mpi_double_precision, sizeof, ierr)
    if (ierr/=0) then
        call f_err_throw('Error in mpi_type_size',&
             err_id=ERR_MPI_WRAPPERS)
    end if
  end function mpitypesize_d0

  function mpitypesize_d1(foo) result(sizeof)
      implicit none
      double precision, dimension(:), intent(in) :: foo
      integer :: sizeof
      sizeof=mpitypesize(1.d0)
  end function mpitypesize_d1

  function mpiinfo(key,val) result(info)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: val
    integer :: info, ierr
    
    call mpi_info_create(info, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_info_create',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
    call mpi_info_set(info, "no_locks", "true", ierr)
    if (ierr/=0) then
       !!call f_err_throw('Error in mpi_info_set, key='//trim(key)//&
       !!     ', value=',trim(val),err_id=ERR_MPI_WRAPPERS)
       call f_err_throw('Error in mpi_info_set, key='//trim(key)//&
            ', value='//trim(val),err_id=ERR_MPI_WRAPPERS)
    end if
    
  end function mpiinfo

  subroutine mpiinfofree(info)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer, intent(inout) :: info
    ! Local variables
    integer :: ierr
    call mpi_info_free(info, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_info_free',&
            err_id=ERR_MPI_WRAPPERS)
   end if
  end subroutine mpiinfofree

#ifdef HAVE_MPI2
  function mpiwindow_d0(size,base,comm) result(window)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer,intent(in) :: size
    double precision,intent(in) :: base
    integer,intent(in) :: comm
    !local variables
    integer :: sizeof,info,ierr
    integer :: window

    sizeof=mpitypesize(base)
    info=mpiinfo("no_locks", "true")

    call mpi_win_create(base, int(size,kind=mpi_address_kind)*int(sizeof,kind=mpi_address_kind), &
         sizeof, info,comm, window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_create',&
            err_id=ERR_MPI_WRAPPERS)
    end if

    call mpiinfofree(info)

    call mpi_win_fence(MPI_MODE_NOPRECEDE, window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_fence',&
            err_id=ERR_MPI_WRAPPERS)
    end if

    
  end function mpiwindow_d0

  subroutine mpi_fenceandfree(window)
    use dictionaries, only: f_err_throw,f_err_define
    ! Calling arguments
    integer,intent(inout) :: window !<window to be synchronized and freed

    ! Local variables
    integer :: ierr

    ! Synchronize the communication
    call mpi_win_fence(0, window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_fence',&
            err_id=ERR_MPI_WRAPPERS)  
    end if
    call mpi_win_free(window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_fence',&
            err_id=ERR_MPI_WRAPPERS)  
    end if
  end subroutine mpi_fenceandfree

  subroutine mpiget_d0(origin,count,target_rank,target_disp,window)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    double precision,intent(inout) :: origin !<fake intent(in)
    integer,intent(in) :: count, target_rank,window
    integer(kind=mpi_address_kind),intent(in) :: target_disp

    ! Local variables
    integer :: ierr

    call mpi_get(origin,count,mpitype(1.d0),target_rank, &
         target_disp,count,mpitype(origin), window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_get',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpiget_d0

  subroutine mpi_get_to_allgatherv_double(sendbuf,sendcount,recvbuf,recvcounts,displs,comm,check_,window_)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_output, only: yaml_toa
    implicit none
    !!double precision,dimension(:),intent(in) :: sendbuf
    !!double precision,dimension(:),intent(inout) :: recvbuf
    double precision,intent(in) :: sendbuf
    double precision,intent(inout) :: recvbuf
    integer,dimension(:),intent(in) :: recvcounts, displs
    integer,intent(in) :: comm, sendcount
    logical,intent(in),optional :: check_
    integer,intent(out),pointer,optional :: window_
    !local variables
    integer :: nproc,jproc,nrecvbuf,ierr
    external :: getall
    logical :: check
    integer,target:: window

    nproc=mpisize(comm)
    nrecvbuf=sum(recvcounts)

    if (present(check_)) then
        check = check_
    else
        check = .false.
    end if

    if (check) then
       !check coherence
       if (any([size(recvcounts),size(displs)] /= nproc)) then
          call f_err_throw("Error in get_to_gatherv, sizes not coherent with communicator"//&
               trim(yaml_toa([size(recvcounts),size(displs), nproc])),&
               err_id=ERR_MPI_WRAPPERS)
          return
       end if
    end if


    if (present(window_)) then
        window_ => window
    end if
    !else
    window = mpiwindow(sendcount,sendbuf,comm)
    !end if


    call getall_d(nproc,recvcounts,displs,window,nrecvbuf,recvbuf)

    if (.not. present(window_)) then
        call mpi_fenceandfree(window)
       !!! Synchronize the communication
       !!call mpi_win_fence(0, window, ierr)
       !!if (ierr/=0) then
       !!   call f_err_throw('Error in mpi_win_fence',&
       !!        err_id=ERR_MPI_WRAPPERS)  
       !!end if
       !!call mpi_win_free(window, ierr)
       !!if (ierr/=0) then
       !!   call f_err_throw('Error in mpi_win_fence',&
       !!        err_id=ERR_MPI_WRAPPERS)  
       !!end if
    end if

  end subroutine mpi_get_to_allgatherv_double
#endif

  
  subroutine mpiiallred_double(sendbuf, recvbuf, ncount, datatype, op, comm, request)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    ! Calling arguments
    integer,intent(in) :: ncount, datatype, op, comm
    double precision,intent(in) :: sendbuf
    double precision,intent(out) :: recvbuf
    integer,intent(out) :: request
    ! Local variables
    integer :: ierr

#ifdef HAVE_MPI3
    call mpi_iallreduce(sendbuf, recvbuf, ncount, datatype, op, comm, request, ierr)
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_IALLREDUCE occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
#else
    call mpi_allreduce(sendbuf, recvbuf, ncount, datatype, op, comm, ierr)
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_ALLREDUCE occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
    request = MPI_REQUEST_NULL
#endif

  end subroutine mpiiallred_double


  subroutine mpiialltoallv_double(sendbuf, sendcounts, senddspls, sendtype, &
             recvbuf, recvcounts, recvdspls, recvtype, comm, request)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    ! Calling arguments
    integer,intent(in) :: sendcounts, senddspls, sendtype, recvcounts, recvdspls, recvtype, comm
    double precision,intent(in) :: sendbuf
    double precision,intent(out) :: recvbuf
    integer,intent(out) :: request
    ! Local variables
    integer :: ierr

#ifdef HAVE_MPI3
    call mpi_ialltoallv(sendbuf, sendcounts, senddspls, sendtype, &
         recvbuf, recvcounts, recvdspls, recvtype, comm, request, ierr)
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_IALLTOALLV occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
#else
    call mpi_alltoallv(sendbuf, sendcounts, senddspls, sendtype, &
         recvbuf, recvcounts, recvdspls, recvtype, comm, ierr)
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_IALLTOALLV occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
    request = MPI_REQUEST_NULL
#endif

  end subroutine mpiialltoallv_double


  subroutine mpiwaitall(ncount, array_of_requests)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    ! Local variables
    integer,intent(in) :: ncount
    integer,dimension(ncount),intent(in) :: array_of_requests
    ! Local variables
    integer :: ierr

    call mpi_waitall(ncount, array_of_requests, MPI_STATUSES_IGNORE, ierr)
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_WAITALL occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if

  end subroutine mpiwaitall


  subroutine mpiwait(request)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    ! Local variables
    integer,intent(in) :: request
    ! Local variables
    integer :: ierr

    if (request /= MPI_REQUEST_NULL) then
       call mpi_wait(request, MPI_STATUSES_IGNORE, ierr)
       if (ierr/=0) then
          call f_err_throw('An error in calling to MPI_WAIT occured',&
               err_id=ERR_MPI_WRAPPERS)
       end if
    end if
  end subroutine mpiwait

end module wrapper_MPI


!> Routine to gather the clocks of all the instances of flib time module
subroutine gather_timings(ndata,nproc,mpi_comm,src,dest)
  use wrapper_MPI
  implicit none
  integer, intent(in) :: ndata !< number of categories of the array
  integer, intent(in) :: nproc,mpi_comm !< number of MPI tasks and communicator
  real(kind=8), dimension(ndata), intent(in) :: src !< total timings of the instance
  real(kind=8), dimension(ndata,nproc), intent(inout) :: dest !< gathered timings 
  call mpigather(sendbuf=src,recvbuf=dest,root=0,comm=mpi_comm)

end subroutine gather_timings

#ifdef HAVE_MPI2
!> used by get_to_allgatherv to pass the good addresses to the mpiget wrapper
subroutine getall_d(nproc,recvcounts,displs,window,nrecvbuffer,recvbuffer)
  use wrapper_MPI, only: mpiget, mpi_address_kind
  implicit none
  integer,intent(in) :: nproc,nrecvbuffer,window
  integer,dimension(0:nproc-1),intent(in) :: recvcounts,displs
  double precision,dimension(nrecvbuffer),intent(out) :: recvbuffer
  ! Local variables
  integer :: jproc, jcount, jst

  do jproc=0,nproc-1
     jcount=recvcounts(jproc)
     jst=displs(jproc)
     if (jcount>0) then
         call mpiget(recvbuffer(jst+1), jcount, jproc, int(0,kind=mpi_address_kind), window)
     end if
  end do

end subroutine getall_d
#endif


!> Activates the nesting for UNBLOCK_COMMS performance case
subroutine bigdft_open_nesting(num_threads)
  use wrapper_mpi
  implicit none
  integer, intent(in) :: num_threads
#ifdef HAVE_MPI_INIT_THREAD
  !$ call OMP_SET_NESTED(.true.) 
  !$ call OMP_SET_MAX_ACTIVE_LEVELS(2)
  !$ call OMP_SET_NUM_THREADS(num_threads)
#else
  integer :: idummy
  write(*,*)'BigDFT_open_nesting is not active!'
  !call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
  stop
  idummy=num_threads
#endif
end subroutine bigdft_open_nesting


!> Activates the nesting for UNBLOCK_COMMS performance case
subroutine bigdft_close_nesting(num_threads)
  use wrapper_mpi
  implicit none
  integer, intent(in) :: num_threads
#ifdef HAVE_MPI_INIT_THREAD
  !$ call OMP_SET_MAX_ACTIVE_LEVELS(1) !redundant
  !$ call OMP_SET_NESTED(.false.) 
  !$ call OMP_SET_NUM_THREADS(num_threads)
#else 
  integer :: idummy
  write(*,*)'BigDFT_close_nesting is not active!'
  stop
  !call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
  idummy=num_threads
#endif
end subroutine bigdft_close_nesting
