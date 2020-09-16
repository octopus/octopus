!! Copyright (C) 2020 Heiko Appel, Kevin Lively
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

module replicas_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use multisystem_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use system_oct_m
  use system_factory_abst_oct_m
  implicit none

  private
  public ::               &
    replicas_t,           &
    replicas_init

  integer, public, parameter ::        &
    UNIFORM_REPLICA             =  1,  &
    GAUSS_REPLICA               =  2,  &
    INPUT_REPLICA               =  3

  type, extends(multisystem_t), abstract :: replicas_t

    integer :: n_replicas = 0
    integer :: replica_distribution = 1
    FLOAT   :: distribution_width = CNST(0.01)
 
  contains
    procedure :: create_system => replicas_create_system
  end type replicas_t

!  interface replicas_t
!    procedure replicas_constructor
!  end interface replicas_t

contains

  recursive subroutine replicas_init(this, namespace, factory)
    class(replicas_t),      intent(inout) :: this
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory

    character(len=128) :: system_name, replica_name
    integer :: system_replicas_default, jj

    PUSH_SUB(replicas_init)

    !%Variable SystemReplicas
    !%Type integer
    !%Section Time-Dependent::Propagation
    !%Description
    !% Number of system replicas
    !%End

    system_replicas_default = 0
    call parse_variable(this%namespace, 'SystemReplicas', system_replicas_default, this%n_replicas)
    write(message(1), '(a,a,a,i6)') 'Namespace: ', trim(this%namespace%get()), ' SystemReplicas:', &
            this%n_replicas
    call messages_info(1)

    call multisystem_init(this, namespace, factory)

!    call this%supported_interactions%add(COPY_REPLICA_DATA)


    !%Variable ReplicaDistribution
    !%Type integer
    !%Section Time-Dependent::Propagation
    !%Description
    !% The distribution to sample the initial conditions of the system replicas from.
    !% If none specific then files containing the initial conditions of each replica must be specified (?)
    !%Option uniform 1
    !%Option gauss 2
    !%Option input 3
    !%End
    !System Replica Distribution
    !For now this just gets the centroid
    call parse_variable(this%namespace, 'ReplicaDistribution', this%replica_distribution, this%replica_distribution)

    call messages_print_var_value(stdout, 'ReplicaDistribution is ', this%replica_distribution)

    !%Variable ReplicaDistributionWidth
    !%Type float
    !%Section Time-Dependent::Propagation
    !%Description
    !% The distribution_width of the distribution distributing the replica
    !% For the gaussian distribtion, this is the standard deviation.
    !%End
    call parse_variable(namespace, 'ReplicaDistributionWidth', this%distribution_width, this%distribution_width)
    call messages_print_var_value(stdout, 'ReplicaDistributionWidth', &
                                  this%distribution_width)


    POP_SUB(replicas_init)
  end subroutine replicas_init

  ! ---------------------------------------------------------------------------------------
  recursive subroutine replicas_create_system(this, system_name, system_type, isys, factory)
    class(replicas_t),      intent(inout) :: this
    character(len=128),           intent(in) :: system_name
    integer,                      intent(in) :: system_type
    integer,                      intent(in) :: isys
    class(system_factory_abst_t), intent(in) :: factory

    integer :: jj
!    type(system_iterator_t) :: iter
    class(system_t), pointer :: sys !, other
    character(len=128) :: replica_name

    PUSH_SUB(replicas_create_system)

    do jj = 0, this%n_replicas
      ! Create system
      if (jj .ne. 0) then
        write(replica_name,'(a,a,i8.8)') trim(system_name), '-', jj
      end if
      ! Create folder to store system files.
      ! Needs to be done before creating the system as this in turn might create subfolders.
      call io_mkdir(replica_name, namespace=this%namespace)

      sys => factory%create(this%namespace, system_name, system_type)
      if (.not. associated(sys)) then
        call messages_input_error(this%namespace, factory%block_name(), 'Unknown system type.')
      end if

      ! Add system to list of systems
      call this%list%add(sys)
    end do

    ! ToDo: needs to be revisited
    ! Check that the system is unique
!    call iter%start(this%list)
!    do while (iter%has_next())
!      other => iter%get_next()
!      if (sys%namespace == other%namespace) then
!        call messages_input_error(this%namespace, factory%block_name(), 'Duplicated system in multisystem', &
!          row=isys-1, column=0)
!      end if
!    end do

    POP_SUB(replicas_create_system)
  end subroutine replicas_create_system

end module replicas_oct_m
