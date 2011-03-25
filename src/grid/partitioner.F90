!! Copyright (C) 2010 X. Andrade
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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: partitioner.F90 6396 2010-03-26 08:51:58Z mjv500 $

#include "global.h"

module partitioner_m
  use batch_m
  use c_pointer_m
  use datasets_m
  use global_m
  use iihash_m
  use index_m
  use io_m
  use mesh_m
  use messages_m
  use mpi_m
  use loct_math_m
  use partition_m
  use parser_m
  use profiling_m
  use stencil_m

  implicit none

  private

  public ::                    &
    partitioner_t,             &
    partitioner_init,          &
    partitioner_perform,       &
    partitioner_end

  type partitioner_t
    type(c_ptr)                :: rng
    integer                    :: npop
    integer                    :: nsteps
  end type partitioner_t

contains

  subroutine partitioner_init(this)
    type(partitioner_t), intent(out) :: this

    call loct_ran_init(this%rng)

    !%Variable MeshPartitionGAPopulation
    !%Type integer
    !%Default 30
    !%Section Execution::Parallelization
    !%Description
    !% The size of the population used for the genetic algorithm used
    !% to optimize the mesh partition. The default is 30.
    !%End
    call parse_integer(datasets_check('MeshPartitionGAPopulation'), 30, this%npop)

    !%Variable MeshPartitionGAMaxSteps
    !%Type integer
    !%Default 1000
    !%Section Execution::Parallelization
    !%Description
    !% The number of steps performed for the genetic algorithm used
    !% to optimize the mesh partition. The default is 1000.
    !%End
    call parse_integer(datasets_check('MeshPartitionGAMaxSteps'), 1000, this%nsteps)

  end subroutine partitioner_init
  
  ! ----------------------------------------------------
  
  subroutine partitioner_perform(this, mesh, stencil, best)
    type(partitioner_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(stencil_t),     intent(in)    :: stencil
    type(partition_t),   intent(inout) :: best

    integer :: ipop, igen, childnum, iparent(1:2), ibest, istep
    type(partition_t), pointer :: parents(:), childs(:), ptmp(:)
    FLOAT, allocatable :: fitness(:), probability(:)
    FLOAT :: fsum, psum, random
    FLOAT, parameter :: mutation_rate = 0.2
    FLOAT :: fbest

    message(1) = "Info: Performing a genetic algorithm optimization of the mesh partitioning"
    call messages_info(1)

    SAFE_ALLOCATE(parents(1:this%npop))
    SAFE_ALLOCATE(childs(1:this%npop))
    
    do ipop = 1, this%npop
      call partition_init(parents(ipop), mesh)
      call partition_init(childs(ipop), mesh)
      call partition_randomize(parents(ipop), this%rng)
    end do

    SAFE_ALLOCATE(fitness(1:this%npop))
    SAFE_ALLOCATE(probability(1:this%npop))
    
    fbest = -M_ONE

    do igen = 1, this%nsteps
      call get_fitness()

      if(abs(maxval(fitness) - fbest) > M_EPSILON) then
        fbest = maxval(fitness)
        write(message(1), '(a,i5,a,e10.4)') "      GA generation: ", igen, " fitness ", fbest 
        call messages_info(1)
      end if

      fsum = sum(fitness)
      
      psum = M_ZERO
      do ipop = 1, this%npop
        probability(ipop) = psum + fitness(ipop)/fsum
        psum = psum + fitness(ipop)/fsum
      end do

      ! elitism: copy the two best partitions
      ibest = maxloc(fitness, dim = 1)
      call partition_copy(parents(ibest), childs(1))
      fitness(ibest) = M_ZERO
      ibest = maxloc(fitness, dim = 1)
      call partition_copy(parents(ibest), childs(2))      
      childnum = 2

      do while(childnum < this%npop)
        
        ! get two parents
        do istep = 1, 2
          random = loct_ran_flat(this%rng, M_ZERO, M_ONE)

          iparent(istep) = this%npop ! we can avoid the last iteration
          do ipop = 1, this%npop - 1
            if(random < probability(ipop)) then 
              iparent(istep) = ipop
              exit
            end if
          end do
        end do
        
        call partition_crossover(parents(iparent(1)), parents(iparent(2)), this%rng, childs(childnum + 1), childs(childnum + 2))

        INCR(childnum, 2)

      end do

      do ipop = 2, this%npop
        random = loct_ran_flat(this%rng, M_ZERO, M_ONE)
        if(random < mutation_rate) call partition_mutate(childs(ipop), this%rng, mesh, stencil) ! my son is a mutant!!!
      end do

      ! now childs become parents (we recycle the parents array)
      ptmp => parents
      parents => childs
      childs => ptmp

    end do
    
    call partition_copy(parents(1), best)

    SAFE_DEALLOCATE_A(fitness)
    SAFE_DEALLOCATE_P(parents)
    SAFE_DEALLOCATE_P(childs)
  contains 

    subroutine get_fitness()
      do ipop = 1, this%npop
        call partition_build(parents(ipop), mesh, stencil)
        fitness(ipop) = partition_quality(parents(ipop))
      end do
      ipop = maxloc(fitness, dim = 1)
    end subroutine get_fitness

  end subroutine partitioner_perform

  ! ------------------------------------------

  subroutine partitioner_end(this)
    type(partitioner_t), intent(inout) :: this

  end subroutine partitioner_end

end module partitioner_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
