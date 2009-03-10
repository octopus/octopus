!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!
!  general module for modelMB particles (eg 4 electrons in 1D equiv to
!  1 in 4D). Also calculate different densities on request.
!
#include "global.h"

module modelMB_particles_m

  use datasets_m
  use global_m
  use grid_m
  use loct_m
  use loct_parser_m
  use messages_m

  implicit none

  private

  public :: modelMB_particles_init,&
            modelMB_particles_destroy,&
            modelMB_particle_t

!==============================================================
!  container type for input vars concerning modelMB particles
!==============================================================
type modelMB_particle_t
   integer :: ndim_modelMB              ! dimensionality of modelMB space each
                                         !  particle lives in

   integer :: ntype_of_particle_modelMB ! number of different types of particles
                                         !  modelMB in MAX_DIM dimensional space

   integer :: nparticle_modelMB         ! number of particles 

   integer :: ndensities_to_calculate

!   %block describe_particles_modelMB
!   label1(char) | particletype1(integer) | mass1 | charge1
!   label2(char) | particletype2(integer) | mass2 | charge2
!   label3(char) | particletype3(integer) | mass3 | charge3
!   ...
!   nparticle_modelMB lines (fixed)
!   %
   character(80), pointer :: labels_particles_modelMB(:)

   integer, pointer :: particletype_modelMB(:)

   FLOAT, pointer :: mass_particle_modelMB(:)

   FLOAT, pointer :: charge_particle_modelMB(:)


!   %block densitiestocalc
!   label1 |  particletokeep1(integer in [1:nparticle_modelMB])
!   label2 |  particletokeep2(integer)
!   label3 |  particletokeep3(integer)
!   ...
!   however many lines wanted (up to nparticle_modelMB)
!   %
   character(80), pointer :: labels_densities(:)

   integer, pointer :: particle_kept_densities(:)

end type modelMB_particle_t


contains


!==============================================================
!  initialization function for modelMB particles information
!==============================================================
subroutine modelMB_particles_init (modelMBparticles,gr)

  implicit none

!args
  type(modelMB_particle_t),intent(inout) :: modelMBparticles
  type(grid_t), intent(in) :: gr

!local vars
  integer :: ipart,ncols,nline
  type(block_t) :: blk

! source code

  call push_sub('states.modelMB_particles_init')


! read in scalar dimensions
  !%Variable NDimModelMB
  !%Type integer
  !%Section States
  !%Default -1
  !%Description
  !% Number of dimensions for modelMB space 
  !% Full Ndim = NDimModelMB*NParticleModelMB
  !%
  !%End
  call loct_parse_int(datasets_check('NDimModelMB'), -1, modelMBparticles%ndim_modelMB)
  call messages_print_var_option(stdout, "NDimModelMB", modelMBparticles%ndim_modelMB)

  !%Variable NParticleModelMB
  !%Type integer
  !%Section States
  !%Default 0
  !%Description
  !% Number of particles in modelMB space 
  !% Full Ndim = NDimModelMB*NParticleModelMB
  !%
  !%End
  call loct_parse_int(datasets_check('NParticleModelMB'), 0, modelMBparticles%nparticle_modelMB)
  call messages_print_var_option(stdout, "NParticleModelMB", modelMBparticles%nparticle_modelMB)

  !%Variable NTypeParticleModelMB
  !%Type integer
  !%Section States
  !%Default 1
  !%Description
  !% Number of different types of particles in modelMB space 
  !%
  !%End
  call loct_parse_int(datasets_check('NTypeParticleModelMB'), 1, modelMBparticles%ntype_of_particle_modelMB)
  call messages_print_var_option(stdout, "NTypeParticleModelMB", modelMBparticles%ntype_of_particle_modelMB)
  if (modelMBparticles%ntype_of_particle_modelMB > modelMBparticles%nparticle_modelMB) then
     message(1) = ' Number of types of modelMB particles should be <= Number of modelMB particles'
     call write_fatal(1)
  end if

  if (modelMBparticles%ndim_modelMB*modelMBparticles%nparticle_modelMB /= gr%sb%dim) then
     message(1) = ' Number of modelMB particles * dimension of modelMB space must be = Ndim'
     call write_fatal(1)
  end if

! read in blocks
  !%Variable DescribeParticlesModelMB
  !%Type block
  !%Section States
  !%Description
  !% Characterization of different modelMB particles in NDIM dimensional space
  !%
  !% <tt>%DescribeParticlesModelMB
  !% <br>&nbsp;&nbsp; proton   | 1 | 1800. | 1.
  !% <br>&nbsp;&nbsp; proton   | 1 | 1800. | 1.
  !% <br>&nbsp;&nbsp; electron | 2 | 1.    | 1.
  !% <br>%</tt>
  !%
  !% would tell octopus that there are presently 3 particles, called proton, proton,
  !% and electron, with types 1, 1, and 2, and corresponding masses and charges.
  !% The label and charge are presently only for informational purposes and
  !% are not checked or used in octopus. The interaction has to take the
  !% actual charge into account.
  !%
  !%End
! allocate stuff
  allocate (modelMBparticles%labels_particles_modelMB(modelMBparticles%nparticle_modelMB))
  allocate (modelMBparticles%particletype_modelMB(modelMBparticles%nparticle_modelMB))
  allocate (modelMBparticles%mass_particle_modelMB(modelMBparticles%nparticle_modelMB))
  allocate (modelMBparticles%charge_particle_modelMB(modelMBparticles%nparticle_modelMB))

! default all particles are electrons
  modelMBparticles%labels_particles_modelMB='electron'
  modelMBparticles%particletype_modelMB=1
  modelMBparticles%mass_particle_modelMB=1.0d0
  modelMBparticles%charge_particle_modelMB=1.0d0
    

  if(loct_parse_block(datasets_check('DescribeParticlesModelMB'), blk)==0) then
      ncols = loct_parse_block_cols(blk, 0)
      if(ncols /= 4 ) then
        call input_error("DescribeParticlesModelMB")
      end if
      nline=loct_parse_block_n(blk)
      if (nline /= modelMBparticles%nparticle_modelMB) then
        call input_error("DescribeParticlesModelMB")
      end if

      do ipart=1,modelMBparticles%nparticle_modelMB
        call loct_parse_block_string(blk, ipart-1, 0, modelMBparticles%labels_particles_modelMB(ipart))
        call loct_parse_block_int   (blk, ipart-1, 1, modelMBparticles%particletype_modelMB(ipart))
        call loct_parse_block_float (blk, ipart-1, 2, modelMBparticles%mass_particle_modelMB(ipart))
        call loct_parse_block_float (blk, ipart-1, 3, modelMBparticles%charge_particle_modelMB(ipart))

        write (message(1),'(a,a)') 'labels_particles_modelMB = ', modelMBparticles%labels_particles_modelMB(ipart)
        write (message(2),'(a,i6)') 'particletype_modelMB = ', modelMBparticles%particletype_modelMB(ipart)
        write (message(3),'(a,E20.10)') 'mass_particle_modelMB = ', modelMBparticles%mass_particle_modelMB(ipart)
        write (message(4),'(a,E20.10)') 'charge_particle_modelMB = ', modelMBparticles%charge_particle_modelMB(ipart)
        call write_info(4)
      end do
      call loct_parse_block_end(blk)
  end if


  !%Variable DensitiestoCalc
  !%Type block
  !%Section States
  !%Description
  !% choice of which particle densities will be calculated and output, in the
  !%  modelMB particles scheme
  !%
  !% <tt>%DensitiestoCalc
  !% <br>&nbsp;&nbsp; proton   | 1
  !% <br>&nbsp;&nbsp; electron | 2
  !% <br>%</tt>
  !%
  !% would ask octopus to calculate the density corresponding to the 1st
  !% particle (whose coordinates correspond to dimensions 1 to ndim_modelMB),
  !% which is an electron, then that corresponding to the 2nd particle
  !% (dimensions ndim_modelMB+1 to 2*ndim_modelMB)
  !%
  !%End
  modelMBparticles%ndensities_to_calculate=0

  if(loct_parse_block(datasets_check('DensitiestoCalc'), blk)==0) then
      ncols = loct_parse_block_cols(blk, 0)
      if(ncols /= 2 ) then
        call input_error("DensitiestoCalc")
      end if
      modelMBparticles%ndensities_to_calculate=loct_parse_block_n(blk)
      if (modelMBparticles%ndensities_to_calculate < 0 .or. &
          modelMBparticles%ndensities_to_calculate > modelMBparticles%nparticle_modelMB) then
        call input_error("DensitiestoCalc")
      end if

      allocate (modelMBparticles%labels_densities(modelMBparticles%ndensities_to_calculate))
      allocate (modelMBparticles%particle_kept_densities(modelMBparticles%ndensities_to_calculate))

      do ipart=1,modelMBparticles%ndensities_to_calculate
        call loct_parse_block_string(blk, ipart-1, 0, modelMBparticles%labels_densities(ipart))
        call loct_parse_block_int(blk, ipart-1, 1, modelMBparticles%particle_kept_densities(ipart))

        write (message(1),'(a,a)') 'labels_densities = ', modelMBparticles%labels_densities(ipart)
        write (message(2),'(a,i6)') 'particle_kept_densities = ', modelMBparticles%particle_kept_densities(ipart)
        call write_info(2)
      end do
      call loct_parse_block_end(blk)
  end if

  call pop_sub()

end subroutine modelMB_particles_init


subroutine modelMB_particles_destroy (modelMBparticles)

  implicit none

!args
  type(modelMB_particle_t),intent(inout) :: modelMBparticles

! source

  call push_sub('states.modelMB_particles_destroy')

  if (associated(modelMBparticles%labels_particles_modelMB)) deallocate(modelMBparticles%labels_particles_modelMB)
  if (associated(modelMBparticles%particletype_modelMB)) deallocate(modelMBparticles%particletype_modelMB)
  if (associated(modelMBparticles%mass_particle_modelMB)) deallocate(modelMBparticles%mass_particle_modelMB)
  if (associated(modelMBparticles%charge_particle_modelMB)) deallocate(modelMBparticles%charge_particle_modelMB)

  if (associated(modelMBparticles%labels_densities)) deallocate(modelMBparticles%labels_densities)
  if (associated(modelMBparticles%particle_kept_densities)) deallocate(modelMBparticles%particle_kept_densities)

  call pop_sub()

end subroutine modelMB_particles_destroy

end module modelMB_particles_m
