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
!  general module for embedded particles (eg 4 electrons in 1D equiv to
!  1 in 4D). Also calculate different densities on request.
!
#include "global.h"

module embedded_particles_m

  use datasets_m
  use global_m
  use grid_m
  use loct_m
  use loct_parser_m
  use messages_m

  implicit none

  private

  public :: embedded_particles_init,&
            embedded_particles_destroy,&
            embedded_particle_t

!==============================================================
!  container type for input vars concerning embedded particles
!==============================================================
type embedded_particle_t
   integer :: ndim_embedded              ! dimensionality of embedded space each
                                         !  particle lives in

   integer :: ntype_of_particle_embedded ! number of different types of particles
                                         !  embedded in MAX_DIM dimensional space

   integer :: nparticle_embedded         ! number of particles 

   integer :: ndensities_to_calculate

!   %block describe_particles_embedded
!   label1(char) | particletype1(integer) | mass1 | charge1
!   label2(char) | particletype2(integer) | mass2 | charge2
!   label3(char) | particletype3(integer) | mass3 | charge3
!   ...
!   nparticle_embedded lines (fixed)
!   %
   character(80), pointer :: labels_particles_embedded(:)

   integer, pointer :: particletype_embedded(:)

   FLOAT, pointer :: mass_particle_embedded(:)

   FLOAT, pointer :: charge_particle_embedded(:)


!   %block densitiestocalc
!   label1 |  particletokeep1(integer in [1:nparticle_embedded])
!   label2 |  particletokeep2(integer)
!   label3 |  particletokeep3(integer)
!   ...
!   however many lines wanted (up to nparticle_embedded)
!   %
   character(80), pointer :: labels_densities(:)

   integer, pointer :: particle_kept_densities(:)

end type embedded_particle_t


contains


!==============================================================
!  initialization function for embedded particles information
!==============================================================
subroutine embedded_particles_init (embeddedparticles,gr)

  implicit none

!args
  type(embedded_particle_t),intent(inout) :: embeddedparticles
  type(grid_t), intent(in) :: gr

!local vars
  integer :: ipart,ncols,nline
  type(block_t) :: blk

! source code

  call push_sub('states.embedded_particles_init')


! read in scalar dimensions
  !%Variable NDimEmbedded
  !%Type integer
  !%Section States
  !%Default -1
  !%Description
  !% Number of dimensions for embedded space 
  !% Full Ndim = NDimEmbedded*NParticleEmbedded
  !%
  !%End
  call loct_parse_int(datasets_check('NDimEmbedded'), -1, embeddedparticles%ndim_embedded)
  call messages_print_var_option(stdout, "NDimEmbedded", embeddedparticles%ndim_embedded)

  !%Variable NParticleEmbedded
  !%Type integer
  !%Section States
  !%Default 0
  !%Description
  !% Number of particles in embedded space 
  !% Full Ndim = NDimEmbedded*NParticleEmbedded
  !%
  !%End
  call loct_parse_int(datasets_check('NParticleEmbedded'), 0, embeddedparticles%nparticle_embedded)
  call messages_print_var_option(stdout, "NParticleEmbedded", embeddedparticles%nparticle_embedded)

  !%Variable NTypeParticleEmbedded
  !%Type integer
  !%Section States
  !%Default 1
  !%Description
  !% Number of different types of particles in embedded space 
  !%
  !%End
  call loct_parse_int(datasets_check('NTypeParticleEmbedded'), 1, embeddedparticles%ntype_of_particle_embedded)
  call messages_print_var_option(stdout, "NTypeParticleEmbedded", embeddedparticles%ntype_of_particle_embedded)
  if (embeddedparticles%ntype_of_particle_embedded > embeddedparticles%nparticle_embedded) then
     message(1) = ' Number of types of embedded particles should be <= Number of embedded particles'
     call write_fatal(1)
  end if

  if (embeddedparticles%ndim_embedded*embeddedparticles%nparticle_embedded /= gr%sb%dim) then
     message(1) = ' Number of embedded particles * dimension of embedded space must be = Ndim'
     call write_fatal(1)
  end if

! read in blocks
  !%Variable DescribeParticlesEmbedded
  !%Type block
  !%Section States
  !%Description
  !% Characterization of different particles embedded in NDIM dimensional space
  !%
  !% <tt>%DescribeParticlesEmbedded
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
  allocate (embeddedparticles%labels_particles_embedded(embeddedparticles%nparticle_embedded))
  allocate (embeddedparticles%particletype_embedded(embeddedparticles%nparticle_embedded))
  allocate (embeddedparticles%mass_particle_embedded(embeddedparticles%nparticle_embedded))
  allocate (embeddedparticles%charge_particle_embedded(embeddedparticles%nparticle_embedded))

! default all particles are electrons
  embeddedparticles%labels_particles_embedded='electron'
  embeddedparticles%particletype_embedded=1
  embeddedparticles%mass_particle_embedded=1.0d0
  embeddedparticles%charge_particle_embedded=1.0d0
    

  if(loct_parse_block(datasets_check('DescribeParticlesEmbedded'), blk)==0) then
      ncols = loct_parse_block_cols(blk, 0)
      if(ncols /= 4 ) then
        call input_error("DescribeParticlesEmbedded")
      end if
      nline=loct_parse_block_n(blk)
      if (nline /= embeddedparticles%nparticle_embedded) then
        call input_error("DescribeParticlesEmbedded")
      end if

      do ipart=1,embeddedparticles%nparticle_embedded
        call loct_parse_block_string(blk, ipart-1, 0, embeddedparticles%labels_particles_embedded(ipart))
        call loct_parse_block_int   (blk, ipart-1, 1, embeddedparticles%particletype_embedded(ipart))
        call loct_parse_block_float (blk, ipart-1, 2, embeddedparticles%mass_particle_embedded(ipart))
        call loct_parse_block_float (blk, ipart-1, 3, embeddedparticles%charge_particle_embedded(ipart))

        write (message(1),'(a,a)') 'labels_particles_embedded = ', embeddedparticles%labels_particles_embedded(ipart)
        write (message(2),'(a,i6)') 'particletype_embedded = ', embeddedparticles%particletype_embedded(ipart)
        write (message(3),'(a,E20.10)') 'mass_particle_embedded = ', embeddedparticles%mass_particle_embedded(ipart)
        write (message(4),'(a,E20.10)') 'charge_particle_embedded = ', embeddedparticles%charge_particle_embedded(ipart)
        call write_info(4)
      end do
      call loct_parse_block_end(blk)
  end if


  !%Variable DensitiestoCalc
  !%Type block
  !%Section States
  !%Description
  !% choice of which particle densities will be calculated and output, in the
  !%  embedded particles scheme
  !%
  !% <tt>%DensitiestoCalc
  !% <br>&nbsp;&nbsp; proton   | 1
  !% <br>&nbsp;&nbsp; electron | 2
  !% <br>%</tt>
  !%
  !% would ask octopus to calculate the density corresponding to the 1st
  !% particle (whose coordinates correspond to dimensions 1 to ndim_embedded),
  !% which is an electron, then that corresponding to the 2nd particle
  !% (dimensions ndim_embedded+1 to 2*ndim_embedded)
  !%
  !%End
  embeddedparticles%ndensities_to_calculate=0

  if(loct_parse_block(datasets_check('DensitiestoCalc'), blk)==0) then
      ncols = loct_parse_block_cols(blk, 0)
      if(ncols /= 2 ) then
        call input_error("DensitiestoCalc")
      end if
      embeddedparticles%ndensities_to_calculate=loct_parse_block_n(blk)
      if (embeddedparticles%ndensities_to_calculate < 0 .or. &
          embeddedparticles%ndensities_to_calculate > embeddedparticles%nparticle_embedded) then
        call input_error("DensitiestoCalc")
      end if

      allocate (embeddedparticles%labels_densities(embeddedparticles%ndensities_to_calculate))
      allocate (embeddedparticles%particle_kept_densities(embeddedparticles%ndensities_to_calculate))

      do ipart=1,embeddedparticles%ndensities_to_calculate
        call loct_parse_block_string(blk, ipart-1, 0, embeddedparticles%labels_densities(ipart))
        call loct_parse_block_int(blk, ipart-1, 1, embeddedparticles%particle_kept_densities(ipart))

        write (message(1),'(a,a)') 'labels_densities = ', embeddedparticles%labels_densities(ipart)
        write (message(2),'(a,i6)') 'particle_kept_densities = ', embeddedparticles%particle_kept_densities(ipart)
        call write_info(2)
      end do
      call loct_parse_block_end(blk)
  end if

  call pop_sub()

end subroutine embedded_particles_init


subroutine embedded_particles_destroy (embeddedparticles)

  implicit none

!args
  type(embedded_particle_t),intent(inout) :: embeddedparticles

! source

  call push_sub('states.embedded_particles_destroy')

  if (associated(embeddedparticles%labels_particles_embedded)) deallocate(embeddedparticles%labels_particles_embedded)
  if (associated(embeddedparticles%particletype_embedded)) deallocate(embeddedparticles%particletype_embedded)
  if (associated(embeddedparticles%mass_particle_embedded)) deallocate(embeddedparticles%mass_particle_embedded)
  if (associated(embeddedparticles%charge_particle_embedded)) deallocate(embeddedparticles%charge_particle_embedded)

  if (associated(embeddedparticles%labels_densities)) deallocate(embeddedparticles%labels_densities)
  if (associated(embeddedparticles%particle_kept_densities)) deallocate(embeddedparticles%particle_kept_densities)

  call pop_sub()

end subroutine embedded_particles_destroy

end module embedded_particles_m
