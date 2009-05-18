!! Copyright (C) 2009 N. Helbig and M. Verstraete
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
!  general module for modelmb particles (eg 4 electrons in 1D equiv to
!  1 in 4D). Also calculate different densities on request.
!
#include "global.h"

module modelmb_particles_m

  use datasets_m
  use global_m
  use grid_m
  use hypercube_m
  use loct_m
  use loct_parser_m
  use messages_m
  use mesh_m
  use profiling_m

  implicit none

  private

  public :: modelmb_particles_nullify,   &
            modelmb_particles_init,      &
            modelmb_particles_end,       &
            modelmb_particles_copy,      &
            modelmb_copy_masses,         &
            modelmb_particle_t

!==============================================================
!  container type for input vars concerning modelmb particles
!==============================================================
type modelmb_particle_t
   integer :: ndim              ! dimensionality of modelmb space each
                                         !  particle lives in

   integer :: ntype_of_particle ! number of different types of particles
                                        !  modelmb in MAX_DIM dimensional space
   integer :: max_particles_per_type    ! max number of different particle

   integer :: nparticle         ! number of particles 

   integer :: ndensities_to_calculate

   !   %block describe_particles_modelmb
   !   label1(char) | particletype1(integer) | mass1 | charge1 | fermion or boson or anyon
   !   label2(char) | particletype2(integer) | mass2 | charge2 | fermion or boson or anyon
   !   label3(char) | particletype3(integer) | mass3 | charge3 | fermion or boson or anyon
   !   ...
   !   nparticle_modelmb lines (fixed)
   !   %
   character(80), pointer :: labels_particles(:)
   character(80), pointer :: bosonfermion(:)

   integer, pointer :: particletype(:)
   integer, pointer :: nparticles_per_type(:)

   integer, pointer :: exchange_symmetry(:,:,:) ! (max_particles_per_type**2, ntype_of_particle)

   FLOAT, pointer :: mass_particle(:)

   FLOAT, pointer :: charge_particle(:)


   !   %block densitiestocalc
   !   label1 |  particletokeep1(integer in [1:nparticle])
   !   label2 |  particletokeep2(integer)
   !   label3 |  particletokeep3(integer)
   !   ...
   !   however many lines wanted (up to nparticle_modelmb)
   !   %
   character(80), pointer :: labels_densities(:)

   integer, pointer :: particle_kept_densities(:)

end type modelmb_particle_t

contains


subroutine modelmb_particles_nullify(this)
  type(modelmb_particle_t), intent(inout) :: this

  nullify(this%labels_particles)
  nullify(this%particletype)
  nullify(this%nparticles_per_type)
  nullify(this%exchange_symmetry)
  nullify(this%bosonfermion)
  nullify(this%mass_particle)
  nullify(this%charge_particle)
  nullify(this%labels_densities)
  nullify(this%particle_kept_densities)
end subroutine modelmb_particles_nullify


!==============================================================
!  initialization function for modelmb particles information
!==============================================================
subroutine modelmb_particles_init (modelmbparticles,gr)
  type(modelmb_particle_t), intent(inout) :: modelmbparticles
  type(grid_t),             intent(in)    :: gr

  integer :: ipart, ncols, nline, itmp, jtmp
  type(block_t) :: blk

  call push_sub('states.modelmb_particles_init')

  ! read in scalar dimensions

  !%Variable NDimModelmb
  !%Type integer
  !%Section States
  !%Default -1
  !%Description
  !% Number of dimensions for modelmb space 
  !% Full Ndim = NDimModelmb*NParticleModelmb
  !%
  !%End
  call loct_parse_int(datasets_check('NDimModelmb'), gr%sb%dim, modelmbparticles%ndim)
  call messages_print_var_option(stdout, "NDimModelmb", modelmbparticles%ndim)

  !%Variable NParticleModelmb
  !%Type integer
  !%Section States
  !%Default 0
  !%Description
  !% Number of particles in modelmb space 
  !% Full Ndim = NDimModelmb*NParticleModelmb
  !%
  !%End
  call loct_parse_int(datasets_check('NParticleModelmb'), 1, modelmbparticles%nparticle)
  call messages_print_var_option(stdout, "NParticleModelmb", modelmbparticles%nparticle)

  !%Variable NTypeParticleModelmb
  !%Type integer
  !%Section States
  !%Default 1
  !%Description
  !% Number of different types of particles in modelmb space 
  !%
  !%End
  call loct_parse_int(datasets_check('NTypeParticleModelmb'), 1, modelmbparticles%ntype_of_particle)
  call messages_print_var_option(stdout, "NTypeParticleModelmb", modelmbparticles%ntype_of_particle)
  if (modelmbparticles%ntype_of_particle > modelmbparticles%nparticle) then
     write (message(1), '(a,2I6)') ' Number of types of modelmb particles should be <= Number of modelmb particles ', &
             modelmbparticles%ntype_of_particle, modelmbparticles%nparticle
     call write_fatal(1)
  end if

  if (modelmbparticles%ndim*modelmbparticles%nparticle /= gr%sb%dim) then
     message(1) = ' Number of modelmb particles * dimension of modelmb space must be = Ndim'
     call write_fatal(1)
  end if

  ! read in blocks

  !%Variable DescribeParticlesModelmb
  !%Type block
  !%Section States
  !%Description
  !% Characterization of different modelmb particles in gr%mesh%sb%dim dimensional space
  !%
  !% <tt>%DescribeParticlesModelmb
  !% <br>&nbsp;&nbsp; proton   | 1 | 1800. | 1. | fermion
  !% <br>&nbsp;&nbsp; proton   | 1 | 1800. | 1. | fermion
  !% <br>&nbsp;&nbsp; electron | 2 | 1.    | 1. | fermion
  !% <br>%</tt>
  !%
  !% would tell octopus that there are presently 3 particles, called proton, proton,
  !% and electron, with types 1, 1, and 2, and corresponding masses and charges.
  !% All particles should be fermions, and this can be later enforced on the spatial
  !% part of the wave functions.
  !% The label and charge are presently only for informational purposes and
  !% are not checked or used in octopus. The interaction has to take the
  !% actual charge into account.
  !%
  !%End
  ! allocate stuff
  SAFE_ALLOCATE (modelmbparticles%labels_particles(1:modelmbparticles%nparticle))
  SAFE_ALLOCATE (modelmbparticles%particletype(1:modelmbparticles%nparticle))
  SAFE_ALLOCATE (modelmbparticles%mass_particle(1:modelmbparticles%nparticle))
  SAFE_ALLOCATE (modelmbparticles%charge_particle(1:modelmbparticles%nparticle))
  SAFE_ALLOCATE (modelmbparticles%bosonfermion(1:modelmbparticles%nparticle))
  SAFE_ALLOCATE (modelmbparticles%nparticles_per_type(1:modelmbparticles%ntype_of_particle))

  ! default all particles are electrons
  modelmbparticles%labels_particles = 'electron'
  modelmbparticles%particletype = 1
  modelmbparticles%mass_particle = 1.0d0
  modelmbparticles%charge_particle = 1.0d0
  modelmbparticles%bosonfermion = 'fermion'
    

  if(loct_parse_block(datasets_check('DescribeParticlesModelmb'), blk)==0) then
      ncols = loct_parse_block_cols(blk, 0)
      if(ncols /= 5 ) then
        call input_error("DescribeParticlesModelmb")
      end if
      nline = loct_parse_block_n(blk)
      if (nline /= modelmbparticles%nparticle) then
        call input_error("DescribeParticlesModelmb")
      end if

      do ipart = 1,modelmbparticles%nparticle
        call loct_parse_block_string(blk, ipart-1, 0, modelmbparticles%labels_particles(ipart))
        call loct_parse_block_int   (blk, ipart-1, 1, modelmbparticles%particletype(ipart))
        call loct_parse_block_float (blk, ipart-1, 2, modelmbparticles%mass_particle(ipart))
        call loct_parse_block_float (blk, ipart-1, 3, modelmbparticles%charge_particle(ipart))
        call loct_parse_block_string(blk, ipart-1, 4, modelmbparticles%bosonfermion(ipart))

        write (message(1),'(a,a)') 'labels_particles_modelmb = ', modelmbparticles%labels_particles(ipart)
        write (message(2),'(a,i6)') 'particletype_modelmb = ', modelmbparticles%particletype(ipart)
        write (message(3),'(a,E20.10)') 'mass_particle_modelmb = ', modelmbparticles%mass_particle(ipart)
        write (message(4),'(a,E20.10)') 'charge_particle_modelmb = ', modelmbparticles%charge_particle(ipart)
        write (message(5),'(a,a)') 'bosonfermion = ', modelmbparticles%bosonfermion(ipart)
        call write_info(5)
      end do
      call loct_parse_block_end(blk)

  end if
    
  modelmbparticles%nparticles_per_type = 0
  do ipart = 1, modelmbparticles%nparticle
    modelmbparticles%nparticles_per_type(modelmbparticles%particletype(ipart)) = & 
      & modelmbparticles%nparticles_per_type(modelmbparticles%particletype(ipart)) + 1
  enddo

  modelmbparticles%max_particles_per_type = maxval(modelmbparticles%nparticles_per_type)
  itmp = modelmbparticles%max_particles_per_type
  jtmp = modelmbparticles%ntype_of_particle
  SAFE_ALLOCATE (modelmbparticles%exchange_symmetry(1:itmp, 1:itmp, 1:jtmp))
  modelmbparticles%exchange_symmetry = 0

  !%Variable DensitiestoCalc
  !%Type block
  !%Section States
  !%Description
  !% choice of which particle densities will be calculated and output, in the
  !%  modelmb particles scheme
  !%
  !% <tt>%DensitiestoCalc
  !% <br>&nbsp;&nbsp; proton   | 1
  !% <br>&nbsp;&nbsp; electron | 2
  !% <br>%</tt>
  !%
  !% would ask octopus to calculate the density corresponding to the 1st
  !% particle (whose coordinates correspond to dimensions 1 to ndim_modelmb),
  !% which is an electron, then that corresponding to the 2nd particle
  !% (dimensions ndim_modelmb+1 to 2*ndim_modelmb)
  !%
  !%End
  modelmbparticles%ndensities_to_calculate = 0

  if(loct_parse_block(datasets_check('DensitiestoCalc'), blk)==0) then
    ncols = loct_parse_block_cols(blk, 0)
    if(ncols /= 2 ) then
      call input_error("DensitiestoCalc")
    end if
    modelmbparticles%ndensities_to_calculate = loct_parse_block_n(blk)
    if (modelmbparticles%ndensities_to_calculate < 0 .or. &
         modelmbparticles%ndensities_to_calculate > modelmbparticles%nparticle) then
      call input_error("DensitiestoCalc")
    end if
    SAFE_ALLOCATE (modelmbparticles%labels_densities(1:modelmbparticles%ndensities_to_calculate))
    SAFE_ALLOCATE (modelmbparticles%particle_kept_densities(1:modelmbparticles%ndensities_to_calculate))
    do ipart = 1,modelmbparticles%ndensities_to_calculate
      call loct_parse_block_string(blk, ipart-1, 0, modelmbparticles%labels_densities(ipart))
      call loct_parse_block_int(blk, ipart-1, 1, modelmbparticles%particle_kept_densities(ipart))
      
      write (message(1),'(a,a)') 'labels_densities = ', modelmbparticles%labels_densities(ipart)
      write (message(2),'(a,i6)') 'particle_kept_densities = ', modelmbparticles%particle_kept_densities(ipart)
      call write_info(2)
    end do
    call loct_parse_block_end(blk)
  else
    nullify(modelmbparticles%labels_densities)
    nullify(modelmbparticles%particle_kept_densities)
  end if

  call pop_sub()

end subroutine modelmb_particles_init


subroutine modelmb_particles_end (modelmbparticles)
  type(modelmb_particle_t),intent(inout) :: modelmbparticles

  call push_sub('states.modelmb_particles_end')

  SAFE_DEALLOCATE_P(modelmbparticles%labels_particles)
  SAFE_DEALLOCATE_P(modelmbparticles%particletype)
  SAFE_DEALLOCATE_P(modelmbparticles%mass_particle)
  SAFE_DEALLOCATE_P(modelmbparticles%charge_particle)
  SAFE_DEALLOCATE_P(modelmbparticles%nparticles_per_type)
  SAFE_DEALLOCATE_P(modelmbparticles%exchange_symmetry)
  SAFE_DEALLOCATE_P(modelmbparticles%bosonfermion)

  SAFE_DEALLOCATE_P(modelmbparticles%labels_densities)
  SAFE_DEALLOCATE_P(modelmbparticles%particle_kept_densities)

  call pop_sub()

end subroutine modelmb_particles_end

subroutine modelmb_particles_copy(modelmb_out, modelmb_in)
  type(modelmb_particle_t), intent(in)  :: modelmb_in
  type(modelmb_particle_t), intent(out) :: modelmb_out

  call push_sub('states.modelmb_particles_copy')

  modelmb_out%ndim = modelmb_in%ndim
  modelmb_out%ntype_of_particle = modelmb_in%ntype_of_particle
  modelmb_out%max_particles_per_type = modelmb_in%max_particles_per_type
  modelmb_out%nparticle = modelmb_in%nparticle
  modelmb_out%ndensities_to_calculate = modelmb_in%ndensities_to_calculate

  call loct_pointer_copy(modelmb_out%labels_particles,modelmb_in%labels_particles)
  call loct_pointer_copy(modelmb_out%particletype,modelmb_in%particletype)
  call loct_pointer_copy(modelmb_out%mass_particle,modelmb_in%mass_particle)
  call loct_pointer_copy(modelmb_out%charge_particle,modelmb_in%charge_particle)
  call loct_pointer_copy(modelmb_out%nparticles_per_type,modelmb_in%nparticles_per_type)
  call loct_pointer_copy(modelmb_out%exchange_symmetry,modelmb_in%exchange_symmetry)
  call loct_pointer_copy(modelmb_out%bosonfermion,modelmb_in%bosonfermion)

  call loct_pointer_copy(modelmb_out%labels_densities,modelmb_in%labels_densities)
  call loct_pointer_copy(modelmb_out%particle_kept_densities,modelmb_in%particle_kept_densities)

  call pop_sub()

end subroutine modelmb_particles_copy

!==============================================================
!  Copy masses for particles to the derivative object
!==============================================================
subroutine modelmb_copy_masses (modelmbparticles,masses)
  type(modelmb_particle_t), intent(in)    :: modelmbparticles
  FLOAT,                    intent(inout) :: masses(MAX_DIM)

  integer :: dimcounter,ipart

  call push_sub('states.modelmb_copy_masses')

  ! copy masses to gr%der%masses
  dimcounter = 0
  do ipart = 1,modelmbparticles%nparticle
    if (abs(modelmbparticles%mass_particle(ipart)-1.0d0) > 1.e-10) then
      masses(dimcounter+1:dimcounter+modelmbparticles%ndim) = modelmbparticles%mass_particle(ipart)
    end if
    dimcounter = dimcounter+modelmbparticles%ndim
  end do

  call pop_sub()

end subroutine modelmb_copy_masses


end module modelmb_particles_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
