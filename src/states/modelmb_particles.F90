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
!  general module for modelMB particles (eg 4 electrons in 1D equiv to
!  1 in 4D). Also calculate different densities on request.
!
#include "global.h"

module modelMB_particles_m

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

  public :: modelMB_particles_nullify,&
            modelMB_particles_init,&
            modelMB_particles_end,&
            modelMB_particles_copy,&
            modelmb_copy_masses,&
            modelMB_particle_t

  public :: modelmb_1part_init, &
            modelmb_1part_nullify, &
            modelmb_1part_end, &
            modelmb_1part_t

!==============================================================
!  container type for input vars concerning modelMB particles
!==============================================================
type modelMB_particle_t
   integer :: ndim_modelMB              ! dimensionality of modelMB space each
                                         !  particle lives in

   integer :: ntype_of_particle_modelMB ! number of different types of particles
                                        !  modelMB in MAX_DIM dimensional space
   integer :: max_particles_per_type    ! max number of different particle

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
   integer, pointer :: nparticles_per_type(:)

   integer, pointer :: exchange_symmetry(:,:,:) ! (max_particles_per_type**2, ntype_of_particle_modelMB)
   character(80), pointer :: bosonfermion(:)

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

!
!  container type for the position and dimensions for 1 particle (out of
!  MAX_DIM/dims)
!
type modelmb_1part_t
  integer :: ndim1part
  integer :: npt_1part
  FLOAT :: vol_elem_1part
  FLOAT, pointer :: origin(:)
  integer, pointer :: enlarge_1part(:)
  integer, pointer :: nr_1part(:,:)
  FLOAT, pointer :: h_1part(:)
  type(hypercube_t) :: hypercube_1part
end type modelmb_1part_t

contains


subroutine modelMB_particles_nullify(this)
  type(modelMB_particle_t), intent(inout) :: this
  nullify(this%labels_particles_modelMB)
  nullify(this%particletype_modelMB)
  nullify(this%nparticles_per_type)
  nullify(this%exchange_symmetry)
  nullify(this%bosonfermion)
  nullify(this%mass_particle_modelMB)
  nullify(this%charge_particle_modelMB)
  nullify(this%labels_densities)
  nullify(this%particle_kept_densities)
end subroutine modelMB_particles_nullify


!==============================================================
!  initialization function for modelMB particles information
!==============================================================
subroutine modelMB_particles_init (modelMBparticles,gr)

  implicit none

!args
  type(modelMB_particle_t),intent(inout) :: modelMBparticles
  type(grid_t), intent(in) :: gr

!local vars
  integer :: ipart, ncols, nline, itmp, jtmp
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
  call loct_parse_int(datasets_check('NDimModelMB'), gr%sb%dim, modelMBparticles%ndim_modelMB)
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
  call loct_parse_int(datasets_check('NParticleModelMB'), 1, modelMBparticles%nparticle_modelMB)
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
     write (message(1), '(a,2I6)') ' Number of types of modelMB particles should be <= Number of modelMB particles ', &
             modelMBparticles%ntype_of_particle_modelMB, modelMBparticles%nparticle_modelMB
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
  !% Characterization of different modelMB particles in gr%mesh%sb%dim dimensional space
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
  SAFE_ALLOCATE (modelMBparticles%labels_particles_modelMB(1:modelMBparticles%nparticle_modelMB))
  SAFE_ALLOCATE (modelMBparticles%particletype_modelMB(1:modelMBparticles%nparticle_modelMB))
  SAFE_ALLOCATE (modelMBparticles%mass_particle_modelMB(1:modelMBparticles%nparticle_modelMB))
  SAFE_ALLOCATE (modelMBparticles%charge_particle_modelMB(1:modelMBparticles%nparticle_modelMB))
  SAFE_ALLOCATE (modelMBparticles%nparticles_per_type(1:modelMBparticles%ntype_of_particle_modelMB))

! default all particles are electrons
  modelMBparticles%labels_particles_modelMB = 'electron'
  modelMBparticles%particletype_modelMB = 1
  modelMBparticles%mass_particle_modelMB = 1.0d0
  modelMBparticles%charge_particle_modelMB = 1.0d0
    

  if(loct_parse_block(datasets_check('DescribeParticlesModelMB'), blk)==0) then
      ncols = loct_parse_block_cols(blk, 0)
      if(ncols /= 4 ) then
        call input_error("DescribeParticlesModelMB")
      end if
      nline = loct_parse_block_n(blk)
      if (nline /= modelMBparticles%nparticle_modelMB) then
        call input_error("DescribeParticlesModelMB")
      end if

      do ipart = 1,modelMBparticles%nparticle_modelMB
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
    
  modelMBparticles%nparticles_per_type = 0
  do ipart = 1, modelmbparticles%nparticle_modelmb
    modelMBparticles%nparticles_per_type(modelMBparticles%particletype_modelMB(ipart)) = & 
      & modelMBparticles%nparticles_per_type(modelMBparticles%particletype_modelMB(ipart)) + 1
  enddo

  modelMBparticles%max_particles_per_type = maxval(modelMBparticles%nparticles_per_type)
  itmp = modelMBparticles%max_particles_per_type
  jtmp = modelMBparticles%ntype_of_particle_modelMB
  SAFE_ALLOCATE (modelMBparticles%exchange_symmetry(1:itmp, 1:itmp, 1:jtmp))
  modelMBparticles%exchange_symmetry = 0
  SAFE_ALLOCATE (modelMBparticles%bosonfermion(1:jtmp))
  modelMBparticles%bosonfermion = 'uncalculated'

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
  modelMBparticles%ndensities_to_calculate = 0

  if(loct_parse_block(datasets_check('DensitiestoCalc'), blk)==0) then
    ncols = loct_parse_block_cols(blk, 0)
    if(ncols /= 2 ) then
      call input_error("DensitiestoCalc")
    end if
    modelMBparticles%ndensities_to_calculate = loct_parse_block_n(blk)
    if (modelMBparticles%ndensities_to_calculate < 0 .or. &
         modelMBparticles%ndensities_to_calculate > modelMBparticles%nparticle_modelMB) then
      call input_error("DensitiestoCalc")
    end if
    SAFE_ALLOCATE (modelMBparticles%labels_densities(1:modelMBparticles%ndensities_to_calculate))
    SAFE_ALLOCATE (modelMBparticles%particle_kept_densities(1:modelMBparticles%ndensities_to_calculate))
    do ipart = 1,modelMBparticles%ndensities_to_calculate
      call loct_parse_block_string(blk, ipart-1, 0, modelMBparticles%labels_densities(ipart))
      call loct_parse_block_int(blk, ipart-1, 1, modelMBparticles%particle_kept_densities(ipart))
      
      write (message(1),'(a,a)') 'labels_densities = ', modelMBparticles%labels_densities(ipart)
      write (message(2),'(a,i6)') 'particle_kept_densities = ', modelMBparticles%particle_kept_densities(ipart)
      call write_info(2)
    end do
    call loct_parse_block_end(blk)
  else
    nullify(modelmbparticles%labels_densities)
    nullify(modelmbparticles%particle_kept_densities)
  end if

  call pop_sub()

end subroutine modelMB_particles_init


subroutine modelMB_particles_end (modelMBparticles)

  implicit none

!args
  type(modelMB_particle_t),intent(inout) :: modelMBparticles

! source

  call push_sub('states.modelMB_particles_end')

  SAFE_DEALLOCATE_P(modelMBparticles%labels_particles_modelMB)
  SAFE_DEALLOCATE_P(modelMBparticles%particletype_modelMB)
  SAFE_DEALLOCATE_P(modelMBparticles%mass_particle_modelMB)
  SAFE_DEALLOCATE_P(modelMBparticles%charge_particle_modelMB)
  SAFE_DEALLOCATE_P(modelMBparticles%nparticles_per_type)
  SAFE_DEALLOCATE_P(modelMBparticles%exchange_symmetry)
  SAFE_DEALLOCATE_P(modelMBparticles%bosonfermion)

  SAFE_DEALLOCATE_P(modelMBparticles%labels_densities)
  SAFE_DEALLOCATE_P(modelMBparticles%particle_kept_densities)

  call pop_sub()

end subroutine modelMB_particles_end

subroutine modelMB_particles_copy(modelmb_out, modelmb_in)

  implicit none

!args
  type(modelMB_particle_t),intent(in) :: modelmb_in
  type(modelMB_particle_t),intent(out) :: modelmb_out


  call push_sub('states.modelMB_particles_copy')

  modelmb_out%ndim_modelMB = modelmb_in%ndim_modelMB
  modelmb_out%ntype_of_particle_modelMB = modelmb_in%ntype_of_particle_modelMB
  modelmb_out%max_particles_per_type = modelmb_in%max_particles_per_type
  modelmb_out%nparticle_modelMB = modelmb_in%nparticle_modelMB
  modelmb_out%ndensities_to_calculate = modelmb_in%ndensities_to_calculate

  call loct_pointer_copy(modelmb_out%labels_particles_modelMB,modelmb_in%labels_particles_modelMB)
  call loct_pointer_copy(modelmb_out%particletype_modelMB,modelmb_in%particletype_modelMB)
  call loct_pointer_copy(modelmb_out%mass_particle_modelMB,modelmb_in%mass_particle_modelMB)
  call loct_pointer_copy(modelmb_out%charge_particle_modelMB,modelmb_in%charge_particle_modelMB)
  call loct_pointer_copy(modelmb_out%nparticles_per_type,modelmb_in%nparticles_per_type)
  call loct_pointer_copy(modelmb_out%exchange_symmetry,modelmb_in%exchange_symmetry)
  call loct_pointer_copy(modelmb_out%bosonfermion,modelmb_in%bosonfermion)

  call loct_pointer_copy(modelmb_out%labels_densities,modelmb_in%labels_densities)
  call loct_pointer_copy(modelmb_out%particle_kept_densities,modelmb_in%particle_kept_densities)

  call pop_sub()

end subroutine modelMB_particles_copy

!==============================================================
!  Copy masses for particles to the derivative object
!==============================================================
subroutine modelmb_copy_masses (modelMBparticles,masses)

  implicit none

!args
  type(modelMB_particle_t),intent(in) :: modelMBparticles
  FLOAT, intent(inout) :: masses(MAX_DIM)

!local
  integer :: dimcounter,ipart

  call push_sub('states.modelmb_copy_masses')

  ! copy masses to gr%der%masses
  dimcounter = 0
  do ipart = 1,modelMBparticles%nparticle_modelMB
    if (abs(modelMBparticles%mass_particle_modelMB(ipart)-1.0d0) > 1.e-10) then
      masses(dimcounter+1:dimcounter+modelMBparticles%ndim_modelMB) = modelMBparticles%mass_particle_modelMB(ipart)
    end if
    dimcounter = dimcounter+modelMBparticles%ndim_modelMB
  end do

  call pop_sub()

end subroutine modelmb_copy_masses


subroutine modelmb_1part_init(this, mesh, ikeeppart, ndim1part, box_offset)
  integer, intent(in) :: ikeeppart, ndim1part
  FLOAT, intent(in) :: box_offset(MAX_DIM)
  type(modelmb_1part_t), intent(out) :: this
  type(mesh_t), intent(in) :: mesh

  !local vars
  integer :: idir, irealdir

  call push_sub('states.modelmb_1part_init')
  
   this%ndim1part = ndim1part

!   get full size of arrays for 1 particle only in ndim_modelmb dimensions
   this%npt_1part = 1
   do idir = 1, ndim1part
     this%npt_1part = this%npt_1part*mesh%idx%ll((ikeeppart-1)*ndim1part+idir)
   end do

!   volume element for the chosen particle
  SAFE_ALLOCATE(this%h_1part(1:ndim1part))
  this%vol_elem_1part = 1.0d0
  do idir = 1,ndim1part
    irealdir = (ikeeppart-1)*ndim1part + idir
    this%vol_elem_1part = this%vol_elem_1part*mesh%h(irealdir)
    this%h_1part(idir) = mesh%h(irealdir)
  end do

!   store start and end positions for the relevant dimensions for this particle
  SAFE_ALLOCATE(this%nr_1part(1:2, 1:ndim1part))
  this%nr_1part(:,:) = mesh%idx%nr(:,(ikeeppart-1)*ndim1part+1:ikeeppart*ndim1part)

!   initialize a hypercube for just this particle
!   NB: hypercube_* presume that enlarge is the same for all dimensions!
  SAFE_ALLOCATE(this%enlarge_1part(1:ndim1part))
  this%enlarge_1part = mesh%idx%enlarge((ikeeppart-1)*ndim1part+1:ikeeppart*ndim1part)
  call hypercube_init(this%hypercube_1part, ndim1part, this%nr_1part, this%enlarge_1part(1))

  ! not always the real origin if the box is shifted, no?
  !  which happens to be my case...
  !  only important for printout, so it is ok
  SAFE_ALLOCATE(this%origin(1:ndim1part))
  do idir = 1,ndim1part
    irealdir = (ikeeppart-1)*ndim1part + idir
    !origin(idir) = (npoints(irealdir)/2)*gr%mesh%h(irealdir)
    this%origin(idir) = box_offset(irealdir)
  end do

  call pop_sub()
end subroutine modelmb_1part_init


subroutine modelmb_1part_nullify(this)
  type(modelmb_1part_t), intent(out) :: this
  call push_sub('states.modelmb_1part_nullify')
  nullify(this%origin)
  nullify(this%enlarge_1part)
  nullify(this%nr_1part)
  nullify(this%h_1part)
  call hypercube_nullify(this%hypercube_1part)
  call pop_sub()
end subroutine modelmb_1part_nullify


subroutine modelmb_1part_end(this)
  type(modelmb_1part_t), intent(inout) :: this
  call push_sub('states.modelmb_1part_end')
  SAFE_DEALLOCATE_P(this%origin)
  SAFE_DEALLOCATE_P(this%enlarge_1part)
  SAFE_DEALLOCATE_P(this%nr_1part)
  SAFE_DEALLOCATE_P(this%h_1part)
  call hypercube_end(this%hypercube_1part)
  call pop_sub()
end subroutine modelmb_1part_end


end module modelMB_particles_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
