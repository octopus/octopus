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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

!>  general module for modelmb particles (e.g. 4 electrons in 1D equiv to
!!  1 in 4D). Also calculate different densities on request.
module modelmb_particles_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m

  implicit none

  private

  public :: modelmb_particles_init,      &
            modelmb_particles_end,       &
            modelmb_particles_copy,      &
            modelmb_copy_masses,         &
            modelmb_particle_t

  !>==============================================================
  !!  container type for input vars concerning modelmb particles
  !!==============================================================
  type modelmb_particle_t
    private
    integer, public :: ndim              !< dimensionality of modelmb space each
                                         !!  particle lives in
    
    integer, public :: ntype_of_particle !< number of different types of particles
                                         !!  modelmb in MAX_DIM dimensional space
    integer :: max_particles_per_type    !< max number of different particle

    integer, public :: nparticle = 0     !< number of particles 

    integer :: ndensities_to_calculate

    !>   %block describe_particles_modelmb
    !!   label1(char) | particletype1(integer) | mass1 | charge1 | fermion or boson or anyon
    !!   label2(char) | particletype2(integer) | mass2 | charge2 | fermion or boson or anyon
    !!   label3(char) | particletype3(integer) | mass3 | charge3 | fermion or boson or anyon
    !!   ...
    !!   nparticle_modelmb lines (fixed)
    !!   %
    character(80), allocatable :: labels_particles(:)

    integer, allocatable, public :: particletype(:)
    integer, allocatable, public :: nparticles_per_type(:)
    integer, allocatable, public :: particles_of_type(:,:)
    integer, allocatable, public :: bosonfermion(:)
    
    integer, allocatable :: exchange_symmetry(:,:,:) !< (max_particles_per_type**2, ntype_of_particle)
    
    FLOAT, allocatable :: mass_particle(:)
    
    FLOAT, allocatable, public :: charge_particle(:)

    !>   %block densitiestocalc
    !!   label1 |  particletokeep1(integer in [1:nparticle])
    !!   label2 |  particletokeep2(integer)
    !!   label3 |  particletokeep3(integer)
    !!   ...
    !!   however many lines wanted (up to nparticle_modelmb)
    !!   %
    character(80), allocatable :: labels_densities(:)
    
    integer, allocatable :: particle_kept_densities(:)

  end type modelmb_particle_t

contains
  
  !>==============================================================
  !!  initialization function for modelmb particles information
  !!==============================================================
  subroutine modelmb_particles_init(this, namespace, space)
    type(modelmb_particle_t), intent(inout) :: this
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    
    integer :: ipart, ncols, nline, itmp, jtmp, npar, ntype
    type(block_t) :: blk
    
    PUSH_SUB(modelmb_particles_init)
    
    ! read in scalar dimensions
    
    !%Variable NParticleModelmb
    !%Type integer
    !%Section States::ModelMB
    !%Default  0
    !%Description
    !% Number of particles in modelmb space. 
    !% Full Ndim = <tt>NDimModelmb</tt>*<tt>NParticleModelmb</tt>
    !%End
    call parse_variable(namespace, 'NParticleModelmb', 0, this%nparticle)

    if (this%nparticle == 0) then
      POP_SUB(modelmb_particles_init)
      return
    end if

    call messages_print_var_value(stdout, "NParticleModelmb", this%nparticle)
    
    !%Variable NDimModelmb
    !%Type integer
    !%Section States::ModelMB
    !%Default 1
    !%Description
    !% Number of dimensions for modelmb space.
    !% Full Ndim = <tt>NDimModelmb</tt>*<tt>NParticleModelmb</tt>
    !%
    !%End
    call parse_variable(namespace, 'NDimModelmb', 1, this%ndim)
    call messages_print_var_value(stdout, "NDimModelmb", this%ndim)
    
    !%Variable NTypeParticleModelmb
    !%Type integer
    !%Section States::ModelMB
    !%Default 1
    !%Description
    !% Number of different types of particles in modelmb space.
    !%End
    call parse_variable(namespace, 'NTypeParticleModelmb', 1, this%ntype_of_particle)
    call messages_print_var_value(stdout, "NTypeParticleModelmb", this%ntype_of_particle)
    if (this%ntype_of_particle > this%nparticle) then
      write (message(1), '(2a,2I6)') ' Number of types of modelmb particles should be <= Number of modelmb particles ', &
        this%ntype_of_particle, this%nparticle
      call messages_fatal(1, namespace=namespace)
    end if
    
    if (this%ndim*this%nparticle /= space%dim) then
      message(1) = ' Number of modelmb particles * dimension of modelmb space must be = Ndim'
      call messages_fatal(1, namespace=namespace)
    end if

    ! read in blocks
    
    !%Variable DescribeParticlesModelmb
    !%Type block
    !%Section States::ModelMB
    !%Description
    !% Characterization of different modelmb particles in space%dim dimensional space.
    !%
    !% <tt>%DescribeParticlesModelmb
    !% <br>&nbsp;&nbsp; "proton"   | 1 | 1800. | 1. | fermion
    !% <br>&nbsp;&nbsp; "proton"   | 1 | 1800. | 1. | fermion
    !% <br>&nbsp;&nbsp; "electron" | 2 | 1.    | 1. | fermion
    !% <br>%</tt>
    !%
    !% would tell <tt>Octopus</tt> that there are presently 3 particles, called proton, proton,
    !% and electron, with types 1, 1, and 2, and corresponding masses and charges.
    !% All particles should be fermions, and this can be later enforced on the spatial
    !% part of the wavefunctions.
    !% The label and charge are presently only for informational purposes and
    !% are not checked or used in <tt>Octopus</tt>. The interaction has to take the
    !% actual charge into account.
    !%
    !%Option fermion 1
    !%  Particle is a fermion.
    !%Option boson 2
    !%  Particle is a boson.
    !%Option anyon 3
    !%  Particle is neither fermion nor boson.
    !%End
    
    ! allocate stuff
    npar = this%nparticle
    ntype = this%ntype_of_particle
    SAFE_ALLOCATE (this%labels_particles(1:npar))
    SAFE_ALLOCATE (this%particletype(1:npar))
    SAFE_ALLOCATE (this%mass_particle(1:npar))
    SAFE_ALLOCATE (this%charge_particle(1:npar))
    SAFE_ALLOCATE (this%bosonfermion(1:npar))
    SAFE_ALLOCATE (this%nparticles_per_type(1:ntype))
    SAFE_ALLOCATE (this%particles_of_type(1:npar, 1:ntype))
    
    ! default all particles are electrons
    this%labels_particles = 'electron'
    this%particletype = 1
    this%mass_particle = M_ONE
    this%charge_particle = M_ONE
    this%bosonfermion = 1 ! set to fermion
    
    
    if(parse_block(namespace, 'DescribeParticlesModelmb', blk) == 0) then
      
      call messages_experimental("Model many-body")
      
      ncols = parse_block_cols(blk, 0)
      if(ncols /= 5 ) then
        call messages_input_error(namespace, "DescribeParticlesModelmb")
      end if
      nline = parse_block_n(blk)
      if (nline /= this%nparticle) then
        call messages_input_error(namespace, "DescribeParticlesModelmb")
      end if
      
      do ipart = 1, this%nparticle
        call parse_block_string(blk, ipart - 1, 0, this%labels_particles(ipart))
        call parse_block_integer(blk, ipart - 1, 1, this%particletype(ipart))
        call parse_block_float(blk, ipart - 1, 2, this%mass_particle(ipart))
        call parse_block_float(blk, ipart - 1, 3, this%charge_particle(ipart))
        call parse_block_integer(blk, ipart - 1, 4, this%bosonfermion(ipart))
        
        write (message(1),'(a,a)') 'labels_particles = ', this%labels_particles(ipart)
        write (message(2),'(a,i6)') 'particletype = ', this%particletype(ipart)
        write (message(3),'(a,E20.10)') 'mass_particle = ', this%mass_particle(ipart)
        write (message(4),'(a,E20.10)') 'charge_particle = ', this%charge_particle(ipart)
        write (message(5),'(a,i6)') 'bosonfermion = ', this%bosonfermion(ipart)
        call messages_info(5)
      end do
      call parse_block_end(blk)
      
    end if
    
    this%nparticles_per_type = 0
    this%particles_of_type = 0
    do ipart = 1, this%nparticle
      this%nparticles_per_type(this%particletype(ipart)) = & 
        this%nparticles_per_type(this%particletype(ipart)) + 1
      this%particles_of_type(this%nparticles_per_type(this%particletype(ipart)), &
        this%particletype(ipart)) = ipart
    end do
    
    this%max_particles_per_type = maxval(this%nparticles_per_type)
    itmp = this%max_particles_per_type
    jtmp = this%ntype_of_particle
    SAFE_ALLOCATE (this%exchange_symmetry(1:itmp, 1:itmp, 1:jtmp))
    this%exchange_symmetry = 0
    
    POP_SUB(modelmb_particles_init)
    
  end subroutine modelmb_particles_init


  subroutine modelmb_particles_end (this)
    type(modelmb_particle_t),intent(inout) :: this
    
    PUSH_SUB(modelmb_particles_end)
    
    SAFE_DEALLOCATE_A(this%labels_particles)
    SAFE_DEALLOCATE_A(this%particletype)
    SAFE_DEALLOCATE_A(this%mass_particle)
    SAFE_DEALLOCATE_A(this%charge_particle)
    SAFE_DEALLOCATE_A(this%nparticles_per_type)
    SAFE_DEALLOCATE_A(this%particles_of_type)
    SAFE_DEALLOCATE_A(this%exchange_symmetry)
    SAFE_DEALLOCATE_A(this%bosonfermion)
    
    SAFE_DEALLOCATE_A(this%labels_densities)
    SAFE_DEALLOCATE_A(this%particle_kept_densities)
    
    POP_SUB(modelmb_particles_end)
  end subroutine modelmb_particles_end
  
  subroutine modelmb_particles_copy(modelmb_out, modelmb_in)
    type(modelmb_particle_t), intent(in)    :: modelmb_in
    type(modelmb_particle_t), intent(inout) :: modelmb_out
    
    PUSH_SUB(modelmb_particles_copy)
    
    call modelmb_particles_end(modelmb_out)

    modelmb_out%ndim = modelmb_in%ndim
    modelmb_out%ntype_of_particle = modelmb_in%ntype_of_particle
    modelmb_out%max_particles_per_type = modelmb_in%max_particles_per_type
    modelmb_out%nparticle = modelmb_in%nparticle
    modelmb_out%ndensities_to_calculate = modelmb_in%ndensities_to_calculate
    
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%labels_particles,modelmb_in%labels_particles)
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%particletype,modelmb_in%particletype)
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%mass_particle,modelmb_in%mass_particle)
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%charge_particle,modelmb_in%charge_particle)
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%nparticles_per_type,modelmb_in%nparticles_per_type)
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%particles_of_type,modelmb_in%particles_of_type)
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%exchange_symmetry,modelmb_in%exchange_symmetry)
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%bosonfermion,modelmb_in%bosonfermion)
    
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%labels_densities,modelmb_in%labels_densities)
    SAFE_ALLOCATE_SOURCE_A(modelmb_out%particle_kept_densities,modelmb_in%particle_kept_densities)
    
    POP_SUB(modelmb_particles_copy)
  end subroutine modelmb_particles_copy

  !>==============================================================
  !!  Copy masses for particles to the derivative object
  !!==============================================================
  subroutine modelmb_copy_masses (this,masses)
    type(modelmb_particle_t), intent(in)    :: this
    FLOAT,                    intent(inout) :: masses(MAX_DIM)
    
    integer :: dimcounter,ipart
    
    PUSH_SUB(modelmb_copy_masses)
    
    ! copy masses to gr%der%masses
    dimcounter = 0
    do ipart = 1,this%nparticle
      if (abs(this%mass_particle(ipart)-1.0d0) > 1.e-10) then
        masses(dimcounter+1:dimcounter+this%ndim) = this%mass_particle(ipart)
      end if
      dimcounter = dimcounter+this%ndim
    end do
    
    POP_SUB(modelmb_copy_masses)
    
  end subroutine modelmb_copy_masses
  
end module modelmb_particles_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
