!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 N. Tancogne-Dejean
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

module ion_interaction_oct_m
  use atom_oct_m
  use comm_oct_m
  use global_oct_m
  use distributed_oct_m
  use lattice_vectors_oct_m
  use loct_math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use periodic_copy_oct_m
  use profiling_oct_m
  use ps_oct_m
  use space_oct_m
  use species_oct_m

  implicit none

  private
  public ::                           &
    ion_interaction_t,                &
    ion_interaction_init,             &
    ion_interaction_end,              &
    ion_interaction_calculate,        &
    ion_interaction_init_parallelization, &
    ion_interaction_test

  type ion_interaction_t
    ! Components are public by default
    FLOAT                      :: alpha

    type(distributed_t) :: dist
  end type ion_interaction_t

  integer, parameter ::            &
    ION_COMPONENT_REAL     = 1,    &
    ION_COMPONENT_SELF     = 2,    &
    ION_COMPONENT_FOURIER  = 3,    &
    ION_NUM_COMPONENTS     = 3
  
contains

  subroutine ion_interaction_init(this, namespace, space, natoms)
    type(ion_interaction_t),      intent(out)   :: this
    type(namespace_t),            intent(in)    :: namespace
    type(space_t),                intent(in)    :: space
    integer,                      intent(in)    :: natoms

    PUSH_SUB(ion_interaction_init)

    !%Variable EwaldAlpha
    !%Type float
    !%Default 0.21
    !%Section Hamiltonian
    !%Description
    !% The value 'Alpha' that controls the splitting of the Coulomb
    !% interaction in the Ewald sum used to calculation the ion-ion
    !% interaction for periodic systems. This value affects the speed
    !% of the calculation, normally users do not need to modify it.
    !%End
    call parse_variable(namespace, 'EwaldAlpha', CNST(0.21), this%alpha)

    call distributed_nullify(this%dist, natoms)

    if(space%periodic_dim == 1) then
      call messages_write('For systems that  are periodic in 1D, interaction between', new_line = .true.)
      call messages_write('ions is assumed to be periodic in 3D. This affects the calculation', new_line = .true.)
      call messages_write('of total energy and forces.')
      call messages_warning(namespace=namespace)
    end if

    
    POP_SUB(ion_interaction_init)
  end subroutine ion_interaction_init

 ! ---------------------------------------------------------  

  subroutine ion_interaction_init_parallelization(this, natoms, mc)
    type(ion_interaction_t),      intent(inout) :: this
    integer,                      intent(in)    :: natoms
    type(multicomm_t),            intent(in)    :: mc

    PUSH_SUB(ion_interaction_init_parallelization)

    !As the code below is not parallelized with any of k-point, states nor domain
    !we can safely parallelize it over atoms
    call distributed_init(this%dist, natoms, mc%master_comm, "Ions")

    POP_SUB(ion_interaction_init_parallelization)
  end subroutine ion_interaction_init_parallelization

  ! ---------------------------------------------------------
  
  subroutine ion_interaction_end(this)
    type(ion_interaction_t), intent(inout) :: this

    PUSH_SUB(ion_interaction_end)

    this%alpha = -CNST(1.0)

    call distributed_end(this%dist)

    POP_SUB(ion_interaction_end)
  end subroutine ion_interaction_end
  
  ! ---------------------------------------------------------
  !> For details about this routine, see
  !! http://octopus-code.org/wiki/Developers:Ion-Ion_interaction
  subroutine ion_interaction_calculate(this, space, latt, atom, natoms, catom, ncatoms, lsize, energy, force, &
    energy_components, force_components)
    type(ion_interaction_t),  intent(inout) :: this
    type(space_t),            intent(in)    :: space
    type(lattice_vectors_t),  intent(in)    :: latt
    type(atom_t),             intent(in)    :: atom(:)
    integer,                  intent(in)    :: natoms
    type(atom_classical_t),   intent(in)    :: catom(:)
    integer,                  intent(in)    :: ncatoms
    FLOAT,                    intent(in)    :: lsize(:)
    FLOAT,                    intent(out)   :: energy
    FLOAT,                    intent(out)   :: force(:, :)
    FLOAT, optional,          intent(out)   :: energy_components(:)
    FLOAT, optional,          intent(out)   :: force_components(:, :, :)
    
    FLOAT, allocatable:: r(:), f(:)
    FLOAT :: rr, dd, zi, zj
    integer :: iatom, jatom, natom, iindex, jindex
    type(species_t), pointer :: spci, spcj
    type(profile_t), save :: ion_ion_prof

    PUSH_SUB(ion_interaction_calculate)
    call profiling_in(ion_ion_prof, "ION_ION_INTERACTION")

    SAFE_ALLOCATE(r(1:space%dim))
    SAFE_ALLOCATE(f(1:space%dim))

    if(present(energy_components)) then
      ASSERT(ubound(energy_components, dim = 1) == ION_NUM_COMPONENTS)
      energy_components = CNST(0.0)
    end if      

    if(present(force_components)) then
      ASSERT(all(ubound(force_components) == (/space%dim, natoms, ION_NUM_COMPONENTS/)))
      force_components = CNST(0.0)
    end if 

    energy = M_ZERO
    force(1:space%dim, 1:natoms) = M_ZERO

    if(space%is_periodic()) then

      ASSERT(ncatoms==0)

      spci => atom(1)%species
      ! This depends on the area, but we should check if it is fully consistent.        
      if( species_type(spci) == SPECIES_JELLIUM_SLAB ) then
        energy = energy + &
          M_PI*species_zval(spci)**2/(M_FOUR*lsize(1)*lsize(2))*(lsize(3) - species_jthick(spci)/M_THREE)
      else
        call ion_interaction_periodic(this, space, latt, atom, natoms, lsize, energy, force, energy_components, force_components)
      end if

    else
      
      natom = natoms + ncatoms

      ! only interaction inside the cell
      do iatom = this%dist%start, this%dist%end
        spci => atom(iatom)%species
        zi = species_zval(spci)

        select case(species_type(spci))
        case(SPECIES_JELLIUM)
          energy = energy + (M_THREE/M_FIVE)*zi**2/species_jradius(spci)
          ! The part depending on the simulation sphere is neglected
          
        case(SPECIES_JELLIUM_SLAB)
          energy = energy - M_PI*zi**2/(M_FOUR*lsize(1)*lsize(2))*species_jthick(spci)/M_THREE
          ! The part depending on the simulation box transverse dimension is neglected
        end select
        
        do jatom = iatom + 1, natoms
          
          spcj => atom(jatom)%species

          r(1:space%dim) = atom(iatom)%x(1:space%dim) - atom(jatom)%x(1:space%dim)

          rr = sqrt(sum(r**2))
          iindex = species_index(spci)
          jindex = species_index(spcj)
          
          zj = species_zval(spcj)
          !the force
          dd = zi*zj/rr
          f(1:space%dim) = (dd/rr**2)*r(1:space%dim)
          force(1:space%dim,iatom) = force(1:space%dim,iatom) + f(1:space%dim)
          force(1:space%dim,jatom) = force(1:space%dim,jatom) - f(1:space%dim)
          !energy
          energy=energy + dd
          
        end do !jatom
      end do !iatom

      do iatom = this%dist%start, this%dist%end

        do jatom = 1, ncatoms
          r(1:space%dim) = atom(iatom)%x(1:space%dim) - catom(jatom)%x(1:space%dim)
          rr = sqrt(sum(r**2))

          zi = species_zval(atom(iatom)%species)
          zj = catom(jatom)%charge
          !the force
          dd = zi*zj/rr
          force(1:space%dim,iatom) = force(1:space%dim,iatom) + (dd/rr*2)*r(1:space%dim)
          !energy
          energy = energy + dd
          
        end do !jatom
      end do !iatom

      call comm_allreduce(this%dist%mpi_grp, energy)
      call comm_allreduce(this%dist%mpi_grp, force)
      
    end if

    SAFE_DEALLOCATE_A(r)
    SAFE_DEALLOCATE_A(f)
    
    call profiling_out(ion_ion_prof)
    
    POP_SUB(ion_interaction_calculate)
  end subroutine ion_interaction_calculate

  ! ---------------------------------------------------------
  
  subroutine ion_interaction_periodic(this, space, latt, atom, natoms, lsize, energy, force, energy_components, force_components)
    type(ion_interaction_t),   intent(in)    :: this
    type(space_t),             intent(in)    :: space
    type(lattice_vectors_t),   intent(in)    :: latt
    type(atom_t),              intent(in)    :: atom(:)
    integer,                   intent(in)    :: natoms
    FLOAT,                     intent(in)    :: lsize(:)
    FLOAT,                     intent(out)   :: energy
    FLOAT,                     intent(out)   :: force(:, :) !< (space%dim, natoms)
    FLOAT, optional,           intent(out)   :: energy_components(:)
    FLOAT, optional,           intent(out)   :: force_components(:, :, :)

    FLOAT :: rr, xi(1:MAX_DIM), zi, zj, ereal, efourier, eself, erfc, rcut, epseudo
    integer :: iatom, jatom, icopy
    type(periodic_copy_t) :: pc
    FLOAT   :: charge
    type(profile_t), save :: prof_short, prof_long
    type(ps_t), pointer :: spec_ps

    PUSH_SUB(ion_interaction_periodic)

    ereal = M_ZERO

    force(1:space%dim, 1:natoms) = M_ZERO

    ! if the system is periodic we have to add the energy of the
    ! interaction with the copies
    
    rcut = CNST(6.0)/this%alpha

    call profiling_in(prof_short, "EWALD_SHORT")
    
    ! the short-range part is calculated directly
    do iatom = this%dist%start, this%dist%end
      if (.not. species_represents_real_atom(atom(iatom)%species)) cycle
      zi = species_zval(atom(iatom)%species)

      call periodic_copy_init(pc, space, latt, lsize, atom(iatom)%x, rcut)
      
      do icopy = 1, periodic_copy_num(pc)
        xi(1:space%dim) = periodic_copy_position(pc, space, latt, lsize, icopy)
        
        do jatom = 1,  natoms
          zj = species_zval(atom(jatom)%species)
          rr = sqrt( sum( (xi(1:space%dim) - atom(jatom)%x(1:space%dim))**2 ) )
          
          if(rr < CNST(1e-5)) cycle
          if(rr > rcut) cycle
          
          erfc = M_ONE - loct_erf(this%alpha*rr)

          ! energy
          ereal = ereal + M_HALF*zj*zi*erfc/rr
          
          ! force
          force(1:space%dim, jatom) = force(1:space%dim, jatom) - &
            zj*zi*(xi(1:space%dim) - atom(jatom)%x(1:space%dim))*&
            (erfc/rr + M_TWO*this%alpha/sqrt(M_PI)*exp(-(this%alpha*rr)**2))/rr**2
        end do

      end do
      
      call periodic_copy_end(pc)
    end do

    call comm_allreduce(this%dist%mpi_grp, ereal)
    call comm_allreduce(this%dist%mpi_grp, force)

    if(present(force_components)) then
      force_components(1:space%dim, 1:natoms, ION_COMPONENT_REAL) = force(1:space%dim, 1:natoms)
    end if

    call profiling_out(prof_short)

    call profiling_in(prof_long, "EWALD_LONG")
    
    ! self-interaction
    eself = M_ZERO
    charge = M_ZERO
    do iatom = this%dist%start, this%dist%end
      zi = species_zval(atom(iatom)%species)
      charge = charge + zi
      eself = eself - this%alpha/sqrt(M_PI)*zi**2
    end do
    call comm_allreduce(this%dist%mpi_grp, eself)
    call comm_allreduce(this%dist%mpi_grp, charge)

    ! Long range part of Ewald sum
    select case (space%periodic_dim)
    case(1)
      ! Not implemented.
      efourier = M_ZERO
      force = M_ZERO
      ! Do not confuse the user and set to zero all the other components
      ereal = M_ZERO
      eself = M_ZERO
      if (present(force_components)) then
        force_components(1:space%dim, 1:natoms, ION_COMPONENT_REAL) = M_ZERO
      end if
    case(2)
      call Ewald_long_2D(this, space, latt, atom, natoms, efourier, force)
    case(3)
      call Ewald_long_3D(this, space, latt, atom, natoms, efourier, force, charge)
    end select


    if(present(energy_components)) then
      energy_components(ION_COMPONENT_REAL) = ereal
      energy_components(ION_COMPONENT_SELF) = eself
      energy_components(ION_COMPONENT_FOURIER) = efourier      
    end if

    if(present(force_components)) then
      force_components(1:space%dim, 1:natoms, ION_COMPONENT_FOURIER) = &
        force(1:space%dim, 1:natoms) - force_components(1:space%dim, 1:natoms, ION_COMPONENT_REAL)
    end if

    energy = ereal + efourier + eself

    ! Warning: The energy contribution of the long range part of the pseudo is
    ! not correctly accounted for in systems periodic in 1D or 2D.
    epseudo = M_ZERO
    if (space%periodic_dim == 3) then
       ! Previously unaccounted G = 0 term from pseudopotentials. 
       ! See J. Ihm, A. Zunger, M.L. Cohen, J. Phys. C 12, 4409 (1979)

       do iatom = this%dist%start, this%dist%end
          if(species_is_ps(atom(iatom)%species)) then
             zi = species_zval(atom(iatom)%species)
             spec_ps => species_ps(atom(iatom)%species)
             epseudo = epseudo + M_PI*zi*&
                  (spec_ps%sigma_erf*sqrt(M_TWO))**2/latt%rcell_volume*charge
          end if
       end do
       call comm_allreduce(this%dist%mpi_grp, epseudo)

       energy = energy + epseudo
    end if
    
    call profiling_out(prof_long)
    
    POP_SUB(ion_interaction_periodic)
  end subroutine ion_interaction_periodic

  ! ---------------------------------------------------------

  subroutine Ewald_long_3D(this, space, latt, atom, natoms, efourier, force, charge)
    type(ion_interaction_t),   intent(in)   :: this
    type(space_t),             intent(in)   :: space
    type(lattice_vectors_t),   intent(in)   :: latt
    type(atom_t),              intent(in)   :: atom(:)
    integer,                   intent(in)   :: natoms
    FLOAT,                     intent(inout) :: efourier
    FLOAT,                     intent(inout) :: force(:, :) !< (space%dim, natoms)
    FLOAT,                     intent(in)   :: charge

    FLOAT :: rcut
    integer :: iatom
    integer :: ix, iy, iz, isph, ss, idim
    FLOAT   :: gg(1:MAX_DIM), gg2, gx
    FLOAT   :: factor
    CMPLX   :: sumatoms, tmp(1:MAX_DIM), aa

    CMPLX, allocatable :: phase(:)

    PUSH_SUB(Ewald_long_3d)

    ASSERT(space%dim == 3)
    ASSERT(space%periodic_dim == 3)

    ! And the long-range part, using an Ewald sum
    SAFE_ALLOCATE(phase(1:natoms))

    ! get a converged value for the cutoff in g
    rcut = huge(rcut)
    do idim = 1, space%dim
      rcut = min(rcut, sum(latt%klattice(1:space%dim, idim)**2))
    end do

    rcut = sqrt(rcut)
    
    isph = ceiling(CNST(9.5)*this%alpha/rcut)

    ! First the G = 0 term (charge was calculated previously)
    efourier = -M_PI*charge**2/(M_TWO*this%alpha**2*latt%rcell_volume)

    do ix = -isph, isph
      do iy = -isph, isph
        do iz = -isph, isph
          
          ss = ix**2 + iy**2 + iz**2
          
          if(ss == 0 .or. ss > isph**2) cycle

          gg(1:space%dim) = ix*latt%klattice(1:space%dim, 1) &
                              + iy*latt%klattice(1:space%dim, 2) &
                              + iz*latt%klattice(1:space%dim, 3)
          gg2 = sum(gg(1:space%dim)**2)

          ! g=0 must be removed from the sum
          if(gg2 < M_EPSILON) cycle
          
          gx = -CNST(0.25)*gg2/this%alpha**2

          if(gx < CNST(-36.0)) cycle

          factor = M_TWO*M_PI/latt%rcell_volume*exp(gx)/gg2

          if(factor < epsilon(factor)) cycle

          sumatoms = M_Z0
          !$omp parallel do private(iatom, gx, aa) reduction(+:sumatoms)
          do iatom = 1, natoms
            gx = sum(gg(1:space%dim)*atom(iatom)%x(1:space%dim))
            aa = species_zval(atom(iatom)%species)*TOCMPLX(cos(gx), sin(gx))
            phase(iatom) = aa 
            sumatoms = sumatoms + aa
          end do
          
          efourier = efourier + TOFLOAT(factor*sumatoms*conjg(sumatoms))
          
          do iatom = 1, natoms
            tmp(1:space%dim) = M_ZI*gg(1:space%dim)*phase(iatom)
            force(1:space%dim, iatom) = force(1:space%dim, iatom) &
              - factor*TOFLOAT(conjg(tmp(1:space%dim))*sumatoms + tmp(1:space%dim)*conjg(sumatoms))
          end do
          
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(phase)

    POP_SUB(Ewald_long_3d)


  end subroutine Ewald_long_3D

  ! ---------------------------------------------------------
  !> In-Chul Yeh and Max L. Berkowitz, J. Chem. Phys. 111, 3155 (1999).
  subroutine Ewald_long_2D(this, space, latt, atom, natoms, efourier, force)
    type(ion_interaction_t),   intent(in)    :: this
    type(space_t),             intent(in)    :: space
    type(lattice_vectors_t),   intent(in)    :: latt
    type(atom_t),              intent(in)    :: atom(:)
    integer,                   intent(in)    :: natoms
    FLOAT,                     intent(inout) :: efourier
    FLOAT,                     intent(inout) :: force(:, :) !< (space%dim, natoms)

    FLOAT :: rcut
    integer :: iatom, jatom
    integer :: ix, iy, ix_max, iy_max, ss
    FLOAT   :: gg(1:MAX_DIM), gg2, gx, gg_abs
    FLOAT   :: factor,factor1,factor2, coeff
    FLOAT   :: dz_max, dz_ij, erfc1, erfc2, tmp_erf
    FLOAT, allocatable :: force_tmp(:,:)

    PUSH_SUB(Ewald_long_2d)

    ASSERT(space%periodic_dim == 2)
    ASSERT(space%dim == 2 .or. space%dim == 3)

    ! And the long-range part, using an Ewald sum

    ! Searching maximum distance
    if (space%dim == 3) then
      dz_max = M_ZERO
      do iatom = 1, natoms
        do jatom = iatom + 1, natoms
          dz_max = max(dz_max, abs(atom(iatom)%x(3) - atom(jatom)%x(3)))
        end do
      end do
    else
      ! For a 2D system, all atoms are on the plane, so the distance is zero
      dz_max = M_ZERO
    end if

    !get a converged value for the cutoff in g
    rcut = M_TWO*this%alpha*CNST(4.6) + M_TWO*this%alpha**2*dz_max
    if(rcut > M_ZERO) then
      do
        if( rcut * dz_max >= 718 ) exit  !Maximum double precision numbber
        erfc1 = M_ONE - loct_erf(this%alpha*dz_max + M_HALF*rcut/this%alpha)
        if( erfc1 * exp(rcut*dz_max) < CNST(1e-10) ) exit
        rcut = rcut * CNST(1.414)
      end do
    end if

    ix_max = ceiling(rcut/sqrt(sum(latt%klattice(1:space%dim, 1)**2)))
    iy_max = ceiling(rcut/sqrt(sum(latt%klattice(1:space%dim, 2)**2)))

    SAFE_ALLOCATE(force_tmp(1:space%dim, 1:natoms))
    force_tmp = M_ZERO

    ! First the G = 0 term (charge was calculated previously)
    efourier = M_ZERO
    factor = M_PI/latt%rcell_volume
    do iatom = this%dist%start, this%dist%end
      do jatom = 1, natoms
        ! efourier
        if (space%dim == 3) then
          dz_ij = atom(iatom)%x(3) - atom(jatom)%x(3)
        else
          dz_ij = M_ZERO
        end if

        tmp_erf = loct_erf(this%alpha*dz_ij)
        factor1 = dz_ij*tmp_erf
        factor2 = exp(-(this%alpha*dz_ij)**2)/(this%alpha*sqrt(M_PI))

        efourier = efourier - factor &
          * species_zval(atom(iatom)%species)*species_zval(atom(jatom)%species) * (factor1 + factor2)

        ! force
        if (iatom == jatom)cycle
        if (abs(tmp_erf) < M_EPSILON) cycle

        if (space%dim == 3) then
          force_tmp(3, iatom) = force_tmp(3, iatom) - (- M_TWO*factor) &
            * species_zval(atom(iatom)%species)*species_zval(atom(jatom)%species) * tmp_erf
        end if

      end do
    end do

    do ix = -ix_max, ix_max
      do iy = -iy_max, iy_max

        ss = ix**2 + iy**2
        if(ss == 0) cycle
          
        gg(1:space%dim) = ix*latt%klattice(1:space%dim, 1) + iy*latt%klattice(1:space%dim, 2)
        gg2 = sum(gg(1:space%dim)**2)

        ! g=0 must be removed from the sum
        if(gg2 < M_EPSILON) cycle
        gg_abs = sqrt(gg2)
        factor = M_HALF*M_PI/(latt%rcell_volume*gg_abs)
          
        do iatom = this%dist%start, this%dist%end
          do jatom = iatom, natoms
            ! efourier
            gx = gg(1)*(atom(iatom)%x(1)-atom(jatom)%x(1)) &
              + gg(2)*(atom(iatom)%x(2)-atom(jatom)%x(2))
            if (space%dim == 3) then
              dz_ij = atom(iatom)%x(3)-atom(jatom)%x(3)
            else
              dz_ij = M_ZERO
            end if

            erfc1 = M_ONE - loct_erf(this%alpha*dz_ij + M_HALF*gg_abs/this%alpha)
            if(abs(erfc1) > M_EPSILON) then
              factor1 = exp(gg_abs*dz_ij)*erfc1
            else
              factor1 = M_ZERO
            end if
            erfc2 = M_ONE - loct_erf(-this%alpha*dz_ij + M_HALF*gg_abs/this%alpha)
            if(abs(erfc2) > M_EPSILON) then
              factor2 = exp(-gg_abs*dz_ij)*erfc2
            else
              factor2 = M_ZERO
            end if

            if(iatom == jatom) then
              coeff = M_ONE
            else
              coeff = M_TWO
            end if

            efourier = efourier &
              + factor * coeff &
              * species_zval(atom(iatom)%species)*species_zval(atom(jatom)%species) &
              * cos(gx)* ( factor1 + factor2)
              
            ! force
            if (iatom == jatom) cycle

            force_tmp(1:2, iatom) = force_tmp(1:2, iatom) &
              - (CNST(-1.0)* M_TWO*factor )* gg(1:2) &
              * species_zval(atom(iatom)%species)*species_zval(atom(jatom)%species) &
              *sin(gx)*(factor1 + factor2)

            force_tmp(1:2, jatom) = force_tmp(1:2, jatom) &
              + (CNST(-1.0)* M_TWO*factor )* gg(1:2) &
              * species_zval(atom(iatom)%species)*species_zval(atom(jatom)%species) &
              *sin(gx)*(factor1 + factor2)

            factor1 = gg_abs*erfc1 &
              - M_TWO*this%alpha/sqrt(M_PI)*exp(-(this%alpha*dz_ij + M_HALF*gg_abs/this%alpha)**2)
            if(abs(factor1) > M_EPSILON) then
              factor1 = factor1*exp(gg_abs*dz_ij)
            else
              factor1 = M_ZERO
            end if

            factor2 = gg_abs*erfc2 &
              - M_TWO*this%alpha/sqrt(M_PI)*exp(-(-this%alpha*dz_ij + M_HALF*gg_abs/this%alpha)**2)
            if(abs(factor2) > M_EPSILON) then
              factor2 = factor2*exp(-gg_abs*dz_ij)
            else
              factor2 = M_ZERO
            end if
             
            if (space%dim == 3) then
              force_tmp(3, iatom) = force_tmp(3, iatom) &
                - M_TWO*factor &
                * species_zval(atom(iatom)%species)*species_zval(atom(jatom)%species) &
                * cos(gx)* ( factor1 - factor2)
              force_tmp(3, jatom) = force_tmp(3, jatom) &
                + M_TWO*factor &
                * species_zval(atom(iatom)%species)*species_zval(atom(jatom)%species) &
                * cos(gx)* ( factor1 - factor2)
            end if

          end do
        end do


      end do
    end do

    call comm_allreduce(this%dist%mpi_grp, efourier)
    call comm_allreduce(this%dist%mpi_grp, force_tmp)

    force = force + force_tmp

    SAFE_DEALLOCATE_A(force_tmp)

    POP_SUB(Ewald_long_2d)
  end subroutine Ewald_long_2D

  ! ---------------------------------------------------------
  
  subroutine ion_interaction_test(space, latt, atom, natoms, catom, ncatoms, lsize, &
                namespace, mc)
    type(space_t),            intent(in)    :: space
    type(lattice_vectors_t),  intent(in)    :: latt
    type(atom_t),             intent(in)    :: atom(:)
    integer,                  intent(in)    :: natoms
    type(atom_classical_t),   intent(in)    :: catom(:)
    integer,                  intent(in)    :: ncatoms
    FLOAT,                    intent(in)    :: lsize(:)
    type(namespace_t),        intent(in)    :: namespace
    type(multicomm_t),        intent(in)    :: mc

    type(ion_interaction_t) :: ion_interaction
    FLOAT :: energy
    FLOAT, allocatable :: force(:, :), force_components(:, :, :)
    FLOAT :: energy_components(1:ION_NUM_COMPONENTS)
    integer :: iatom, idir
    
    PUSH_SUB(ion_interaction_test)

    call ion_interaction_init(ion_interaction, namespace, space, natoms)
    call ion_interaction_init_parallelization(ion_interaction, natoms, mc)

    SAFE_ALLOCATE(force(1:space%dim, 1:natoms))
    SAFE_ALLOCATE(force_components(1:space%dim, 1:natoms, ION_NUM_COMPONENTS))
    
    call ion_interaction_calculate(ion_interaction, space, latt, atom, natoms, catom, ncatoms, lsize, energy, force, &
      energy_components = energy_components, force_components = force_components)

    call messages_write('Ionic energy        =')
    call messages_write(energy, fmt = '(f20.10)')
    call messages_info()

    call messages_write('Real space energy   =')
    call messages_write(energy_components(ION_COMPONENT_REAL), fmt = '(f20.10)')
    call messages_info()

    call messages_write('Self energy         =')
    call messages_write(energy_components(ION_COMPONENT_SELF), fmt = '(f20.10)')
    call messages_info()

    call messages_write('Fourier energy      =')
    call messages_write(energy_components(ION_COMPONENT_FOURIER), fmt = '(f20.10)')
    call messages_info()
    
    call messages_info()

    do iatom = 1, natoms
      call messages_write('Ionic force         atom')
      call messages_write(iatom)
      call messages_write(' =')
      do idir = 1, space%dim
        call messages_write(force(idir, iatom), fmt = '(f20.10)')
      end do
      call messages_info()

      call messages_write('Real space force    atom')
      call messages_write(iatom)
      call messages_write(' =')
      do idir = 1, space%dim
        call messages_write(force_components(idir, iatom, ION_COMPONENT_REAL), fmt = '(f20.10)')
      end do
      call messages_info()

      call messages_write('Fourier space force atom')
      call messages_write(iatom)
      call messages_write(' =')
      do idir = 1, space%dim
        call messages_write(force_components(idir, iatom, ION_COMPONENT_FOURIER), fmt = '(f20.10)')
      end do
      call messages_info()

      call messages_info()
    end do

    SAFE_DEALLOCATE_A(force)
    SAFE_DEALLOCATE_A(force_components)
    
    call ion_interaction_end(ion_interaction)
    
    POP_SUB(ion_interaction_test)
  end subroutine ion_interaction_test
    
end module ion_interaction_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
