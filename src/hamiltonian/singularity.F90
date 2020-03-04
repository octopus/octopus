!! Copyright (C) 2019 N. Tancogne-Dejean
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

module singularity_oct_m
  use comm_oct_m
  use distributed_oct_m
  use global_oct_m
  use kpoints_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                          &
    singularity_t,             &
    singularity_nullify,       &
    singularity_init,          &
    singularity_end

  integer, public, parameter ::        &
    SINGULARITY_NONE                  = 0, &
    SINGULARITY_GENERAL               = 1, &
    SINGULARITY_GYGI                  = 2, &
    SINGULARITY_SPHERE                = 3

   
  type singularity_t
    !For the treatment of the singularity in solids
    integer :: coulomb_singularity
    FLOAT, allocatable :: Fk(:)
    FLOAT :: FF
  end type singularity_t
 
contains

  subroutine singularity_nullify(this)
    type(singularity_t), intent(out) :: this

    PUSH_SUB(singularity_nullify)

    this%coulomb_singularity = 0

    POP_SUB(singularity_nullify)
  end subroutine singularity_nullify
 
  subroutine singularity_init(this, namespace, st, sb)
    type(singularity_t),       intent(inout) :: this
    type(namespace_t),         intent(in)    :: namespace
    type(states_elec_t),       intent(in)    :: st
    type(simul_box_t),         intent(in)    :: sb

    integer :: default

    PUSH_SUB(singularity_init)


    if(.not.allocated(this%Fk)) then
      SAFE_ALLOCATE(this%Fk(st%d%kpt%start:st%d%kpt%end))
      this%Fk(st%d%kpt%start:st%d%kpt%end) = M_ZERO
      this%FF = M_ZERO
    end if

    if(.not.allocated(this%Fk) .and. sb%periodic_dim > 0) then

      !%Variable HFSingularity
      !%Type integer
      !%Default general
      !%Section Hamiltonian::XC
      !%Description
      !% (Experimental) This variable selects the method used for the treatment of the 
      !% singularity of the Coulomb potential in Hatree-Fock and hybrid-functional DFT calculations.
      !% This shoulbe be only applied for periodic systems and is only
      !% used for FFT kernels of the Poisson solvers.
      !%Option none 0
      !% The singularity is replaced by zero.
      !%Option general 1
      !% The general treatment of the singularity, as described in Carrier et al, PRB 75 205126 (2007).
      !% This is the default option
      !%Option fcc 2
      !% The treatment of the singulariy as described in Gygi and Baldereschi, PRB 34, 4405 (1986).
      !% This is formaly valid for cubic systems only.
      !%Option spherical_bz 3
      !% The divergence in q=0 is treated analytically assuming a spherical Brillouin zone
      !%End

      default = SINGULARITY_NONE
      if(sb%dim == 3) default = SINGULARITY_GENERAL
      
      call parse_variable(namespace, 'HFSingularity', default, this%coulomb_singularity)
      call messages_print_var_option(stdout,  'HFSingularity', this%coulomb_singularity)

      if(this%coulomb_singularity /= SINGULARITY_NONE) then
        call singularity_correction(this, namespace, st, sb)
      end if
    end if

    POP_SUB(singularity_init)
  end subroutine singularity_init

  subroutine singularity_end(this)
    type(singularity_t), intent(inout) :: this

    PUSH_SUB(singularity_end)

    this%coulomb_singularity = -1
    
    POP_SUB(singularity_end)
  end subroutine singularity_end

  !This routine implements the general tratment of the singularity for periodic solids,
  !as described in Carrier et al. PRB 75, 205126 (2007)
  subroutine singularity_correction(this, namespace, st, sb)
    type(singularity_t),       intent(inout) :: this
    type(namespace_t),         intent(in)    :: namespace
    type(states_elec_t),       intent(in)    :: st
    type(simul_box_t),         intent(in)    :: sb

    integer :: ik, ik2, ikpoint, Nk, Nsteps
    integer :: ikx, iky, ikz, istep, kpt_start, kpt_end
    FLOAT :: length
    FLOAT :: kpoint(1:MAX_DIM), qpoint(1:MAX_DIM)
    FLOAT :: kvol_element
    FLOAT :: energy
    type(distributed_t) :: dist_kpt
    type(profile_t), save :: prof
    FLOAT, parameter :: SINGUL_CNST = 7.7955541794415 !The constant is 4*pi*(3/(4*pi))^1/3
 
    PUSH_SUB(singularity_correction)

    call profiling_in(prof, "COULOMB_SINGULARITY")

    !At the moment this is only implemented in 3D.
    ASSERT(sb%dim == 3)

    call distributed_nullify(dist_kpt, 0)
    kpt_start = st%d%kpt%start
    kpt_end = st%d%kpt%end

    if(.not.st%d%kpt%parallel) then
      call distributed_init(dist_kpt, st%d%nik, mpi_world%comm, "singularity")
      kpt_start = dist_kpt%start
      kpt_end = dist_kpt%end
    end if

    do ik = kpt_start, kpt_end
      ikpoint = states_elec_dim_get_kpoint_index(st%d, ik)
      kpoint = M_ZERO
      kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, ikpoint, absolute_coordinates = .false.) 

      this%Fk(ik) = M_ZERO

      do ik2 = 1, sb%kpoints%full%npoints
        qpoint = M_ZERO
        qpoint(1:sb%dim) = kpoint(1:sb%dim) - sb%kpoints%full%red_point(1:sb%dim, ik2)
 
        if(all(abs(qpoint(1:sb%dim))< CNST(1e-6))) cycle

        this%Fk(ik) = this%Fk(ik) + aux_funct(qpoint) * sb%kpoints%full%weight(ik2)
      end do
      this%Fk(ik) = this%Fk(ik)*CNST(4.0)*M_PI/sb%rcell_volume
    end do

    if(dist_kpt%parallel) then
      call comm_allreduce(dist_kpt%mpi_grp%comm, this%Fk)
    end if
    call distributed_end(dist_kpt)


    if(this%coulomb_singularity == SINGULARITY_GENERAL) then 
      !%Variable HFSingularityNk
      !%Type integer
      !%Default 60
      !%Section Hamiltonian::XC
      !%Description
      !% Number of k-point used (total number of k-points) is (2*Nk+1)^3) in the numerical integration
      !% of the auxiliary function f(q). See PRB 75, 205126 (2007) for more details. 
      !% Only for HFSingularity=general.
      !%End
      call parse_variable(namespace, 'HFSingularityNk', 60, Nk)
      if(abs(Nk/M_THREE-nint(Nk/M_THREE)) > M_EPSILON) then
        message(1) = 'HFSingularity_Nk must be a multiple of 3.'
        call messages_fatal(1, namespace=namespace)
      end if

      !%Variable HFSingularityNsteps
      !%Type integer
      !%Default 7
      !%Section Hamiltonian::XC
      !%Description
      !% Number of grid refinement steps in the numerical integration of the auxiliary function f(q).
      !% See PRB 75, 205126 (2007) for more details. Only for HFSingularity=general.
      !%End
      call parse_variable(namespace, 'HFSingularityNsteps', 7, Nsteps)

      this%FF = M_ZERO
      length = M_ONE
      kvol_element = (M_ONE/(M_TWO*Nk+M_ONE))**3*((M_TWO*M_PI)**3)/sb%rcell_volume
      do istep = 1, Nsteps

        do ikx = 0, Nk
          qpoint(1) = ikx/(M_TWO*Nk)*length

          do iky = -Nk, Nk
            qpoint(2) = iky/(M_TWO*Nk)*length
 
            do ikz = -Nk, Nk
              qpoint(3) = ikz/(M_TWO*Nk)*length

              if(abs(ikx)<=Nk/3 .and. abs(iky)<=Nk/3 .and. abs(ikz)<=Nk/3) cycle

              this%FF = this%FF + aux_funct(qpoint)*kvol_element
            end do
          end do
        end do
        if(istep<Nsteps) then
          length = length/M_THREE
          kvol_element = kvol_element / CNST(27.0)
        end if
      end do

      !We have a factor two because we used the fact that f(q)=f(-q)
      !We multiply by 4*pi/((2*pi)^3)
      this%FF = this%FF*CNST(8.0)*M_PI/((M_TWO*M_PI)**3)
      !The remaining part is treated as a spherical BZ
      this%FF = this%FF + SINGUL_CNST*(sb%rcell_volume)**(CNST(2.0/3.0))/M_PI/sb%rcell_volume*length

    else if(this%coulomb_singularity == SINGULARITY_GYGI) then
      !See Eq. (7) of PRB 34, 4405 (1986)
      !Here we use the fact that the fcc volume is a^3/4
      this%FF = CNST(4.423758)*(sb%rcell_volume*CNST(4.0))**(CNST(2.0/3.0))/M_PI/sb%rcell_volume

    else
      !The constant is 4*pi*(3/(4*pi))^1/3
      !We multiply by 4*pi/(2*pi^3)
      this%FF = SINGUL_CNST*(sb%rcell_volume)**(CNST(2.0/3.0))/M_PI/sb%rcell_volume
    end if

    if(debug%info) then
      energy = M_ZERO
      do ik = st%d%kpt%start, st%d%kpt%end
        energy = energy + this%Fk(ik)*st%d%kweights(ik)
      end do

      if(st%d%kpt%parallel) then
        call comm_allreduce(st%d%kpt%mpi_grp%comm, energy) 
      end if

      write(message(1), '(a,f12.6,a,a,a)') 'Debug: Singularity energy ', &
              units_from_atomic(units_out%energy, (energy-this%FF)*st%qtot/st%smear%el_per_state), &
              ' [',trim(units_abbrev(units_out%energy)),']'
      call messages_info(1)
    end if

    call profiling_out(prof)
    POP_SUB(singularity_correction)

  contains
    
    FLOAT function aux_funct(qq) result(ff)
      FLOAT,   intent(in) :: qq(1:MAX_DIM)
     
      FLOAT :: half_a, qq_abs(1:MAX_DIM)

      PUSH_SUB(singularity_correction.aux_funct)

      if(this%coulomb_singularity == SINGULARITY_GENERAL) then
        !See Eq. (16) of PRB 75, 205126 (2007)
        ff = (M_TWO*M_PI)**2/(M_TWO*(                                                              &
           (M_TWO*sin(qq(1)*M_PI)*sin(qq(1)*M_PI)*dot_product(sb%klattice(1:3,1),sb%klattice(1:3,1))  &
         +sin(qq(1)*M_TWO*M_PI)*sin(qq(2)*M_TWO*M_PI)*dot_product(sb%klattice(1:3,1),sb%klattice(1:3,2))) &
          +(M_TWO*sin(qq(2)*M_PI)*sin(qq(2)*M_PI)*dot_product(sb%klattice(1:3,2),sb%klattice(1:3,2))  &
         +sin(qq(2)*M_TWO*M_PI)*sin(qq(3)*M_TWO*M_PI)*dot_product(sb%klattice(1:3,2),sb%klattice(1:3,3))) &
          +(M_TWO*sin(qq(3)*M_PI)*sin(qq(3)*M_PI)*dot_product(sb%klattice(1:3,3),sb%klattice(1:3,3))  &
         +sin(qq(3)*M_TWO*M_PI)*sin(qq(1)*M_TWO*M_PI)*dot_product(sb%klattice(1:3,3),sb%klattice(1:3,1)))))
      else
        half_a = M_HALF*(sb%rcell_volume*CNST(4.0))**(CNST(1.0/3.0))
        call kpoints_to_absolute(sb%klattice, qq, qq_abs, 3)
        !See Eq. (6) of PRB 34, 4405 (1986)
        ff = (half_a)**2/(M_THREE-cos(qq_abs(1)*half_a)*cos(qq_abs(2)*half_a) &
                            -cos(qq_abs(1)*half_a)*cos(qq_abs(3)*half_a)         &
                            -cos(qq_abs(3)*half_a)*cos(qq_abs(2)*half_a))
      end if

      POP_SUB(singularity_correction.aux_funct)
    end function aux_funct
  end subroutine singularity_correction

end module singularity_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
