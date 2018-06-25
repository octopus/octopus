!! Copyright (C) 2015 X. Andrade
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

module hirshfeld_oct_m
  use derivatives_oct_m
  use messages_oct_m
  use geometry_oct_m
  use global_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use periodic_copy_oct_m
  use profiling_oct_m
  use ps_oct_m
  use species_pot_oct_m
  use states_oct_m
  use species_oct_m
  use splines_oct_m
  use mpi_oct_m !for writing
  use parser_oct_m !for the variable
 
  implicit none

  private
  public ::                       &
    hirshfeld_t,                  &
    hirshfeld_init,               &
    hirshfeld_end,                &
    hirshfeld_charge,             &
    hirshfeld_volume_ratio,       &
    hirshfeld_density_derivative, &
    hirshfeld_position_derivative
  
  type hirshfeld_t

    private
    type(mesh_t),     pointer     :: mesh
    type(geometry_t), pointer     :: geo
    type(states_t),   pointer     :: st
    FLOAT,            pointer     :: total_density(:)  !< (mesh%np)
    FLOAT,            pointer     :: free_volume(:)    !< (natoms)
    FLOAT,            pointer     :: free_vol_r3(:,:)  !< (natoms,mesh%np)

  end type hirshfeld_t

  real(8), parameter, public :: TOL_HIRSHFELD = CNST(1e-12)
    
contains

  subroutine hirshfeld_init(this, mesh, geo, st)
    type(hirshfeld_t),         intent(out)   :: this
    type(mesh_t),      target, intent(in)    :: mesh
    type(geometry_t),  target, intent(in)    :: geo
    type(states_t),    target, intent(in)    :: st
    
    integer :: iatom, ip, isp
    FLOAT :: rr, pos(1:MAX_DIM), rmax
    FLOAT, allocatable :: atom_density(:, :), atom_density_acc(:)
    type(ps_t), pointer :: ps
    type(periodic_copy_t) :: pp
    INTEGER :: icell
    type(profile_t), save :: prof
    PUSH_SUB(hirshfeld_init)

    call profiling_in(prof, "HIRSHFELD_INIT")

    this%mesh => mesh
    this%geo  => geo
    this%st   => st

    SAFE_ALLOCATE(this%total_density(1:this%mesh%np))
    SAFE_ALLOCATE(this%free_volume(1:geo%natoms))
    SAFE_ALLOCATE(this%free_vol_r3(1:geo%natoms,1:this%mesh%np))
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    SAFE_ALLOCATE(atom_density_acc(1:this%mesh%np))

    this%total_density = CNST(0.0)
  
    do iatom = 1, geo%natoms
      ps => species_ps(this%geo%atom(iatom)%species)
      atom_density_acc(1:this%mesh%np) = M_ZERO

      rmax = CNST(0.0)
      do isp = 1, this%st%d%nspin
        rmax = max(rmax, spline_cutoff_radius(ps%density(isp), ps%projectors_sphere_threshold))
      end do

      call periodic_copy_init(pp, this%mesh%sb, this%geo%atom(iatom)%x, rmax)

      do icell = 1, periodic_copy_num(pp)
        pos(1:this%mesh%sb%dim) = periodic_copy_position(pp, this%mesh%sb, icell) 
        !We get the non periodized density
        !We need to do it to have the r^3 correctly computed for periodic systems
        call species_atom_density_np(this%mesh, this%mesh%sb, this%geo%atom(iatom), pos, this%st%d%nspin, atom_density)

        forall(ip = 1:this%mesh%np) this%total_density(ip) = this%total_density(ip) + sum(atom_density(ip, 1:st%d%nspin))
      
        do ip = 1, this%mesh%np
          rr = sqrt(sum((this%mesh%x(ip, 1:this%mesh%sb%dim) - pos(1:this%mesh%sb%dim))**2))
          atom_density_acc(ip) = atom_density_acc(ip) + sum(atom_density(ip, 1:this%st%d%nspin))*rr**3  
        end do
      end do
      call periodic_copy_end(pp)
      this%free_volume(iatom) = dmf_integrate(this%mesh, atom_density_acc(:))
      this%free_vol_r3(iatom,:) = atom_density_acc(:)
    end do

    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(atom_density_acc)

    call profiling_out(prof)

    POP_SUB(hirshfeld_init)    
  end subroutine hirshfeld_init

  ! ------------------------------------------------
  
  subroutine hirshfeld_end(this)
    type(hirshfeld_t), intent(inout) :: this

    PUSH_SUB(hirshfeld_end)

    SAFE_DEALLOCATE_P(this%total_density)
    SAFE_DEALLOCATE_P(this%free_volume)
    SAFE_DEALLOCATE_P(this%free_vol_r3)

    nullify(this%mesh)
    nullify(this%geo)
    nullify(this%st)
    
    POP_SUB(hirshfeld_end)    
  end subroutine hirshfeld_end

  ! -----------------------------------------------
  
  subroutine hirshfeld_charge(this, iatom, density, charge)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(in)    :: density(:, :)
    FLOAT,                     intent(out)   :: charge

    integer :: ip
    FLOAT :: dens_ip
    FLOAT, allocatable :: atom_density(:, :), hirshfeld_density(:)
    type(profile_t), save :: prof
    
    PUSH_SUB(hirshfeld_charge)

    call profiling_in(prof, "HIRSHFELD_CHARGE")

    ASSERT(associated(this%total_density))
    
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    SAFE_ALLOCATE(hirshfeld_density(1:this%mesh%np))
    
    call species_atom_density(this%mesh, this%mesh%sb, this%geo%atom(iatom), this%st%d%nspin, atom_density)

    do ip = 1, this%mesh%np
      dens_ip = sum(atom_density(ip, 1:this%st%d%nspin))
      if(abs(dens_ip) > TOL_HIRSHFELD) then
        hirshfeld_density(ip) = sum(density(ip, 1:this%st%d%nspin))*dens_ip/this%total_density(ip)
      else
        hirshfeld_density(ip) = CNST(0.0)
      end if
    end do

    charge = dmf_integrate(this%mesh, hirshfeld_density)
    
    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(hirshfeld_density)

    call profiling_out(prof)
    
    POP_SUB(hirshfeld_charge)
  end subroutine hirshfeld_charge

  ! -----------------------------------------------    

  subroutine hirshfeld_volume_ratio(this, iatom, density, volume_ratio)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(in)    :: density(:, :)
    FLOAT,                     intent(out)   :: volume_ratio

    integer :: ip
    FLOAT :: dens_ip, rr
    FLOAT, allocatable :: atom_density(:, :), hirshfeld_density(:)

    type(profile_t), save :: prof
    
    PUSH_SUB(hirshfeld_volume_ratio)

    call profiling_in(prof, "HIRSHFELD_VOLUME_RATIO")

    ASSERT(associated(this%total_density))
    
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    SAFE_ALLOCATE(hirshfeld_density(1:this%mesh%np))
    

    do ip = 1, this%mesh%np
      if( this%total_density(ip) > TOL_HIRSHFELD) then
        hirshfeld_density(ip) = this%free_vol_r3(iatom,ip)*sum(density(ip, 1:this%st%d%nspin))/this%total_density(ip)
      else
        hirshfeld_density(ip) = CNST(0.0)
      end if
    end do
    
    volume_ratio = dmf_integrate(this%mesh, hirshfeld_density)/this%free_volume(iatom)
    
    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(hirshfeld_density)

    call profiling_out(prof)
    
    POP_SUB(hirshfeld_volume_ratio)
  end subroutine hirshfeld_volume_ratio
  
  ! -----------------------------------------------

  subroutine hirshfeld_density_derivative(this, iatom, ddensity)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(out)   :: ddensity(:)

    integer :: ip
    FLOAT :: dens_ip, rr
    FLOAT, allocatable :: atom_density(:, :)
    type(profile_t), save :: prof   
 
    PUSH_SUB(hirshfeld_density_derivative)

    call profiling_in(prof, "HIRSHFELD_DENSITY_DER")

    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    
    do ip = 1, this%mesh%np

      if(abs(this%total_density(ip)) > TOL_HIRSHFELD) then
        ddensity(ip) = this%free_vol_r3(iatom,ip)/(this%total_density(ip)*this%free_volume(iatom))
      else
        ddensity(ip) = CNST(0.0)
      end if
    end do

    SAFE_DEALLOCATE_A(atom_density)
    
    call profiling_out(prof)
    
    POP_SUB(hirshfeld_density_derivative)
  end subroutine hirshfeld_density_derivative

  ! -----------------------------------------------

  subroutine hirshfeld_position_derivative(this, der, iatom, jatom, density, dposition)
    type(hirshfeld_t),         intent(in)    :: this
    type(derivatives_t),       intent(in)    :: der
    integer,                   intent(in)    :: iatom
    integer,                   intent(in)    :: jatom
    FLOAT,                     intent(in)    :: density(:, :)
    FLOAT,                     intent(out)   :: dposition(:)

    integer :: ip, idir, icell, jcell, isp, ipp
    FLOAT :: atom_dens, atom_der,rri, rrj, tdensity, pos_i(1:MAX_DIM), pos_j(1:MAX_DIM), rmax_i, rmax_j, &
             VDW_TS_extraterm_i, VDW_TS_extraterm_j, ttt, rij
    FLOAT, allocatable :: grad(:, :), atom_density(:, :), atom_derivative(:, :)
    type(periodic_copy_t) :: pp_i, pp_j
    type(ps_t), pointer :: ps_i, ps_j
    type(profile_t), save :: prof, prof2
    FLOAT :: tmp, xxi(1:MAX_DIM), xxj(1:MAX_DIM)

    PUSH_SUB(hirshfeld_position_derivative)


    call profiling_in(prof, "HIRSHFELD_POSITION_DER")

    if(this%mesh%sb%periodic_dim > 0) then ! periodic case

      !if(mpi_grp_is_root(mpi_world)) then
        !call profiling_in(prof, "HIRSHFELD_POSITION_DER")

        SAFE_ALLOCATE(grad(1:this%mesh%np, 1:this%mesh%sb%dim))
        SAFE_ALLOCATE(atom_derivative(1:this%mesh%np, 1:this%st%d%nspin))
        SAFE_ALLOCATE(atom_density(1:this%mesh%np, 1:this%st%d%nspin))

        dposition(1:this%mesh%sb%dim) = M_ZERO
        grad(1:this%mesh%np, 1:this%mesh%sb%dim) = M_ZERO

        ps_i => species_ps(this%geo%atom(iatom)%species)
        ps_j => species_ps(this%geo%atom(jatom)%species)

        rmax_i = CNST(0.0)
        rmax_j = CNST(0.0)
        do isp = 1, this%st%d%nspin
          rmax_i = max(rmax_i, spline_cutoff_radius(ps_i%density(isp), ps_i%projectors_sphere_threshold))
          rmax_j = max(rmax_j, spline_cutoff_radius(ps_j%density_der(isp), ps_j%projectors_sphere_threshold))
        end do

        !if(mpi_grp_is_root(mpi_world)) then
        !  print*, 'rmaxi rmaxj', rmax_i, rmax_j
        !end if


        !%Variable VDW_TS_extraterm_i
        !%Type float
        !%Default 1
        !%Section Hamiltonian::XC
        !%Description
        !% to be converged  
        !%End
        call parse_variable('VDW_TS_extraterm_i', CNST(1.0), VDW_TS_extraterm_i)

        !%Variable VDW_TS_extraterm_j
        !%Type float
        !%Default 1
        !%Section Hamiltonian::XC
        !%Description
        !% to be converged                                                                         
        !%End
        call parse_variable('VDW_TS_extraterm_j', CNST(1.0), VDW_TS_extraterm_j)
    
        call periodic_copy_init(pp_i, this%mesh%sb, this%geo%atom(iatom)%x, VDW_TS_extraterm_i*rmax_i)

        do icell = 1, periodic_copy_num(pp_i)
          atom_density(1:this%mesh%np, 1:this%st%d%nspin) = M_ZERO
          pos_i(1:this%mesh%sb%dim) = periodic_copy_position(pp_i, this%mesh%sb, icell)

          !print*, 'icell, pos_i, pos initial', icell, pos_i(1:this%mesh%sb%dim),  this%geo%atom(iatom)%x

          !We get the non periodized density
          !We need to do it to have the r^3 correctly computed for periodic systems
          !call species_atom_density_np(this%mesh, this%mesh%sb, this%geo%atom(iatom), &
          !       this%geo%atom(iatom)%x, this%st%d%nspin, atom_density(1:this%mesh%np, 1:this%st%d%nspin))
          call species_atom_density_np(this%mesh, this%mesh%sb, this%geo%atom(iatom), &
                 pos_i, this%st%d%nspin, atom_density(1:this%mesh%np, 1:this%st%d%nspin))

          !ttt = 0
          !do ip = 1, this%mesh%np
          !  ttt = max(ttt, abs(sum(atom_density(ip, 1:this%st%d%nspin))))
          !end do

          !print*, 'condition passed, icell:', icell,  ttt, abs(sum(atom_density(1:this%mesh%np, 1:this%st%d%nspin)))
          !if(ttt < 0.0000000000000001) cycle
          
          call periodic_copy_init(pp_j, this%mesh%sb, pos_i, VDW_TS_extraterm_j*(rmax_j+rmax_i))

          !call profiling_in(prof2, "HIRSHFELD_POSITION_DER2")

          do jcell = 1, periodic_copy_num(pp_j)

            pos_j(1:this%mesh%sb%dim) = periodic_copy_position(pp_j, this%mesh%sb, jcell) + (this%geo%atom(jatom)%x(1:this%mesh%sb%dim)-this%geo%atom(iatom)%x(1:this%mesh%sb%dim))
            rij =  sqrt(sum((pos_i(1:this%mesh%sb%dim)-pos_j(1:this%mesh%sb%dim))**2))

            if(rij < 1.00000000001*(rmax_j+rmax_i)) then
             
              print*, ' icell, jcell, pos_i, pos_j:', icell, jcell, pos_i, pos_j

              atom_derivative(1:this%mesh%np, 1:this%st%d%nspin) = M_ZERO
              call species_atom_density_derivative_np(this%mesh, this%mesh%sb, this%geo%atom(jatom), &
                   pos_j, this%st%d%spin_channels, &
                   atom_derivative(1:this%mesh%np, 1:this%st%d%nspin))
              !ttt = 0
              !do ip = 1, this%mesh%np
              !  ttt= max(ttt, abs(sum(atom_derivative(ip, 1:this%st%d%nspin))))
              !end do          
              !if(mpi_grp_is_root(mpi_world)) then
              !  print*, 'contribution de la celulle jcell', ttt
              !end if

              !if(ttt < 0.00000000001) cycle
              !print*, 'condition passed, icell, jcell:', icell, jcell

              do ip = 1, this%mesh%np
                if(this%total_density(ip)< TOL_HIRSHFELD) cycle

                xxi(1:this%mesh%sb%dim) = this%mesh%x(ip, 1:this%mesh%sb%dim) - pos_i(1:this%mesh%sb%dim)
                rri = sqrt(sum(xxi(1:this%mesh%sb%dim)**2))
                tdensity = sum(density(ip, 1:this%st%d%nspin))
                atom_dens = sum(atom_density(ip, 1:this%st%d%nspin))

                tmp = rri**3*atom_dens*tdensity/this%total_density(ip)**2

                xxj(1:this%mesh%sb%dim) = this%mesh%x(ip, 1:this%mesh%sb%dim) - pos_j(1:this%mesh%sb%dim)
                rrj = sqrt(sum(xxj(1:this%mesh%sb%dim)**2))
                atom_der = sum(atom_derivative(ip, 1:this%st%d%nspin))

                if(rrj > TOL_HIRSHFELD) then
                  do idir = 1, this%mesh%sb%dim
                    grad(ip, idir) = grad(ip, idir) - tmp*atom_der*xxj(idir)/rrj
                  end do
                end if

                if(iatom == jatom .and. sum(abs(pos_i(1:this%mesh%sb%dim)- pos_j(1:this%mesh%sb%dim))) < 0.000001) then
                  !Only if we really have the same atoms
                  !if(all(abs(pos_i(1:this%mesh%sb%dim)-this%geo%atom(iatom)%x(1:this%mesh%sb%dim)) < TOL_HIRSHFELD)) then
                  print*, 'meme atomes'
                  do idir = 1, this%mesh%sb%dim
                    grad(ip, idir) = grad(ip, idir) + (CNST(3.0)*rri*atom_dens + rri**2*atom_der)&
                                        *tdensity/this%total_density(ip)*xxi(idir)
                  end do
                  !end if
                end if
              end do
           ! else 
           !   print*, 'impossible car rij trop grand'
           !   atom_derivative(1:this%mesh%np, 1:this%st%d%nspin) = M_ZERO
           !   call species_atom_density_derivative_np(this%mesh, this%mesh%sb, this%geo%atom(jatom), &
           !        pos_j, this%st%d%spin_channels, &
           !        atom_derivative(1:this%mesh%np, 1:this%st%d%nspin))

           !   do ip = 1, this%mesh%np
           !     if(this%total_density(ip)< TOL_HIRSHFELD) cycle

           !     xxi(1:this%mesh%sb%dim) = this%mesh%x(ip, 1:this%mesh%sb%dim) - pos_i(1:this%mesh%sb%dim)
           !     rri = sqrt(sum(xxi(1:this%mesh%sb%dim)**2))
           !     tdensity = sum(density(ip, 1:this%st%d%nspin))
           !     atom_dens = sum(atom_density(ip, 1:this%st%d%nspin))

           !     tmp = rri**3*atom_dens*tdensity/this%total_density(ip)**2

           !     xxj(1:this%mesh%sb%dim) = this%mesh%x(ip, 1:this%mesh%sb%dim) - pos_j(1:this%mesh%sb%dim)
           !     rrj = sqrt(sum(xxj(1:this%mesh%sb%dim)**2))
           !     atom_der = sum(atom_derivative(ip, 1:this%st%d%nspin))


           !     if(atom_der > 0.0000000000001 .AND. atom_dens > 0.00000000001) then
           !       print*, 'probleme'
           !    end if

           !    if(rrj > TOL_HIRSHFELD) then
           !      do idir = 1, this%mesh%sb%dim
           !     !!!grad(ip, idir) = grad(ip, idir) - tmp*atom_der*xxj(idir)/rrj
           !      end do
           !    end if

           !     !if(iatom == jatom .and. icell == jcell) then
           !    if(iatom == jatom .and. sum(abs(pos_i(1:this%mesh%sb%dim)- pos_j(1:this%mesh%sb%dim))) < 0.000001) then

           !       !Only if we really have the same atoms
           !      if(all(abs(pos_i(1:this%mesh%sb%dim)-this%geo%atom(iatom)%x(1:this%mesh%sb%dim)) < TOL_HIRSHFELD)) then
           !        do idir = 1, this%mesh%sb%dim
           !        !!!grad(ip, idir) = grad(ip, idir) + (CNST(3.0)*rri*atom_dens + rri**2*atom_der)&
           !        !!!                    *tdensity/this%total_density(ip)*xxi(idir)
           !        end do
           !      end if
           !    end if
           !  end do
            end if
          end do
          call periodic_copy_end(pp_j)
        end do
        call periodic_copy_end(pp_i)
        !print *,'calculate the integral'
        do idir = 1, this%mesh%sb%dim
          !dposition(idir) = dmf_integrate(this%mesh, grad(1:this%mesh%np, idir))/this%free_volume(iatom)
          dposition(idir) = sum(grad(1:this%mesh%np, idir))/this%free_volume(iatom)
        end do

        SAFE_DEALLOCATE_A(atom_density)
        SAFE_DEALLOCATE_A(atom_derivative)
        SAFE_DEALLOCATE_A(grad)

        !call profiling_out(prof)
      !end if

    else ! Non periodic case  

      SAFE_ALLOCATE(grad(1:this%mesh%np, 1:this%mesh%sb%dim))

      grad(1:this%mesh%np, 1:this%mesh%sb%dim) = M_ZERO
      dposition(1:this%mesh%sb%dim) = CNST(0.0)
      grad(1:this%mesh%np, 1:this%mesh%sb%dim) = M_ZERO
      ps_i => species_ps(this%geo%atom(iatom)%species)
      ps_j => species_ps(this%geo%atom(jatom)%species)

      rmax_i = CNST(0.0)
      rmax_j = CNST(0.0)
      do isp = 1, this%st%d%nspin
        rmax_i = max(rmax_i, spline_cutoff_radius(ps_i%density(isp), ps_i%projectors_sphere_threshold))
        rmax_j = max(rmax_j, spline_cutoff_radius(ps_j%density_der(isp), ps_j%projectors_sphere_threshold))
      end do
      
      pos_i(1:this%mesh%sb%dim) = this%geo%atom(iatom)%x(1:this%mesh%sb%dim)
      pos_j(1:this%mesh%sb%dim) = this%geo%atom(jatom)%x(1:this%mesh%sb%dim)
      rij =  sqrt(sum((pos_i(1:this%mesh%sb%dim)-pos_j(1:this%mesh%sb%dim))**2))

      if(rij < 1*(rmax_i+rmax_j)) then ! Overlap between atom_derivative and atom_density

        SAFE_ALLOCATE(atom_density(1:this%mesh%np, 1:this%st%d%nspin))
        SAFE_ALLOCATE(atom_derivative(1:this%mesh%np, 1:this%st%d%nspin))

        atom_derivative(1:this%mesh%np, 1:this%st%d%nspin) = M_ZERO
        atom_density(1:this%mesh%np, 1:this%st%d%nspin) = M_ZERO
        grad(1:this%mesh%np, 1:this%mesh%sb%dim) = M_ZERO


        !We get the non periodized density
        !We need to do it to have the r^3 correctly computed for periodic systems
        call species_atom_density_derivative_np(this%mesh, this%mesh%sb, this%geo%atom(jatom), &
                                    pos_j, this%st%d%spin_channels, &
                                    atom_derivative(1:this%mesh%np, 1:this%st%d%nspin))
     

        !We get the non periodized density
        !We need to do it to have the r^3 correctly computed for periodic systems
        call species_atom_density_np(this%mesh, this%mesh%sb, this%geo%atom(iatom), &
               pos_i, this%st%d%nspin, atom_density)


        !ttt = 0
        !do ip = 1, this%mesh%np
        !  ttt = max(ttt, sum(atom_density(ip, 1:this%st%d%nspin)))
        !end do

        !if(mpi_grp_is_root(mpi_world)) then
        !  print*, ttt
        !end if


        !call profiling_in(prof2, "HIRSHFELD_POSITION_DER2")
        do ip = 1, this%mesh%np
          if(this%total_density(ip)< TOL_HIRSHFELD) cycle

          xxi(1:this%mesh%sb%dim) = this%mesh%x(ip, 1:this%mesh%sb%dim) - pos_i(1:this%mesh%sb%dim)
          rri = sqrt(sum(xxi(1:this%mesh%sb%dim)**2))
          tdensity = sum(density(ip, 1:this%st%d%nspin))
          atom_dens = sum(atom_density(ip, 1:this%st%d%nspin))

          tmp = rri**3*atom_dens*tdensity/this%total_density(ip)**2

          xxj(1:this%mesh%sb%dim) = this%mesh%x(ip, 1:this%mesh%sb%dim) - pos_j(1:this%mesh%sb%dim)
          rrj = sqrt(sum(xxj(1:this%mesh%sb%dim)**2))
          atom_der = sum(atom_derivative(ip, 1:this%st%d%nspin))

          if(rrj > TOL_HIRSHFELD) then
            do idir = 1, this%mesh%sb%dim
              grad(ip, idir) = grad(ip, idir) - tmp*atom_der*xxj(idir)/rrj
            end do
          end if

          if(iatom == jatom)then
            !Only if we really have the same atoms
            !if(all(abs(pos_i(1:this%mesh%sb%dim)-this%geo%atom(iatom)%x(1:this%mesh%sb%dim)) < TOL_HIRSHFELD)) then
            do idir = 1, this%mesh%sb%dim
              grad(ip, idir) = grad(ip, idir) + (CNST(3.0)*rri*atom_dens + rri**2*atom_der)&
                                    *tdensity/this%total_density(ip)*xxi(idir)
            end do
            !end if
          end if
        end do

        !call profiling_out(prof2)

        SAFE_DEALLOCATE_A(atom_density)
        SAFE_DEALLOCATE_A(atom_derivative)
      end if
 
      do idir = 1, this%mesh%sb%dim
        dposition(idir) = dmf_integrate(this%mesh, grad(1:this%mesh%np, idir))/this%free_volume(iatom)
      end do

      SAFE_DEALLOCATE_A(grad)
    end if

    call profiling_out(prof)

    POP_SUB(hirshfeld_position_derivative)
  end subroutine hirshfeld_position_derivative
  
end module hirshfeld_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
