!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module system
use global
use mesh
use geometry
use states

implicit none

type system_type
  FLOAT                      :: val_charge ! the charge of the valence electrons (necessary to initialize states)
  type(mesh_type)            :: m          ! the mesh
  type(geometry_type)        :: geo        ! the geometry
  type(states_type), pointer :: st         ! the states
  type(output_type)          :: outp
end type system_type

contains

  subroutine system_init(s)
    type(system_type), intent(out) :: s

    FLOAT :: def_h, def_rsize

    call push_sub('system_init')
    
    ! initialize the stuff related to the mesh
    call geometry_init_xyz(s%geo)
    call geometry_init_species(s%geo, s%val_charge, def_h, def_rsize)
    call mesh_init(s%m, s%geo, def_h, def_rsize)
    call functions_init(s%m)
    
    ! initialize the other stuff
    allocate(s%st)
    call states_init(s%st, s%m, s%val_charge, s%geo%nlcc)
    call output_init(s%outp)
    
    call pop_sub()
  end subroutine system_init

  subroutine system_end(s)
    type(system_type), intent(inout) :: s
    
    call push_sub('system_end')
    
    if(associated(s%st)) then
      call states_end(s%st)
      deallocate(s%st); nullify(s%st)
    end if
    call geometry_end(s%geo)
    call functions_end(s%m)
    call mesh_end(s%m)
    
    call pop_sub()
  end subroutine system_end

  !! The reason why the next subroutines are in system is that they need to use
  !! both mesh and atom. Previously they were inside atom, but this leads to
  !! circular dependencies of the modules

  subroutine atom_get_wf(m, atom, l, lm, ispin, psi)
    type(mesh_type), intent(IN)  :: m             ! the mesh in which the psi is defined
    type(atom_type), intent(IN)  :: atom          !
    integer,         intent(in)  :: l, lm, ispin  ! quantum numbers of psi
    R_TYPE,          intent(out) :: psi(m%np)     ! the atomic wavefunction
    
    integer :: j, ll
    FLOAT :: x(3), a(3), r, p, ylm
    type(loct_spline_type), pointer :: s

    call push_sub('atom_get_wf')
    
    a = atom%x
    if(atom%spec%local) then
      ! add a couple of harmonic oscilator functions
    else
      s => atom%spec%ps%Ur(l, ispin)

      ll = atom%spec%ps%conf%l(l)
      do j = 1, m%np
        call mesh_r(m, j, r, x=x, a=a)
        p = loct_splint(s, r)
        ylm = loct_ylm(x(1), x(2), x(3), ll, lm)
        psi(j) = p * ylm
      end do
    end if
 
    call pop_sub()
  end subroutine atom_get_wf

  function atom_density(m, atom, spin_channels) result(rho)
    type(mesh_type), intent(IN) :: m
    type(atom_type), intent(IN) :: atom
    integer,         intent(in) :: spin_channels
    FLOAT                       :: rho(m%np, spin_channels)

    integer :: opt, i, in_points, k, n
    FLOAT :: r
    R_TYPE :: psi1
    type(specie_type), pointer :: s

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    s => atom%spec

    ! select what kind of density we should build
    if(conf%dim==3) then
      select case(s%label(1:5))
      case('usdef')
        opt = 1
      case('jelli', 'point')
        opt = 2
      case default
        opt = 3
      end select
    else
      opt = 1      
    end if

    rho = M_ZERO
    ! build density...
    select case (opt)
    case (1) ! ...from userdef
      rho(1:m%np, 1:spin_channels) = real(s%Z_val, PRECISION) &
                                     /(m%np*m%vol_pp*real(spin_channels, PRECISION))
    case (2) ! ...from jellium
      in_points = 0
      do i = 1, m%np
        call mesh_r(m, i, r, a=atom%x)
        if(r <= s%jradius) then
          in_points = in_points + 1
        end if
      end do
    
      if(in_points > 0) then
        do i = 1, m%np
          call mesh_r(m, i, r, a=atom%x)
          if(r <= s%jradius) then
            rho(i, 1:spin_channels) = real(s%z_val, PRECISION)/ &
                                      (real(in_points*spin_channels, PRECISION)*m%vol_pp)
          end if
        end do
      end if

    case (3) ! ...from pseudopotential
      ! the outer loop sums densities over atoms in neighbour cells
      do k = 1,3**conf%periodic_dim
        do i = 1, m%np
          call mesh_r(m, i, r, a=atom%x+m%shift(k,:))
          do n = 1, s%ps%conf%p
            if(r >= r_small) then
              select case(spin_channels)
              case(1)
                psi1 = loct_splint(s%ps%Ur(n, 1), r)
                rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)  
              case(2)
                psi1 = loct_splint(s%ps%Ur(n, 1), r)
                rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
                rho(i, 2) = rho(i, 2) + s%ps%conf%occ(n, 2)*psi1*psi1 /(M_FOUR*M_PI)
              end select
            end if
          end do
        end do
      end do
    end select

  end function atom_density

  ! builds a density which is the sum of the atomic densities
  subroutine system_guess_density(m, geo, qtot, nspin, spin_channels, rho)
    type(mesh_type),     intent(IN)  :: m
    type(geometry_type), intent(IN)  :: geo
    FLOAT,               intent(in)  :: qtot  ! the total charge of the system
    integer,             intent(in)  :: nspin, spin_channels
    FLOAT,               intent(out) :: rho(m%np, nspin)
    
    integer :: ia, is, gd_spin, i
    integer, save :: iseed = 321
    FLOAT :: r, rnd, phi, theta
    FLOAT, allocatable :: atom_rho(:,:)
    
    call push_sub('guess_density')
    
    if (spin_channels == 1) then
      gd_spin = 1
    else
      call loct_parse_int("GuessDensitySpin", 2, gd_spin)
      if (gd_spin < 0 .or. gd_spin > 3) then
        write(message(1),'(a,i1,a)') "Input: '",gd_spin ,"' is not a valid GuessDensitySpin"
        message(2) = '(GuessDensitySpin = 1 | 2 | 3)'
        call write_fatal(2)
      end if
    end if
    
    rho = M_ZERO
    select case (gd_spin)
    case (1) ! Spin-unpolarized
      do ia = 1, geo%natoms
        rho(1:m%np, 1:1) = rho(1:m%np, 1:1) + atom_density(m, geo%atom(ia), 1)
      end do
      if (spin_channels == 2) then
        rho(1:m%np, 1) = M_HALF*rho(1:m%np, 1)
        rho(1:m%np, 2) = rho(1:m%np, 1)
      end if
      
    case (2) ! Spin-polarized
      do ia = 1, geo%natoms
        rho(1:m%np, 1:2) = rho(1:m%np, 1:2) + atom_density(m, geo%atom(ia), 2)
      end do
      
    case (3) ! Random oriented spins
      allocate(atom_rho(m%np, 2))
      do ia = 1, geo%natoms
        atom_rho = atom_density(m, geo%atom(ia), 2)
        
        if (nspin == 2) then
          call quickrnd(iseed, rnd)
          rnd = rnd - M_HALF
          if (rnd > M_ZERO) then
            rho(1:m%np, 1:2) = rho(1:m%np, 1:2) + atom_rho(1:m%np, 1:2)
          else
            rho(1:m%np, 1) = rho(1:m%np, 1) + atom_rho(1:m%np, 2)
            rho(1:m%np, 2) = rho(1:m%np, 2) + atom_rho(1:m%np, 1)
          end if
        elseif (nspin == 4) then
          call quickrnd(iseed, phi)
          call quickrnd(iseed, theta)
          phi = phi*M_TWO*M_PI
          theta = theta*M_PI
          rho(1:m%np, 1) = rho(1:m%np, 1) + cos(theta/M_TWO)**2*atom_rho(1:m%np, 1) &
                                          + sin(theta/M_TWO)**2*atom_rho(1:m%np, 2)
          rho(1:m%np, 2) = rho(1:m%np, 2) + sin(theta/M_TWO)**2*atom_rho(1:m%np, 1) &
                                          + cos(theta/M_TWO)**2*atom_rho(1:m%np, 2)
          rho(1:m%np, 3) = rho(1:m%np, 3) + cos(theta/M_TWO)*sin(theta/M_TWO)*cos(phi)* &
                                            (atom_rho(1:m%np, 1) - atom_rho(1:m%np, 2))
          rho(1:m%np, 4) = rho(1:m%np, 4) - cos(theta/M_TWO)*sin(theta/M_TWO)*sin(phi)* &
                                            (atom_rho(1:m%np, 1) - atom_rho(1:m%np, 2))
        end if
      end do
      
    end select

    ! we now renormalize the density (necessary if we have a charged system)
    r = M_ZERO
    do is = 1, spin_channels
      r = r + dmf_integrate(m, rho(:, is))
    end do

    write(message(1),'(a,f13.6)')'Info: Unnormalized total charge = ', r
    call write_info(1)

    r = qtot/r
    rho = r*rho
    r = M_ZERO
    do is = 1, spin_channels
      r = r + dmf_integrate(m, rho(:, is))
    end do

    write(message(1),'(a,f13.6)')'Info: Renormalized total charge = ', r
    call write_info(1)
    
    call pop_sub()

  end subroutine system_guess_density
  
end module system
