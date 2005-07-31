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
!!
!! $Id$

#include "global.h"

module system
  use global
  use messages
  use syslabels
  use lib_oct
  use lib_oct_parser
  use lib_oct_gsl_spline
  use lib_basic_alg
  use math
  use mesh
  use mesh_function
  use functions
  use grid
  use specie
  use geometry
  use states
  use v_ks
  use hamiltonian
  use output

  implicit none

  private
  public :: system_type,     &
       system_init,          &
       system_end,           &
       atom_get_wf,          &
       atom_density,         &
       system_guess_density, &
       dsystem_h_setup,      &
       zsystem_h_setup


  type system_type
     type(grid_type),     pointer :: gr         ! the mesh
     type(states_type),   pointer :: st         ! the states
     type(v_ks_type)              :: ks         ! the Kohn-Sham potentials
     type(output_type)            :: outp       ! the output
  end type system_type

contains

  !----------------------------------------------------------
  subroutine system_init(s)
    type(system_type), intent(out) :: s

    call push_sub('system_init')

    allocate(s%gr)
    call grid_init(s%gr)

    ! initialize the other stuff
    allocate(s%st)
    call states_init(s%st, s%gr)
    call output_init(s%gr%sb, s%outp)

    call v_ks_init(s%gr, s%ks, s%st%d)

    call pop_sub()
  end subroutine system_init


  !----------------------------------------------------------
  subroutine system_end(s)
    type(system_type), intent(inout) :: s

    call push_sub('system_end')

    call v_ks_end(s%ks)

    if(associated(s%st)) then
       call states_end(s%st)
       deallocate(s%st); nullify(s%st)
    end if

    call grid_end(s%gr)
    deallocate(s%gr);  nullify(s%gr)

    call pop_sub()
  end subroutine system_end

  !! The reason why the next subroutines are in system is that they need to use
  !! both mesh and atom. Previously they were inside atom, but this leads to
  !! circular dependencies of the modules

  !----------------------------------------------------------
  subroutine atom_get_wf(m, atom, l, lm, ispin, psi)
    type(mesh_type), intent(in)  :: m             ! the mesh in which the psi is defined
    type(atom_type), intent(in)  :: atom          !
    integer,         intent(in)  :: l, lm, ispin  ! quantum numbers of psi
    R_TYPE,          intent(out) :: psi(:)        ! psi(m%np) the atomic wavefunction

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


  !----------------------------------------------------------
  subroutine atom_density(m, sb, atom, spin_channels, rho)
    type(mesh_type),      intent(in)    :: m
    type(simul_box_type), intent(in)    :: sb
    type(atom_type),      intent(in)    :: atom
    integer,              intent(in)    :: spin_channels
    FLOAT,                intent(inout) :: rho(:, :) ! (m%np, spin_channels)

    integer :: i, in_points, k, n
    FLOAT :: r
    R_TYPE :: psi1, psi2
    type(specie_type), pointer :: s

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    s => atom%spec
    rho = M_ZERO

    ! build density ...
    select case (s%type)
    case (SPEC_USDEF) ! ... from userdef
       do i = 1, spin_channels
          rho(1:m%np, i) = real(s%Z_val, PRECISION) /  &
               (m%np*m%vol_pp(1:m%np)*real(spin_channels, PRECISION))
       end do

    case (SPEC_POINT, SPEC_JELLI) ! ... from jellium
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
                rho(i, 1:spin_channels) = real(s%z_val, PRECISION) /   &
                     (m%vol_pp(i)*real(in_points*spin_channels, PRECISION))
             end if
          end do
       end if

    case (SPEC_PS_TM2,SPEC_PS_HGH) ! ...from pseudopotential
       ! the outer loop sums densities over atoms in neighbour cells
       do k = 1, 3**sb%periodic_dim
          do i = 1, m%np
             call mesh_r(m, i, r, a=atom%x + sb%shift(k,:))
             r = max(r, r_small)
             do n = 1, s%ps%conf%p
                select case(spin_channels)
                case(1)
                   psi1 = loct_splint(s%ps%Ur(n, 1), r)
                   rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
                case(2)
                   psi1 = loct_splint(s%ps%Ur(n, 1), r)
                   psi2 = loct_splint(s%ps%Ur(n, 2), r)
                   rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
                   rho(i, 2) = rho(i, 2) + s%ps%conf%occ(n, 2)*psi2*psi2 /(M_FOUR*M_PI)
                end select
             end do
          end do
       end do
    end select

  end subroutine atom_density


  ! builds a density which is the sum of the atomic densities
  subroutine system_guess_density(m, sb, geo, qtot, nspin, spin_channels, rho)
    type(mesh_type),      intent(in)  :: m
    type(simul_box_type), intent(in)  :: sb
    type(geometry_type),  intent(in)  :: geo
    FLOAT,                intent(in)  :: qtot  ! the total charge of the system
    integer,              intent(in)  :: nspin, spin_channels
    FLOAT,                intent(out) :: rho(:, :)

    integer :: ia, is, gmd_opt, i
    integer, save :: iseed = 321
    integer(POINTER_SIZE) :: blk
    FLOAT :: r, rnd, phi, theta, mag(3)
    FLOAT, allocatable :: atom_rho(:,:)

    call push_sub('guess_density')

    if (spin_channels == 1) then
       gmd_opt = 1
    else
       call loct_parse_int(check_inp('GuessMagnetDensity'), 2, gmd_opt)
       if (gmd_opt < 1 .or. gmd_opt > 4) then
          write(message(1),'(a,i1,a)') "Input: '",gmd_opt ,"' is not a valid GuessMagnetDensity"
          message(2) = '(GuessMagnetDensity = 1 | 2 | 3 | 4)'
          call write_fatal(2)
       end if
    end if

    rho = M_ZERO
    select case (gmd_opt)
    case (1) ! Paramagnetic
       allocate(atom_rho(m%np, 1))
       do ia = 1, geo%natoms
          call atom_density(m, sb, geo%atom(ia), 1, atom_rho(1:m%np, 1:1))
          rho(1:m%np, 1:1) = rho(1:m%np, 1:1) + atom_rho(1:m%np, 1:1)
       end do
       if (spin_channels == 2) then
          rho(1:m%np, 1) = M_HALF*rho(1:m%np, 1)
          rho(1:m%np, 2) = rho(1:m%np, 1)
       end if

    case (2) ! Ferromagnetic
       allocate(atom_rho(m%np, 2))
       do ia = 1, geo%natoms
          call atom_density(m, sb, geo%atom(ia), 2, atom_rho(1:m%np, 1:2))
          call lalg_axpy(m%np, 2, M_ONE, atom_rho, rho)
       end do

    case (3) ! Random oriented spins
       allocate(atom_rho(m%np, 2))
       do ia = 1, geo%natoms
          call atom_density(m, sb, geo%atom(ia), 2, atom_rho)

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

    case (4) ! User-defined

       if(loct_parse_block(check_inp('AtomsMagnetDirection'), blk) < 0) then
          message(1) = "AtomsMagnetDirection block is not defined "
          call write_fatal(1)
       end if

       if (loct_parse_block_n(blk) /= geo%natoms) then
          message(1) = "AtomsMagnetDirection block has the wrong number of rows"
          call write_fatal(1)
       end if

       allocate(atom_rho(m%np, 2))
       do ia = 1, geo%natoms
          call atom_density(m, sb, geo%atom(ia), 2, atom_rho)

          if (nspin == 2) then
             call loct_parse_block_float(blk, ia-1, 0, mag(1))
             if (mag(1) > M_ZERO) then
                rho(1:m%np, 1:2) = rho(1:m%np, 1:2) + atom_rho(1:m%np, 1:2)
             else
                rho(1:m%np, 1) = rho(1:m%np, 1) + atom_rho(1:m%np, 2)
                rho(1:m%np, 2) = rho(1:m%np, 2) + atom_rho(1:m%np, 1)
             end if

          elseif (nspin == 4) then
             do i = 1, 3
                call loct_parse_block_float(blk, ia-1, i-1, mag(i))
                if (abs(mag(i)) < CNST(1.0e-20)) mag(i) = M_ZERO
             end do

             theta = acos(mag(3)/sqrt(dot_product(mag, mag)))
             if (mag(1) == M_ZERO) then
                if (mag(2) == M_ZERO) then
                   phi = M_ZERO
                elseif (mag(2) < M_ZERO) then
                   phi = M_PI*M_TWOTHIRD
                elseif (mag(2) > M_ZERO) then
                   phi = M_PI*M_HALF
                end if
             else
                if (mag(2) < M_ZERO) then
                   phi = M_TWO*M_PI - acos(mag(1)/sin(theta)/sqrt(dot_product(mag, mag)))
                elseif (mag(2) >= M_ZERO) then
                   phi = acos(mag(1)/sin(theta)/sqrt(dot_product(mag, mag)))
                end if
             end if

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

       call loct_parse_block_end(blk)

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

    deallocate(atom_rho)
    call pop_sub()
  end subroutine system_guess_density


  !----------------------------------------------------------
  subroutine dsystem_h_setup(sys, h)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h

    call push_sub('hamiltonian_setup')

    call states_fermi(sys%st, sys%gr%m)
    call dstates_calc_dens(sys%st, sys%gr%m%np, sys%st%rho)

    call dv_ks_calc(sys%gr, sys%ks, h, sys%st, calc_eigenval=.true.) ! get potentials
    call states_fermi(sys%st, sys%gr%m)                            ! occupations
    call hamiltonian_energy(h, sys%st, sys%gr%geo%eii, -1)            ! total energy

    call pop_sub()
  end subroutine dsystem_h_setup


  !----------------------------------------------------------
  subroutine zsystem_h_setup(sys, h)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h

    call push_sub('hamiltonian_setup')

    call states_fermi(sys%st, sys%gr%m)
    call zstates_calc_dens(sys%st, sys%gr%m%np, sys%st%rho)

    call zv_ks_calc(sys%gr, sys%ks, h, sys%st, calc_eigenval=.true.) ! get potentials
    call states_fermi(sys%st, sys%gr%m)                            ! occupations
    call hamiltonian_energy(h, sys%st, sys%gr%geo%eii, -1)            ! total energy

    call pop_sub()
  end subroutine zsystem_h_setup

end module system
