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

module hamiltonian_m
  use global_m
  use messages_m
  use datasets_m
  use profiling_m
  use units_m
  use lib_oct_parser_m
  use lib_basic_alg_m
  use functions_m
  use mesh_m
  use grid_m
  use simul_box_m
  use mesh_function_m
  use geometry_m
  use specie_m
  use states_m
  use external_pot_m
  use output_m
  use mpi_m
  use varinfo_m

  implicit none

  private
  public ::                &
    hamiltonian_t,         &
    hamiltonian_init,      &
    hamiltonian_end,       &
    hamiltonian_energy,    &
    hamiltonian_output,    &
    hamiltonian_span,      &
    dhamiltonian_eigenval, &
    zhamiltonian_eigenval, &
    dhpsi,                 &
    dvlpsi,                &
    dvnlpsi,               &
    dmagnus,               &
    dkinetic,              &
    zhpsi,                 &
    zvlpsi,                &
    zvnlpsi,               &
    zmagnus,               &
    zkinetic

#ifdef COMPLEX_WFNS
  public :: zso
#endif


  type hamiltonian_t
    ! The Hamiltonian must know what are the "dimensions" of the spaces,
    ! in order to be able to operate on the states.
    type(states_dim_t), pointer :: d

    integer :: reltype            ! type of relativistic correction to use

    FLOAT, pointer :: vhartree(:)
    FLOAT, pointer :: Vhxc(:,:)   ! xc potential + hartree potential
    FLOAT, pointer :: axc(:,:,:)  ! xc vector-potential divided by c

    ! the energies (total, ion-ion, exchange, correlation)
    FLOAT :: etot, eii, ex, ec, exc_j, epot

    logical :: ip_app             ! independent particle approximation, or not.
    ! copied from sys%ks

    type(epot_t) :: ep         ! handles the external potential

    ! gauge
    integer :: gauge              ! in which gauge shall we work in

    ! absorbing boundaries
    integer  :: ab                ! do we have absorbing boundaries?
    FLOAT :: ab_width             ! width of the absorbing boundary
    FLOAT :: ab_height            ! height of the absorbing boundary
    FLOAT, pointer :: ab_pot(:)   ! where we store the ab potential

    ! Spectral range
    FLOAT :: spectral_middle_point
    FLOAT :: spectral_half_span

    ! Kinetic Cutoff
    FLOAT :: cutoff

  end type hamiltonian_t

  integer, public, parameter :: &
    NOREL      = 0,             &
    SPIN_ORBIT = 1,             &
    APP_ZORA   = 2,             &
    ZORA       = 3

  integer, public, parameter :: &
    LENGTH     = 1,             &
    VELOCITY   = 2

  integer, public, parameter :: &
    NO_ABSORBING        = 0,    &
    IMAGINARY_ABSORBING = 1,    &
    MASK_ABSORBING      = 2


contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_init(h, gr, states_dim, ip_app)
    type(hamiltonian_t), intent(out)   :: h
    type(grid_t),        intent(inout) :: gr
    type(states_dim_t),  pointer       :: states_dim
    logical,                intent(in)    :: ip_app

    integer :: i, j, n
    FLOAT :: d(3)

    call push_sub('h.hamiltonian_init')

    ! make a couple of local copies
    h%ip_app = ip_app
    h%d => states_dim

    ! initialize variables
    h%epot = M_ZERO
    h%ex = M_ZERO; h%ec = M_ZERO
    h%etot = M_ZERO

    ! allocate potentials and density of the cores
    ! In the case of spinors, vxc_11 = h%vxc(:, 1), vxc_22 = h%vxc(:, 2), Re(vxc_12) = h%vxc(:. 3);
    ! Im(vxc_12) = h%vxc(:, 4)
    ALLOCATE(h%vhartree(NP),        NP)
    ALLOCATE(h%Vhxc(NP, h%d%nspin), NP*h%d%nspin)
    h%Vhxc = M_ZERO
    if (h%d%cdft) then
      ALLOCATE(h%axc(NP, NDIM, h%d%nspin), NP*NDIM*h%d%nspin)
      h%axc = M_ZERO
    else
      nullify(h%axc)
    end if

    call epot_init(h%ep, gr)

    !%Variable RelativisticCorrection
    !%Type integer
    !%Default non_relativistic
    !%Section Hamiltonian
    !%Description
    !% The default value means that <i>no</i> relativistic correction is used. To
    !% include spin-orbit coupling turn <tt>RelativisticCorrection</tt> to <tt>spin_orbit</tt> 
    !% (this will only work when using an executable compiled for complex wave-functions,
    !% and if <tt>SpinComponents</tt> has been set to <tt>non_collinear</tt>, which ensures
    !% the use of spinors).
    !%Option non_relativistic 0
    !% No relativistic corrections.
    !%Option spin_orbit 1
    !% Spin-Orbit.
    !%Option app_zora 2
    !% Approximated ZORA (Not implemented)
    !%Option zora 3
    !% ZORA (Not implemented)
    !%End
    call loct_parse_int(check_inp('RelativisticCorrection'), NOREL, h%reltype)
    if(.not.varinfo_valid_option('RelativisticCorrection', h%reltype)) call input_error('RelativisticCorrection')

#ifdef COMPLEX_WFNS
    call messages_print_var_option(stdout, "RelativisticCorrection", h%reltype)

    ! This is temporary...
    if(h%reltype > SPIN_ORBIT) then
      message(1) = 'Error: ZORA corrections not working yet. Visit us soon.'
      call write_fatal(1)
    end if
#else
    if(h%reltype .ne. NOREL) then
      message(1) = "Cannot apply relativistic corrections with an executable compiled"
      message(2) = "for real wavefunctions."
      call write_fatal(2)
    end if
#endif

    !%Variable TDGauge
    !%Type integer
    !%Default length
    !%Section Time Dependent
    !%Description
    !% In which gauge to treat the laser field.
    !%Option length 1
    !% Length gauge.
    !%Option velocity 2
    !% Velocity gauge.
    !%End
    call loct_parse_int(check_inp('TDGauge'), LENGTH, h%gauge)
    if(.not.varinfo_valid_option('TDGauge', h%gauge)) call input_error('TDGauge')
    call messages_print_var_option(stdout, "TDGauge", h%gauge)

    !%Variable AbsorbingBoundaries
    !%Type integer
    !%Default no_absorbing
    !%Section Time Dependent::Absorbing Boundaries
    !%Description
    !% To improve the quality of the spectra by avoiding the formation of 
    !% standing density waves, one can make the boundaries of the simulation 
    !% box absorbing.
    !%Option no_absorbing 0
    !% No absorbing boundaries.
    !%Option sin2 1
    !% A <math>\sin^2</math> imaginary potential is added at the boundaries.
    !%Option mask 2
    !% A mask is applied to the wave-functions at the boundaries.
    !%End
    call loct_parse_int(check_inp('AbsorbingBoundaries'), NO_ABSORBING, h%ab)
    if(.not.varinfo_valid_option('AbsorbingBoundaries', h%ab)) call input_error('AbsorbingBoundaries')
    call messages_print_var_option(stdout, "AbsorbingBoundaries", h%ab)
    
    nullify(h%ab_pot)

    absorbing_boundaries: if(h%ab.ne.NO_ABSORBING) then
      !%Variable ABWidth
      !%Type float
      !%Default 0.4 a.u.
      !%Section Time Dependent::Absorbing Boundaries
      !%Description
      !% Width of the region used to apply the absorbing boundaries.
      !%End
      call loct_parse_float(check_inp('ABWidth'), CNST(0.4)/units_inp%length%factor, h%ab_width)
      h%ab_width  = h%ab_width * units_inp%length%factor
      if(h%ab == 1) then
        !%Variable ABHeight
        !%Type float
        !%Default -0.2 a.u.
        !%Section Time Dependent::Absorbing Boundaries
        !%Description 
        !% When <tt>AbsorbingBoundaries == sin2</tt>, is the height of the imaginary potential.
        !%End
        call loct_parse_float(check_inp('ABHeight'), -CNST(0.2)/units_inp%energy%factor, h%ab_height)
        h%ab_height = h%ab_height * units_inp%energy%factor
      else
        h%ab_height = M_ONE
      end if

      ! generate boundary potential...
      ALLOCATE(h%ab_pot(NP), NP)
      h%ab_pot = M_ZERO
      do i = 1, NP
        call mesh_inborder(gr%m, i, n, d, h%ab_width)
        if(n>0) then
          do j = 1, n
            h%ab_pot(i) = h%ab_pot(i) + h%ab_height * sin(d(j)*M_PI/(M_TWO*h%ab_width))**2
          end do
        end if
        if(abs(h%ab_pot(i)) > abs(h%ab_height)) h%ab_pot(i) = h%ab_height
      end do

    end if absorbing_boundaries

    ! Cutoff applied to the kinetic term.
    ! it is used *both* in the calculation of the derivatives and in the split operator method
    call loct_parse_float(check_inp('KineticCutoff'), -M_ONE, h%cutoff)
    if(h%cutoff > M_ZERO) then
      h%cutoff = h%cutoff * units_inp%energy%factor
      write(message(1),'(a,f7.2,a)') 'Info: The kinetic operator will have a cutoff of',&
        h%cutoff/units_out%energy%factor, units_out%energy%abbrev
      write(message(2),'(a)')        '      (only if DerivativesSpace = 1 is set)'
      call write_info(2)
    end if

    call pop_sub()
  end subroutine hamiltonian_init


  ! ---------------------------------------------------------
  subroutine hamiltonian_end(h, gr)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),        intent(in)    :: gr

    call push_sub('h.hamiltonian_end')

    if(associated(h%vhartree)) then
      deallocate(h%vhartree)
      nullify(h%vhartree)
    end if
    if(associated(h%vhxc)) then
      deallocate(h%vhxc)
      nullify(h%vhxc)
    end if
    if(associated(h%axc)) then
      deallocate(h%axc)
      nullify(h%axc)
    end if

    call epot_end(h%ep, gr%sb, gr%geo)

    if(associated(h%ab_pot)) then
      deallocate(h%ab_pot); nullify(h%ab_pot)
    end if

    if(associated(h%d)) then
      nullify(h%d)
    end if

    call pop_sub()
  end subroutine hamiltonian_end


  ! ---------------------------------------------------------
  ! This subroutine calculates the total energy of the system. Basically, it
  ! adds up the KS eigenvalues, and then it substracts the whatever double
  ! counts exist (see TDDFT theory for details).
  subroutine hamiltonian_energy(h, st, eii, iunit)
    type(hamiltonian_t), intent(inout) :: h
    type(states_t),      intent(in)    :: st
    FLOAT,                  intent(in)    :: eii
    integer,                intent(in)    :: iunit

    FLOAT :: e
#ifdef HAVE_MPI
    FLOAT :: s
    integer :: err
#endif

    call push_sub('h.hamiltonian_energy')

    e = states_eigenvalues_sum(st)

#ifdef HAVE_MPI
    if(st%parallel_in_states) then
      call MPI_ALLREDUCE(e, s, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, err)
      e = s
    end if
#endif

    h%eii    = eii
    h%etot   = e + eii + h%epot + h%ex + h%ec

    if (iunit > 0) then
      write(message(1), '(6x,a, f15.8)')'Ion-ion     = ', h%eii  / units_out%energy%factor
      write(message(2), '(6x,a, f15.8)')'Eigenvalues = ', e      / units_out%energy%factor
      write(message(3), '(6x,a, f15.8)')'Potentials  = ', h%epot / units_out%energy%factor
      write(message(4), '(6x,a, f15.8)')'Exchange    = ', h%ex   / units_out%energy%factor
      write(message(5), '(6x,a, f15.8)')'Correlation = ', h%ec   / units_out%energy%factor
      write(message(6), '(6x,a, f15.8)')'Total       = ', h%etot / units_out%energy%factor
      call write_info(6, iunit)
    end if

    call pop_sub()
  end subroutine hamiltonian_energy


  ! ---------------------------------------------------------
  subroutine hamiltonian_span(h, delta, emin)
    type(hamiltonian_t), intent(inout) :: h
    FLOAT,                  intent(in)    :: delta, emin

    call push_sub('h.hamiltonian_span')

    h%spectral_middle_point = ((M_Pi**2/(2*delta**2)) + emin)/M_TWO
    h%spectral_half_span    = ((M_Pi**2/(2*delta**2)) - emin)/M_TWO

    call pop_sub()
  end subroutine hamiltonian_span


  ! ---------------------------------------------------------
  subroutine hamiltonian_output(h, m, sb, dir, outp)
    type(hamiltonian_t), intent(in) :: h
    type(mesh_t),        intent(in) :: m
    type(simul_box_t),   intent(in) :: sb
    character(len=*),       intent(in) :: dir
    type(output_t),      intent(in) :: outp

    integer :: is, err
    character(len=80) :: fname
    FLOAT :: u

    call push_sub('h.hamiltonian_output')

    u = units_out%energy%factor
    if(iand(outp%what, output_potential).ne.0) then
      call doutput_function(outp%how, dir, "v0", m, sb, h%ep%vpsl, u, err)

      if(h%ep%classic_pot > 0) then
        call doutput_function(outp%how, dir, "vc", m, sb, h%ep%Vclassic, u, err)
      end if

      if(.not.h%ip_app) then
        do is = 1, min(h%d%ispin, 2)
          write(fname, '(a,i1)') 'vhxc-', is
          call doutput_function(outp%how, dir, fname, m, sb, h%Vhxc(:, is), u, err)
        end do
      end if
    end if

    call pop_sub()
  end subroutine hamiltonian_output

#include "undef.F90"
#include "real.F90"
#include "h_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "h_inc.F90"

#if defined(COMPLEX_WFNS)
#include "h_so.F90"
#endif

end module hamiltonian_m
