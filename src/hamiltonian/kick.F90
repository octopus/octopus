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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module kick_oct_m
  use iso_c_binding
  use global_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use loct_math_oct_m
  use math_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pcm_eom_oct_m
  use pcm_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use species_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use symm_op_oct_m
  use symmetries_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::               &
    kick_t,               &
    kick_init,            &
    kick_copy,            &
    kick_end,             &
    kick_read,            &
    kick_write,           &
    kick_apply,           &
    kick_function_get,    &
    kick_get_type


  integer, public, parameter ::        &
    KICK_FUNCTION_DIPOLE        = 0,   &
    KICK_FUNCTION_MULTIPOLE     = 1,   &
    KICK_FUNCTION_USER_DEFINED  = 2

  integer, public, parameter ::    &
    KICK_DENSITY_MODE        = 0,  &
    KICK_SPIN_MODE           = 1,  &
    KICK_SPIN_DENSITY_MODE   = 2,  &
    KICK_MAGNON_MODE         = 3

  integer, public, parameter ::    &
    QKICKMODE_NONE           = 0,  &
    QKICKMODE_EXP            = 1,  &
    QKICKMODE_COS            = 2,  &
    QKICKMODE_SIN            = 3,  &
    QKICKMODE_BESSEL         = 4


  type kick_t
    ! Components are public by default

    !> Dimensions
    integer           :: dim
    !> The time which the kick is applied (normally, this is zero)
    FLOAT             :: time
    !> The strength, and strength "mode".
    integer, private  :: delta_strength_mode
    FLOAT             :: delta_strength
    !> In case we use a normal dipole kick:
    FLOAT             :: pol(MAX_DIM, MAX_DIM)
    integer           :: pol_dir
    integer           :: pol_equiv_axes
    FLOAT             :: wprime(MAX_DIM)
    FLOAT             :: easy_axis(MAX_DIM)
    !> In case we have a general multipolar kick,
    !! the form of this "kick" will be (atomic units):
    !! \f[
    !! V(\vec{r}) = sum_{i=1}^{n\_multipoles}
    !!                weight(i) * (e^2 / a_0^(l+1)) * r^l(i) * Y_{l(i),m(i)} (\vec{r})
    !! \f]
    !! which has units of energy; if we include the time-dependence (delta function):
    !! \f[
    !! V(\vec{r}) = sum_{i=1}^{n\_multipoles}
    !!                 weight(i) * (\hbar / a_0^l) * r^l(i) * Y_{l(i),m(i)} (\vec{r}) * \delta(t)
    !! \f]
    integer              :: n_multipoles
    integer, allocatable :: l(:), m(:)
    FLOAT,   allocatable :: weight(:)
    integer              :: nqmult(1:MAX_DIM)
    integer              :: nqvec
    FLOAT,   allocatable :: qvector(:,:)
    FLOAT                :: trans_vec(MAX_DIM,2)
    FLOAT                :: qlength
    integer              :: qkick_mode
    integer              :: qbessel_l, qbessel_m
    !> In case we use a general function
    integer              :: function_mode
    character(len=200), private:: user_defined_function
  end type kick_t

contains

  ! ---------------------------------------------------------
  subroutine kick_init(kick, namespace, space, kpoints, nspin)
    type(kick_t),      intent(out) :: kick
    type(namespace_t), intent(in)  :: namespace
    type(space_t),     intent(in)  :: space
    type(kpoints_t),   intent(in)  :: kpoints
    integer,           intent(in)  :: nspin

    type(block_t) :: blk
    integer :: n_rows, irow, idir, iop, iq, iqx, iqy, iqz
    FLOAT :: norm, dot
    FLOAT :: qtemp(1:MAX_DIM)
    integer :: periodic_dim

    PUSH_SUB(kick_init)

    kick%dim = space%dim
    periodic_dim = space%periodic_dim

    !%Variable TDDeltaKickTime
    !%Type float
    !%Default 0.0
    !%Section Time-Dependent::Response
    !%Description
    !% The delta-perturbation that can be applied by making use of the <tt>TDDeltaStrength</tt> variable,
    !% can be applied at the time given by this variable. Usually, this time is zero, since one wants
    !% to apply the delta pertubation or "kick" at a system at equilibrium, and no other time-dependent
    !% external potential is used. However, one may want to apply a kick on top of a laser field,
    !% for example.
    !%End
    call parse_variable(namespace, 'TDDeltaKickTime', M_ZERO, kick%time, units_inp%time)
    if(kick%time > M_ZERO) then
      call messages_experimental('TDDeltaKickTime > 0')
    end if

    !%Variable TDDeltaStrength
    !%Type float
    !%Default 0
    !%Section Time-Dependent::Response
    !%Description
    !% When no laser is applied, a delta (in time) perturbation with
    !% strength <tt>TDDeltaStrength</tt> can be applied. This is used to
    !% calculate, <i>e.g.</i>, the linear optical spectra. If the ions are
    !% allowed to move, the kick will affect them also.
    !% The electric field is <math>-(\hbar k / e) \delta(t)</math> for a dipole with
    !% zero wavevector, where <i>k</i> = <tt>TDDeltaStrength</tt>, which causes
    !% the wavefunctions instantaneously to acquire a phase <math>e^{ikx}</math>.
    !% The unit is inverse length.
    !%End
    call parse_variable(namespace, 'TDDeltaStrength', M_ZERO, kick%delta_strength, units_inp%length**(-1))

    kick%function_mode = KICK_FUNCTION_DIPOLE

    if(abs(kick%delta_strength) <= M_EPSILON) then
      kick%delta_strength_mode = 0
      kick%pol_equiv_axes = 0
      kick%pol = M_ZERO
      do idir = 1, kick%dim
        kick%pol(idir, idir) = M_ONE
      end do
      kick%pol_dir = 0
      kick%wprime = M_ZERO
      kick%n_multipoles = 0
      kick%qkick_mode = QKICKMODE_NONE
      kick%easy_axis(1:MAX_DIM) = M_ZERO
      POP_SUB(kick_init)
      return
    end if

    !%Variable TDDeltaStrengthMode
    !%Type integer
    !%Default kick_density
    !%Section Time-Dependent::Response
    !%Description
    !% When calculating the density response via real-time propagation,
    !% one needs to perform an initial kick on the KS system, at
    !% time zero. Depending on what kind of response property one wants to obtain,
    !% this kick may be done in several modes. For use to calculate triplet excitations,
    !% see MJT Oliveira, A Castro, MAL Marques, and A Rubio, <i>J. Nanoscience and Nanotechnology</i> <b>8</b>, 3392 (2008).
    !%Option kick_density 0
    !% The total density of the system is perturbed. This mode is appropriate for
    !% electric dipole response, as for optical absorption.
    !%Option kick_spin 1
    !% The individual spin densities are perturbed oppositely. Note that this mode
    !% is only possible if the run is done in spin-polarized mode, or with spinors.
    !% This mode is appropriate for the paramagnetic dipole response, which can couple
    !% to triplet excitations.
    !%Option kick_spin_and_density 2
    !% A combination of the two above. Note that this mode
    !% is only possible if the run is done in spin-polarized mode, or with spinors.
    !% This mode is intended for use with symmetries to obtain both of the responses
    !% at once, at described in the reference above.
    !%Option kick_magnon 3
    !% Rotates the magnetization. Only works for spinors.
    !% Can be used in a supercell or my making use of the generalized Bloch theorem.
    !% In the later case (see <tt>SpiralBoundaryConditions</tt>) spin-orbit coupling cannot be used.
    !%End
    call parse_variable(namespace, 'TDDeltaStrengthMode', KICK_DENSITY_MODE, kick%delta_strength_mode)
    select case (kick%delta_strength_mode)
    case (KICK_DENSITY_MODE)
    case (KICK_SPIN_MODE, KICK_SPIN_DENSITY_MODE)
    case (KICK_MAGNON_MODE)
      if(nspin /= SPINORS) call messages_input_error(namespace, 'TDDeltaStrengthMode', 'Magnon kick is incompatible with spinors')
    case default
      call messages_input_error(namespace, 'TDDeltaStrengthMode', 'Unknown mode')
    end select

    if(parse_is_defined(namespace, 'TDDeltaUserDefined')) then

      kick%function_mode = KICK_FUNCTION_USER_DEFINED
      kick%n_multipoles = 0

      !%Variable TDDeltaUserDefined
      !%Type string
      !%Section Time-Dependent::Response
      !%Description
      !% By default, the kick function will be a dipole. This will change if (1) the variable
      !% <tt>TDDeltaUserDefined</tt> is present in the inp file, or (2) if the block <tt>TDKickFunction</tt>
      !% is present in the <tt>inp</tt> file. If both are present in the <tt>inp</tt> file, the <tt>TDKickFunction</tt>
      !% block will be ignored. The value of <tt>TDDeltaUserDefined</tt> should be a string describing
      !% the function that is going to be used as delta perturbation.
      !%End
      call parse_variable(namespace, 'TDDeltaUserDefined', '0', kick%user_defined_function)

      !%Variable TDKickFunction
      !%Type block
      !%Section Time-Dependent::Response
      !%Description
      !% If the block <tt>TDKickFunction</tt> is present in the input file, and the variable
      !% <tt>TDDeltaUserDefined</tt> is not present in the input file, the kick function to
      !% be applied at time zero of the time-propagation will not be a "dipole" function
      !% (<i>i.e.</i> <math>\phi \rightarrow e^{ikx} \phi</math>, but a general multipole in the form <math>r^l Y_{lm}(r)</math>.
      !%
      !% Each line has three columns: integers <i>l</i> and <i>m</i> that defines the
      !% multipole, and a weight. Any number of lines may be given, and the kick will be the sum of those
      !% multipoles with the given weights.
      !%
      !% This feature allows calculation of quadrupole, octupole, etc., response functions.
      !%End
    else if(parse_block(namespace, 'TDKickFunction', blk) == 0) then

      kick%function_mode = KICK_FUNCTION_MULTIPOLE
      n_rows = parse_block_n(blk)
      kick%n_multipoles = n_rows
      SAFE_ALLOCATE(     kick%l(1:n_rows))
      SAFE_ALLOCATE(     kick%m(1:n_rows))
      SAFE_ALLOCATE(kick%weight(1:n_rows))
      do irow = 1, n_rows
        call parse_block_integer(blk, irow - 1, 0, kick%l(irow))
        call parse_block_integer(blk, irow - 1, 1, kick%m(irow))
        call parse_block_float(blk, irow - 1, 2, kick%weight(irow))
        if( (kick%l(irow) < 0) .or. (abs(kick%m(irow)) > abs(kick%l(irow))) ) then
          call messages_input_error(namespace, 'TDkickFunction')
        end if
      end do

    else

      kick%function_mode = KICK_FUNCTION_DIPOLE
      kick%n_multipoles = 0

      ! Find out how many equivalent axes we have...
      !%Variable TDPolarizationEquivAxes
      !%Type integer
      !%Default 0
      !%Section Time-Dependent::Response::Dipole
      !%Description
      !% Defines how many of the <tt>TDPolarization</tt> axes are equivalent. This information is stored in a file and then
      !% used by <tt>oct-propagation_spectrum</tt> to rebuild the full polarizability tensor from just the
      !% first <tt>TDPolarizationEquivAxes</tt> directions. This variable is also used by <tt>CalculationMode = vdw</tt>.
      !%End
      call parse_variable(namespace, 'TDPolarizationEquivAxes', 0, kick%pol_equiv_axes)

      !%Variable TDPolarizationDirection
      !%Type integer
      !%Section Time-Dependent::Response::Dipole
      !%Description
      !% When a delta potential is included in a time-dependent run, this
      !% variable defines in which direction the field will be applied
      !% by selecting one of the lines of <tt>TDPolarization</tt>. In a
      !% typical run (without using symmetry), the <tt>TDPolarization</tt> block
      !% would contain the three Cartesian unit vectors (the default
      !% value), and one would make 3 runs varying
      !% <tt>TDPolarization</tt> from 1 to 3.
      !% If one is using symmetry,  <tt>TDPolarization</tt> should run only from 1
      !% to <tt>TDPolarizationEquivAxes</tt>.
      !%End

      call parse_variable(namespace, 'TDPolarizationDirection', 0, kick%pol_dir)

      if(kick%delta_strength_mode /= KICK_MAGNON_MODE) then
        if(kick%pol_dir < 1 .or. kick%pol_dir > kick%dim) call messages_input_error(namespace, 'TDPolarizationDirection')
      end if

      !%Variable TDPolarization
      !%Type block
      !%Section Time-Dependent::Response::Dipole
      !%Description
      !% The (real) polarization of the delta electric field. Normally
      !% one needs three perpendicular polarization directions to calculate a
      !% spectrum (unless symmetry is used).
      !% The format of the block is:
      !%
      !% <tt>%TDPolarization
      !% <br>&nbsp;&nbsp;pol1x | pol1y | pol1z
      !% <br>&nbsp;&nbsp;pol2x | pol2y | pol2z
      !% <br>&nbsp;&nbsp;pol3x | pol3y | pol3z
      !% <br>%</tt>
      !%
      !% <tt>Octopus</tt> uses both this block and the variable
      !% <tt>TDPolarizationDirection</tt> to determine the polarization
      !% vector for the run. For example, if
      !% <tt>TDPolarizationDirection=2</tt> the polarization <tt>(pol2x,
      !% pol2y, pol2z)</tt> would be used.
      !% These directions may not be in periodic directions.
      !%
      !% The default value for <tt>TDPolarization</tt> is the three
      !% Cartesian unit vectors (1,0,0), (0,1,0), and (0,0,1).
      !%
      !% Note that the directions do not necessarily need to be perpendicular
      !% when symmetries are used.
      !%
      !% WARNING: If you want to obtain the cross-section tensor, the
      !% <tt>TDPolarization</tt> block must be exactly the same for the run in
      !% each direction. The direction must be selected by the
      !% <tt>TDPolarizationDirection</tt> variable.
      !%
      !%End

      ! Default basis is the Cartesian unit vectors.
      ! FIXME: Here the symmetry of the system should be analyzed, and the polarization
      ! basis built accordingly.
      kick%pol(:, :) = M_ZERO
      do idir = 1, kick%dim
        kick%pol(idir, idir) = M_ONE
      end do
      if(parse_block(namespace, 'TDPolarization', blk)==0) then
        n_rows = parse_block_n(blk)

        if(n_rows < kick%dim) call messages_input_error(namespace, 'TDPolarization', 'There should be one line per dimension')
        
        do irow = 1, n_rows
          do idir = 1, kick%dim
            call parse_block_float(blk, irow - 1, idir - 1, kick%pol(idir, irow))
          end do
        end do
        call parse_block_end(blk)
      end if

      ! Normalize
      do idir = 1, kick%dim
        kick%pol(1:kick%dim, idir) = kick%pol(1:kick%dim, idir) / sqrt(sum(kick%pol(1:kick%dim, idir)**2))
      end do

      if(kick%delta_strength_mode /= KICK_MAGNON_MODE) then
        if(any(abs(kick%pol(1:periodic_dim, :)) > M_EPSILON)) then
          message(1) = "Kick cannot be applied in a periodic direction. Use GaugeVectorField instead."
          call messages_fatal(1, namespace=namespace)
        end if
      end if

      !%Variable TDPolarizationWprime
      !%Type block
      !%Section Time-Dependent::Response::Dipole
      !%Description
      !% This block is needed only when
      !% <tt>TDPolarizationEquivAxes</tt> is set to 3.  In such a case,
      !% the three directions (<i>pol1</i>, <i>pol2</i>, and <i>pol3</i>) defined in
      !% the <tt>TDPolarization</tt> block should be related by symmetry
      !% operations. If <i>A</i> is the symmetry operation that takes you
      !% from <i>pol1</i> to <i>pol2</i>, then <tt>TDPolarizationWprime</tt>
      !% should be set to the direction defined by <i>A</i><math>^{-1}</math><i>pol3</i>.
      !% For more information see MJT Oliveira
      !% <i>et al.</i>, <i>J. Nanoscience and Nanotechnology</i> <b>8</b>,
      !% 3392 (2008).
      !%End
      if(parse_block(namespace, 'TDPolarizationWprime', blk)==0) then
        do idir = 1, kick%dim
          call parse_block_float(blk, 0, idir - 1, kick%wprime(idir))
        end do
        kick%wprime(1:kick%dim) = kick%wprime(1:kick%dim) / sqrt(sum(kick%wprime(1:kick%dim)**2))
        call parse_block_end(blk)
      else
        kick%wprime(1:kick%dim-1) = M_ZERO
        kick%wprime(kick%dim) = M_ONE
      end if
    end if

    ! for non-dipole, it is more complicated to check whether it is actually in the periodic direction
    if(periodic_dim > 0 .and. kick%delta_strength_mode /= KICK_MAGNON_MODE) then
      message(1) = "Kicks cannot be applied correctly in periodic directions."
      call messages_warning(1, namespace=namespace)
    end if

    !%Variable TDMomentumTransfer
    !%Type block
    !%Section Time-Dependent::Response
    !%Description
    !% Momentum-transfer vector for the calculation of the dynamic structure factor.
    !% When this variable is set, a non-dipole field is applied, and an output file
    !% <tt>ftchd</tt> is created (it contains the Fourier transform of the charge density
    !% at each time). The type of the applied external field can be set by
    !% an optional last number. Possible options are <tt>qexp</tt> (default), <tt>qcos</tt>,
    !% <tt>qsin</tt>, or <tt>qcos+qsin</tt>. In the formulae below,
    !% <math>\vec{q}</math> is the momentum-transfer vector.
    !%Option qexp 1
    !% External field is <math>e^{i \vec{q} \cdot \vec{r}}</math>.
    !%Option qcos 2
    !% External field is <math>\cos \left( i \vec{q} \cdot \vec{r} \right)</math>.
    !%Option qsin 3
    !% External field is <math>\sin \left( i \vec{q} \cdot \vec{r} \right)</math>.
    !%Option qbessel 4
    !% External field is <math>j_l \left( \vec{q} \cdot \vec{r} \right) Y_{lm} \left(\vec{r} \right)</math>.
    !% In this case, the block has to include two extra values (<i>l</i> and <i>m</i>).
    !%End

    if(parse_block(namespace, 'TDMomentumTransfer', blk)==0) then
      kick%nqvec = 1
      SAFE_ALLOCATE(kick%qvector(1:MAX_DIM,1))
      do idir = 1, MAX_DIM
        call parse_block_float(blk, 0, idir - 1, kick%qvector(idir,1))
        kick%qvector(idir,1) = units_to_atomic(unit_one / units_inp%length, kick%qvector(idir,1))
      end do

      ! Read the calculation mode (exp, cos, sin, or bessel)
      if(parse_block_cols(blk, 0) > MAX_DIM) then

        call parse_block_integer(blk, 0, idir - 1, kick%qkick_mode)

        ! Read l and m if bessel mode (j_l*Y_lm) is used
        if(kick%qkick_mode == QKICKMODE_BESSEL .and. parse_block_cols(blk, 0) == MAX_DIM+3) then
          call parse_block_integer(blk, 0, idir + 0, kick%qbessel_l)
          call parse_block_integer(blk, 0, idir + 1, kick%qbessel_m)
        else
          kick%qbessel_l = 0
          kick%qbessel_m = 0
        end if

      else
        kick%qkick_mode = QKICKMODE_EXP
      end if

      call parse_block_end(blk)

      if(kpoints%use_symmetries) then
        do iop = 1, symmetries_number(kpoints%symm)
          if(iop == symmetries_identity_index(kpoints%symm)) cycle
          if(.not. symm_op_invariant_cart(kpoints%symm%ops(iop), kick%qvector(:,1), CNST(1e-5))) then
            message(1) = "The TDMomentumTransfer breaks (at least) one of the symmetries used to reduce the k-points."
            message(2) = "Set SymmetryBreakDir equal to TDMomemtumTransfer."
            call messages_fatal(2, namespace=namespace)
          end if
        end do
      end if

    else
      kick%qkick_mode = QKICKMODE_NONE
      kick%nqvec = 1
      SAFE_ALLOCATE(kick%qvector(1:MAX_DIM,1))
      kick%qvector(:,1) = M_ZERO
    end if

    kick%qlength = sqrt(sum(kick%qvector(:,1)**2))

    if(kick%delta_strength_mode == KICK_MAGNON_MODE) then
      !%Variable TDEasyAxis
      !%Type block
      !%Section Time-Dependent::Response::Dipole
      !%Description
      !% For magnon kicks only.
      !% This variable defines the direction of the easy axis of the crystal.
      !% The magnetization is kicked in the plane transverse to this vector
      !%End
      if(parse_block(namespace, 'TDEasyAxis', blk)==0) then
        n_rows = parse_block_n(blk)

        do idir = 1, 3
          call parse_block_float(blk, 0, idir - 1, kick%easy_axis(idir))
        end do
        norm = sqrt(sum(kick%easy_axis(1:3)**2))
        if(norm < CNST(1e-9)) then
          message(1) = "TDEasyAxis norm is too small."
          call messages_fatal(1, namespace=namespace)
        end if
        kick%easy_axis(1:3) = kick%easy_axis(1:3)/norm
        call parse_block_end(blk)
      else
        message(1) = "For magnons, the variable TDEasyAxis must be defined."
        call messages_fatal(1, namespace=namespace)
      end if

      !We first two vectors defining a basis in the transverse plane
      !For this we take two vectors not align with the first one
      !and we perform a Gram-Schmidt orthogonalization
      kick%trans_vec(1,1) = -kick%easy_axis(2)
      kick%trans_vec(2,1) = M_TWO*kick%easy_axis(3)
      kick%trans_vec(3,1) = M_THREE*kick%easy_axis(1)

      dot = sum(kick%easy_axis(1:3)*kick%trans_vec(1:3,1))
      kick%trans_vec(1:3,1) = kick%trans_vec(1:3,1) - dot*kick%easy_axis(1:3)
      norm = sum(kick%trans_vec(1:3,1)**2)
      kick%trans_vec(1:3,1) = kick%trans_vec(1:3,1)/sqrt(norm)

      !To get a direct basis, the last vector is obtained by the cross product
      kick%trans_vec(1,2) = kick%easy_axis(2) * kick%trans_vec(3,1) - kick%easy_axis(3) * kick%trans_vec(2,1)
      kick%trans_vec(2,2) = kick%easy_axis(3) * kick%trans_vec(1,1) - kick%easy_axis(1) * kick%trans_vec(3,1)
      kick%trans_vec(3,2) = kick%easy_axis(1) * kick%trans_vec(2,1) - kick%easy_axis(2) * kick%trans_vec(1,1)

      !The perturbation direction is defined as
      !cos(q.r)*uvec + sin(q.r)*vvec


      if(parse_is_defined(namespace, 'TDMomentumTransfer') &
            .and. parse_is_defined(namespace, 'TDMultipleMomentumTransfer')) then
        message(1) = "TDMomentumTransfer and TDMultipleMomentumTransfer cannot be defined at the same time."
        call messages_fatal(1, namespace=namespace)
      end if

      if(parse_is_defined(namespace, 'TDMultipleMomentumTransfer')) then

        kick%qkick_mode = QKICKMODE_EXP

        !%Variable TDMultipleMomentumTransfer
        !%Type block
        !%Section Time-Dependent::Response
        !%Description
        !% For magnon kicks only.
        !% A simple way to specify momentum-transfer vectors for the calculation of
        !% the magnetization dynamics. This variable should be used for a supercell.
        !% For each reciprocal lattice vectors, the code will kick the original magnetization
        !% using all the multiples of it.
        !% The syntax reads:
        !%
        !% <tt>%TDMultipleMomentumTransfer
        !% <br>&nbsp;&nbsp;N_x | N_y | N_z
        !% <br>%</tt>
        !%
        !% and will include the (2N_x+1)*(2N_y+1)*(2N_z+1) multiples vectors of the reciprocal
        !% lattice vectors of the current cell.
        !%End
        if(parse_block(namespace, 'TDMultipleMomentumTransfer', blk) /= 0) then
          write(message(1),'(a)') 'Internal error while reading TDMultipleMomentumTransfer.'
          call messages_fatal(1, namespace=namespace)
        end if

        do idir = 1, 3
          call parse_block_integer(blk, 0, idir-1, kick%nqmult(idir))
        end do

        call parse_block_end(blk)


        kick%nqvec = (2*kick%nqmult(1)+1)*(2*kick%nqmult(2)+1)*(2*kick%nqmult(3)+1)
        !qvector has been allocated by default to a null vector before
        SAFE_DEALLOCATE_A(kick%qvector)
        SAFE_ALLOCATE(kick%qvector(1:MAX_DIM, 1:kick%nqvec))
        iq = 0
        do iqx = -kick%nqmult(1), kick%nqmult(1)
          do iqy = -kick%nqmult(2), kick%nqmult(2)
            do iqz = -kick%nqmult(3), kick%nqmult(3)
              iq = iq + 1
              qtemp(1:3) = (/iqx, iqy, iqz/)
              call kpoints_to_absolute(kpoints%latt, qtemp, kick%qvector(1:3, iq))

              !Checking symmetries for all G vectors
              if(kpoints%use_symmetries) then
                do iop = 1, symmetries_number(kpoints%symm)
                  if(iop == symmetries_identity_index(kpoints%symm)) cycle
                  if(.not. symm_op_invariant_cart(kpoints%symm%ops(iop), kick%qvector(:,iq), CNST(1e-5))) then
                    message(1) = "The TDMultipleMomentumTransfer breaks (at least) one " &
                                      // "of the symmetries used to reduce the k-points."
                    message(2) = "Set SymmetryBreakDir accordingly."
                    call messages_fatal(2, namespace=namespace)
                  end if
                end do
              end if
            end do
          end do
        end do

      end if

    else
      kick%easy_axis(1:MAX_DIM) = M_ZERO
    end if

    if(kick%delta_strength_mode == KICK_MAGNON_MODE .and. kick%qkick_mode /= QKICKMODE_EXP) then
      message(1) = "For magnons, the kick mode must be exponential."
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(kick_init)
  end subroutine kick_init


  ! ---------------------------------------------------------
  subroutine kick_copy(kick_out, kick_in)
    type(kick_t), intent(inout) :: kick_out
    type(kick_t), intent(in)    :: kick_in

    PUSH_SUB(kick_copy)

    kick_out%dim = kick_in%dim

    kick_out%time = kick_in%time

    kick_out%delta_strength_mode = kick_in%delta_strength_mode
    kick_out%delta_strength = kick_in%delta_strength

    !> In case we use a normal dipole kick:
    kick_out%pol(1:MAX_DIM, 1:MAX_DIM) = kick_in%pol(1:MAX_DIM, 1:MAX_DIM)
    kick_out%pol_dir = kick_in%pol_dir
    kick_out%pol_equiv_axes = kick_in%pol_equiv_axes
    kick_out%wprime(1:MAX_DIM) = kick_in%wprime(1:MAX_DIM)

    !> In case we have a general multipolar kick,
    kick_out%n_multipoles = kick_in%n_multipoles
    if (kick_out%n_multipoles > 0) then
      SAFE_ALLOCATE(kick_out%l(1:kick_out%n_multipoles))
      SAFE_ALLOCATE(kick_out%m(1:kick_out%n_multipoles))
      SAFE_ALLOCATE(kick_out%weight(1:kick_out%n_multipoles))
      kick_out%l = kick_in%l
      kick_out%m = kick_in%m
      kick_out%weight = kick_in%weight
    end if
    kick_out%nqvec = kick_in%nqvec
    SAFE_ALLOCATE(kick_out%qvector(1:MAX_DIM, 1:kick_in%nqvec))
    kick_out%qvector(1:MAX_DIM, 1:kick_in%nqvec) = kick_in%qvector(1:MAX_DIM, 1:kick_in%nqvec)
    kick_out%qlength = kick_in%qlength
    kick_out%qkick_mode = kick_in%qkick_mode
    kick_out%qbessel_l = kick_in%qbessel_l
    kick_out%qbessel_m = kick_in%qbessel_m

    !> In case we use a general function
    kick_out%function_mode = kick_in%function_mode
    kick_out%user_defined_function = kick_in%user_defined_function

    kick_out%easy_axis(1:MAX_DIM) = kick_in%easy_axis(1:MAX_DIM)

    POP_SUB(kick_copy)
  end subroutine kick_copy

  ! ---------------------------------------------------------
  subroutine kick_end(kick)
    type(kick_t), intent(inout) :: kick

    PUSH_SUB(kick_end)

    kick%delta_strength_mode = 0
    kick%dim = 0
    kick%pol_equiv_axes = 0
    kick%pol = M_ZERO
    kick%pol_dir = 0
    kick%wprime = M_ZERO
    if (kick%n_multipoles > 0) then
      SAFE_DEALLOCATE_A(kick%l)
      SAFE_DEALLOCATE_A(kick%m)
      SAFE_DEALLOCATE_A(kick%weight)
    end if
    kick%n_multipoles = 0
    kick%qkick_mode = QKICKMODE_NONE
    SAFE_DEALLOCATE_A(kick%qvector)
    kick%easy_axis(1:MAX_DIM) = M_ZERO

    POP_SUB(kick_end)
  end subroutine kick_end

  ! ---------------------------------------------------------
  subroutine kick_read(kick, iunit, namespace)
    type(kick_t),      intent(inout) :: kick
    integer,           intent(in)    :: iunit
    type(namespace_t), intent(in)    :: namespace

    integer :: idir, im, ierr
    character(len=100) :: line

    PUSH_SUB(kick_read)

    kick%function_mode = -1

    read(iunit, '(15x,i2)')     kick%delta_strength_mode
    read(iunit, '(15x,f18.12)') kick%delta_strength
    read(iunit, '(15x,i2)')     kick%dim
    read(iunit, '(a)') line
    if(index(line,'defined') /= 0) then
      kick%function_mode = KICK_FUNCTION_USER_DEFINED
      ! "# User defined: "
      read(line,'(16x,a)') kick%user_defined_function
    elseif(index(line,'multipole') /= 0) then
      kick%function_mode = KICK_FUNCTION_MULTIPOLE
      ! "# N multipoles "
      read(line, '(15x,i3)') kick%n_multipoles
      SAFE_ALLOCATE(     kick%l(1:kick%n_multipoles))
      SAFE_ALLOCATE(     kick%m(1:kick%n_multipoles))
      SAFE_ALLOCATE(kick%weight(1:kick%n_multipoles))
      do im = 1, kick%n_multipoles
        ! "# multipole    "
        read(iunit, '(15x,2i3,f18.12)') kick%l(im), kick%m(im), kick%weight(im)
      end do
    else
      kick%function_mode = KICK_FUNCTION_DIPOLE
      kick%n_multipoles = 0
      backspace(iunit)

      do idir = 1, kick%dim
        read(iunit, '(15x,99f18.12)') kick%pol(1:kick%dim, idir)
      end do
      read(iunit, '(15x,i2)')      kick%pol_dir
      read(iunit, '(15x,i2)')      kick%pol_equiv_axes
      read(iunit, '(15x,99f18.12)') kick%wprime(1:kick%dim)
    end if
    if(kick%delta_strength_mode == KICK_MAGNON_MODE) then
      read(iunit, '(15x,i3)') kick%nqvec
      SAFE_ALLOCATE(kick%qvector(1:MAX_DIM, 1:kick%nqvec))
      do im = 1, kick%nqvec
        read(iunit, '(15x,3f18.12)') kick%qvector(1:3, im)
      end do
      read(iunit, '(15x,3f18.12)')   kick%easy_axis(1:3)
      read(iunit, '(15x,3f18.12)')   kick%trans_vec(1:3,1)
      read(iunit, '(15x,3f18.12)')   kick%trans_vec(1:3,2)
    end if
    read(iunit, '(15x,f18.12)', iostat = ierr) kick%time

    if(ierr /= 0) then
      kick%time = M_ZERO
      backspace(iunit)
    end if

    if (kick%function_mode < 0) then
      message(1) = "No kick could be read from file."
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(kick_read)
  end subroutine kick_read


  ! ---------------------------------------------------------
  subroutine kick_write(kick, iunit, out)
    type(kick_t),          intent(in)    :: kick
    integer,    optional,  intent(in)    :: iunit
    type(c_ptr), optional, intent(inout) :: out

    integer :: idir, im
    character(len=120) :: aux

    PUSH_SUB(kick_write)

    if(present(iunit)) then
      write(iunit, '(a15,i1)')      '# kick mode    ', kick%delta_strength_mode
      write(iunit, '(a15,f18.12)')  '# kick strength', kick%delta_strength
      write(iunit, '(a15,i2)')      '# dim          ', kick%dim
      ! if this were to be read by humans, we would want units_from_atomic(units_out%length**(-1))
      if(kick%function_mode  ==  KICK_FUNCTION_USER_DEFINED) then
        write(iunit,'(a15,1x,a)')     '# User defined:', trim(kick%user_defined_function)
      elseif(kick%n_multipoles > 0) then
        write(iunit, '(a15,i3)')    '# N multipoles ', kick%n_multipoles
        do im = 1, kick%n_multipoles
          write(iunit, '(a15,2i3,f18.12)') '# multipole    ', kick%l(im), kick%m(im), kick%weight(im)
        end do
      else
        do idir = 1, kick%dim
          write(iunit, '(a6,i1,a8,99f18.12)') '# pol(', idir, ')       ', kick%pol(1:kick%dim, idir)
        end do
        write(iunit, '(a15,i1)')      '# direction    ', kick%pol_dir
        write(iunit, '(a15,i1)')      '# Equiv. axes  ', kick%pol_equiv_axes
        write(iunit, '(a15,99f18.12)') '# wprime       ', kick%wprime(1:kick%dim)
      end if
      write(iunit, '(a15,f18.12)') "# kick time    ", kick%time

    else if(present(out)) then
      write(aux, '(a15,i2)')      '# kick mode    ', kick%delta_strength_mode
      call write_iter_string(out, aux)
      call write_iter_nl(out)
      write(aux, '(a15,f18.12)')  '# kick strength', kick%delta_strength
      call write_iter_string(out, aux)
      call write_iter_nl(out)
      write(aux, '(a15,i2)')      '# dim          ', kick%dim
      call write_iter_string(out, aux)
      call write_iter_nl(out)
      if(kick%function_mode  ==  KICK_FUNCTION_USER_DEFINED) then
        write(aux,'(a15,1x,a)')     '# User defined:', trim(kick%user_defined_function)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
      elseif(kick%n_multipoles > 0) then
        write(aux, '(a15,i3)')      '# N multipoles ', kick%n_multipoles
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        do im = 1, kick%n_multipoles
          write(aux, '(a15,2i3,f18.12)') '# multipole    ', kick%l(im), kick%m(im), kick%weight(im)
          call write_iter_string(out, aux)
          call write_iter_nl(out)
        end do
      else
        do idir = 1, kick%dim
          write(aux, '(a6,i1,a8,99f18.12)') '# pol(', idir, ')       ', kick%pol(1:kick%dim, idir)
          call write_iter_string(out, aux)
          call write_iter_nl(out)
        end do
        write(aux, '(a15,i2)')      '# direction    ', kick%pol_dir
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,i2)')      '# Equiv. axes  ', kick%pol_equiv_axes
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,99f18.12)') '# wprime       ', kick%wprime(1:kick%dim)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
      end if
      if(present(out) .and. kick%delta_strength_mode == KICK_MAGNON_MODE) then
        write(aux, '(a15,i3)')      '# N q-vectors  ', kick%nqvec
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        do im = 1, kick%nqvec
          write(aux, '(a15,3f18.12)') '# q-vector     ', kick%qvector(1:3, im)
          call write_iter_string(out, aux)
          call write_iter_nl(out)
        end do
        write(aux, '(a15,3f18.12)')   '# Easy axis    ', kick%easy_axis(1:3)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,3f18.12)')   '# Trans. dir 1 ', kick%trans_vec(1:3,1)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,3f18.12)')   '# Trans. dir 2 ', kick%trans_vec(1:3,2)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
      end if
      write(aux, '(a15,f18.12)') "# kick time    ", kick%time
      call write_iter_string(out, aux)
      call write_iter_nl(out)

    end if

    POP_SUB(kick_write)
  end subroutine kick_write


  ! ---------------------------------------------------------
  !
  subroutine kick_function_get(space, mesh, kick, kick_function, iq, to_interpolate)
    type(space_t),        intent(in)    :: space
    type(mesh_t),         intent(in)    :: mesh
    type(kick_t),         intent(in)    :: kick
    CMPLX,                intent(out)   :: kick_function(:)
    integer,              intent(in)    :: iq
    logical, optional,    intent(in)    :: to_interpolate

    integer :: ip, im
    FLOAT   :: xx(space%dim)
    FLOAT   :: rkick, ikick, rr, ylm

    integer :: np

    PUSH_SUB(kick_function_get)

    np = mesh%np
    if(present(to_interpolate)) then
      if(to_interpolate) np = mesh%np_part
    end if

    if(abs(kick%qlength) > M_EPSILON .or. kick%delta_strength_mode == KICK_MAGNON_MODE) then ! q-vector is set
      ASSERT(space%dim == 3)

      select case (kick%qkick_mode)
        case (QKICKMODE_COS)
          write(message(1), '(a,3F9.5,a)') 'Info: Using cos(q.r) field with q = (', kick%qvector(1:3, iq), ')'
        case (QKICKMODE_SIN)
          write(message(1), '(a,3F9.5,a)') 'Info: Using sin(q.r) field with q = (', kick%qvector(1:3, iq), ')'
        case (QKICKMODE_SIN + QKICKMODE_COS)
          write(message(1), '(a,3F9.5,a)') 'Info: Using sin(q.r)+cos(q.r) field with q = (', kick%qvector(1:3, iq), ')'
        case (QKICKMODE_EXP)
          write(message(1), '(a,3F9.5,a)') 'Info: Using exp(iq.r) field with q = (', kick%qvector(1:3, iq), ')'
        case (QKICKMODE_BESSEL)
          write(message(1), '(a,I2,a,I2,a,F9.5)') 'Info: Using j_l(qr)*Y_lm(r) field with (l,m)= (', &
            kick%qbessel_l, ",", kick%qbessel_m,') and q = ', kick%qlength
        case default
           write(message(1), '(a,3F9.6,a)') 'Info: Unknown field type!'
      end select
      call messages_info(1)

      kick_function = M_z0
      do ip = 1, np
        xx = mesh%x(ip, :)
        select case (kick%qkick_mode)
          case (QKICKMODE_COS)
            kick_function(ip) = kick_function(ip) + cos(sum(kick%qvector(1:3, iq) * xx(1:3)))
          case (QKICKMODE_SIN)
            kick_function(ip) = kick_function(ip) + sin(sum(kick%qvector(1:3, iq) * xx(1:3)))
          case (QKICKMODE_SIN+QKICKMODE_COS)
            kick_function(ip) = kick_function(ip) + sin(sum(kick%qvector(1:3, iq) * xx(1:3)))
          case (QKICKMODE_EXP)
            kick_function(ip) = kick_function(ip) + exp(M_zI * sum(kick%qvector(1:3, iq) * xx(1:3)))
          case (QKICKMODE_BESSEL)
            call grylmr(mesh%x(ip, 1), mesh%x(ip, 2), mesh%x(ip, 3), kick%qbessel_l, kick%qbessel_m, ylm)
            kick_function(ip) = kick_function(ip) + loct_sph_bessel(kick%qbessel_l, kick%qlength*sqrt(sum(xx(:)**2)))*ylm
        end select
      end do

    else
      if(kick%function_mode  ==  KICK_FUNCTION_USER_DEFINED) then

        kick_function = M_z0
        do ip = 1, np
          call mesh_r(mesh, ip, rr, coords = xx)
            rkick = M_ZERO; ikick = M_ZERO
          call parse_expression(rkick, ikick, space%dim, xx, rr, M_ZERO, trim(kick%user_defined_function))
            kick_function(ip) = rkick
        end do

      elseif(kick%n_multipoles > 0) then

        kick_function = M_z0
        do im = 1, kick%n_multipoles
          do ip = 1, np
            call mesh_r(mesh, ip, rr, coords = xx)
            call loct_ylm(1, xx(1), xx(2), xx(3), kick%l(im), kick%m(im), ylm)
              kick_function(ip) = kick_function(ip) + kick%weight(im) * (rr**kick%l(im)) * ylm
          end do
        end do
      else
        do ip = 1, np
          kick_function(ip) = sum(mesh%x(ip, 1:space%dim) * &
            kick%pol(1:space%dim, kick%pol_dir))
        end do
      end if
    end if

    POP_SUB(kick_function_get)
  end subroutine kick_function_get


  ! ---------------------------------------------------------
  !
  subroutine kick_pcm_function_get(space, mesh, kick, psolver, pcm, kick_pcm_function)
    type(space_t),        intent(in)    :: space
    type(mesh_t),         intent(in)    :: mesh
    type(kick_t),         intent(in)    :: kick
    type(poisson_t),      intent(in)    :: psolver
    type(pcm_t),          intent(inout) :: pcm
    CMPLX,                intent(out)   :: kick_pcm_function(:)

    CMPLX, allocatable :: kick_function_interpolate(:)
    FLOAT, allocatable :: kick_function_real(:)

    PUSH_SUB(kick_pcm_function_get)

    kick_pcm_function = M_ZERO
    if ( pcm%localf ) then
      SAFE_ALLOCATE(kick_function_interpolate(1:mesh%np_part))
      kick_function_interpolate = M_ZERO
      call kick_function_get(space, mesh, kick, kick_function_interpolate, 1, to_interpolate = .true.)
      SAFE_ALLOCATE(kick_function_real(1:mesh%np_part))
      kick_function_real = DREAL(kick_function_interpolate)
      if ( pcm%kick_like ) then
        ! computing kick-like polarization due to kick
        call pcm_calc_pot_rs(pcm, mesh, psolver, kick = kick%delta_strength * kick_function_real, kick_time = .true.)
      else if ( .not.pcm%kick_like .and. pcm%which_eps == PCM_DEBYE_MODEL) then
        ! computing the kick-like part of polarization due to kick for Debye dielectric model
        pcm%kick_like = .true.
        call pcm_calc_pot_rs(pcm, mesh, psolver, kick = kick%delta_strength * kick_function_real, kick_time = .true.)
        pcm%kick_like = .false.
      else if ( .not.pcm%kick_like .and. pcm%which_eps == PCM_DRUDE_MODEL ) then
        POP_SUB(kick_pcm_function_get)
        return
      end if
      kick_pcm_function = pcm%v_kick_rs / kick%delta_strength
    end if

    POP_SUB(kick_pcm_function_get)
  end subroutine kick_pcm_function_get


  ! ---------------------------------------------------------
  !> Applies the delta-function electric field \f$ E(t) = E_0 \Delta(t) \f$
  !! where \f$ E_0 = \frac{- k \hbar}{e} \f$ k = kick\%delta_strength.
  subroutine kick_apply(space, mesh, st, ions_dyn, ions, kick, psolver, kpoints, pcm)
    type(space_t),         intent(in)    :: space
    type(mesh_t),          intent(in)    :: mesh
    type(states_elec_t),   intent(inout) :: st
    type(ion_dynamics_t),  intent(in)    :: ions_dyn
    type(ions_t),          intent(inout) :: ions
    type(kick_t),          intent(in)    :: kick
    type(poisson_t),       intent(in)    :: psolver
    type(kpoints_t),       intent(in)    :: kpoints
    type(pcm_t), optional, intent(inout) :: pcm

    integer :: iqn, ist, idim, ip, ispin, iatom
    CMPLX   :: cc(2), kick_value
    CMPLX, allocatable :: kick_function(:), psi(:, :)

    CMPLX, allocatable :: kick_pcm_function(:)
    integer :: ns, iq
    FLOAT :: uvec(MAX_DIM), vvec(MAX_DIM), Gvec(MAX_DIM,MAX_DIM)
    FLOAT :: xx(MAX_DIM), rr

    PUSH_SUB(kick_apply)

    ! The wavefunctions at time delta t read
    ! psi(delta t) = psi(t) exp(i k x)
    delta_strength: if(kick%delta_strength /= M_ZERO) then

      SAFE_ALLOCATE(kick_function(1:mesh%np))
      if(kick%delta_strength_mode /= KICK_MAGNON_MODE .or. kick%nqvec == 1) then
        call kick_function_get(space, mesh, kick, kick_function, 1)
      end if

      ! PCM - computing polarization due to kick
      if( present(pcm) ) then
        SAFE_ALLOCATE(kick_pcm_function(1:mesh%np))
        call kick_pcm_function_get(space, mesh, kick, psolver, pcm, kick_pcm_function)
        kick_function = kick_function + kick_pcm_function
      end if

      write(message(1),'(a,f11.6)') 'Info: Applying delta kick: k = ', kick%delta_strength
      select case (kick%function_mode)
      case (KICK_FUNCTION_DIPOLE)
        message(2) = "Info: kick function: dipole."
      case (KICK_FUNCTION_MULTIPOLE)
        message(2) = "Info: kick function: multipoles."
      case (KICK_FUNCTION_USER_DEFINED)
        message(2) = "Info: kick function: user defined function."
      end select
      select case (kick%delta_strength_mode)
      case (KICK_DENSITY_MODE)
        message(3) = "Info: Delta kick mode: Density mode"
      case (KICK_SPIN_MODE)
        message(3) = "Info: Delta kick mode: Spin mode"
      case (KICK_SPIN_DENSITY_MODE)
        message(3) = "Info: Delta kick mode: Density + Spin modes"
      end select
      call messages_info(3)

      ns = 1
      if(st%d%nspin == 2) ns = 2

      SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))

      if(kick%delta_strength_mode /= KICK_MAGNON_MODE) then

        do iqn = st%d%kpt%start, st%d%kpt%end
          do ist = st%st_start, st%st_end
            call states_elec_get_state(st, mesh, ist, iqn, psi)

            select case (kick%delta_strength_mode)
            case (KICK_DENSITY_MODE)
              do idim = 1, st%d%dim
                do ip = 1, mesh%np
                  psi(ip, idim) = exp(M_zI*kick%delta_strength*kick_function(ip))*psi(ip, idim)
                end do
              end do

            case (KICK_SPIN_MODE)
              ispin = st%d%get_spin_index(iqn)
              do ip = 1, mesh%np
                kick_value = M_zI*kick%delta_strength*kick_function(ip)

                cc(1) = exp(kick_value)
                cc(2) = exp(-kick_value)

                select case (st%d%ispin)
                case (SPIN_POLARIZED)
                  psi(ip, 1) = cc(ispin)*psi(ip, 1)
                case (SPINORS)
                  psi(ip, 1) = cc(1)*psi(ip, 1)
                  psi(ip, 2) = cc(2)*psi(ip, 2)
                end select
              end do

            case (KICK_SPIN_DENSITY_MODE)
              do ip = 1, mesh%np
                kick_value = M_zI*kick%delta_strength*kick_function(ip)
                cc(1) = exp(M_TWO*kick_value)

                select case (st%d%ispin)
                case (SPIN_POLARIZED)
                  if(is_spin_up(iqn)) then
                    psi(ip, 1) = cc(1)*psi(ip, 1)
                  end if
                case (SPINORS)
                  psi(ip, 1) = cc(1)*psi(ip, 1)
                end select
              end do
            end select

            call states_elec_set_state(st, mesh, ist, iqn, psi)

          end do
        end do

      else
        ASSERT(st%d%ispin==SPINORS)

        if(kick%nqvec == 1) then
          !The perturbation direction is defined as
          !cos(q.r)*uvec + sin(q.r)*vvec
          uvec(1:3) = kick%trans_vec(1:3,1)
          vvec(1:3) = kick%trans_vec(1:3,2)

          do iqn = st%d%kpt%start, st%d%kpt%end, ns
            do ist = st%st_start, st%st_end

              call states_elec_get_state(st, mesh, ist, iqn, psi)

              do ip = 1, mesh%np

                cc(1) = psi(ip, 1)
                cc(2) = psi(ip, 2)

                !First part: 1I*cos(\lambda)
                psi(ip, 1) = cos(kick%delta_strength)* cc(1)
                psi(ip, 2) = cos(kick%delta_strength)* cc(2)

                !We now add -i sin(\lambda) u.\sigma
                !           (u_z      u_x-i*u_y)            (v_z         v_x-i*v_y)
                ! =cos(q.r) (                  )  + sin(q.r)(                     )
                !           (u_x+i*u_y  -u_z   )            (v_x+i*v_y   -v_z     )
                psi(ip, 1) = psi(ip, 1) -M_zI*sin(kick%delta_strength)*( TOFLOAT(kick_function(ip)) &
                                  * (uvec(3)*cc(1) + (uvec(1)-M_zI*uvec(2))*cc(2)) &
                       + aimag(kick_function(ip)) * (vvec(3)*cc(1) + (vvec(1)-M_zI*vvec(2))*cc(2)))
                psi(ip, 2) = psi(ip, 2) -M_zI*sin(kick%delta_strength)*( TOFLOAT(kick_function(ip)) &
                                  * (-uvec(3)*cc(2) + (uvec(1)+M_zI*uvec(2))*cc(1)) &
                       + aimag(kick_function(ip)) * (-vvec(3)*cc(2) + (vvec(1)+M_zI*vvec(2))*cc(1)))

              end do

              call states_elec_set_state(st, mesh, ist, iqn, psi)

            end do
          end do

        else ! Multi-q kick

           call kpoints_to_absolute(kpoints%latt, (/M_ONE,M_ZERO,M_ZERO/), Gvec(1:3, 1))
           call kpoints_to_absolute(kpoints%latt, (/M_ZERO,M_ONE,M_ZERO/), Gvec(1:3, 2))
           call kpoints_to_absolute(kpoints%latt, (/M_ZERO,M_ZERO,M_ONE/), Gvec(1:3, 3))

           kick_function = M_ONE
           do ip = 1, mesh%np
             call mesh_r(mesh, ip, rr, coords = xx)
             do iq = 1, 3
               if(kick%nqmult(iq) == 0) cycle
               if(abs(sin(M_HALF*sum(Gvec(1:3, iq) * xx(1:3)))) <= M_EPSILON) cycle

               kick_function(ip) = kick_function(ip)*sin(M_HALF*(2*kick%nqmult(iq)+1) &
                     *sum(Gvec(1:3, iq) * xx(1:3)))/sin(M_HALF*sum(Gvec(1:3, iq) * xx(1:3)))
             end do
             kick_function(ip) = kick_function(ip)*kick%delta_strength
           end do

           do iqn = st%d%kpt%start, st%d%kpt%end, ns
            do ist = st%st_start, st%st_end

              call states_elec_get_state(st, mesh, ist, iqn, psi)

              do ip = 1, mesh%np

                cc(1) = psi(ip, 1)
                cc(2) = psi(ip, 2)

                !   (cos(F) + in_x sin(F)                   sin(F)(u_y (u_x-iu_y)/(1+u_z) - iu_z))
                ! = (                                                                            )
                !   (-sin(F)(u_y (u_x+iu_y)/(1+u_z)+iu_z)   cos(F) - in_x sin(F)                 )

                psi(ip, 1) = (cos(kick_function(ip))+M_zI*kick%easy_axis(1)*sin(kick_function(ip)))*cc(1) &
                        +sin(kick_function(ip))*(kick%easy_axis(2)*(kick%easy_axis(1) &
                        -M_zI*kick%easy_axis(2))/(1+kick%easy_axis(3))-M_zI*kick%easy_axis(3))*cc(2)
                psi(ip, 2) =-sin(kick_function(ip))*(kick%easy_axis(2)*(kick%easy_axis(1) &
                        +M_zI*kick%easy_axis(2))/(1+kick%easy_axis(3))+M_zI*kick%easy_axis(3))*cc(1) &
                        + (cos(kick_function(ip))-m_zI*kick%easy_axis(1)*sin(kick_function(ip)))*cc(2)

              end do

              call states_elec_set_state(st, mesh, ist, iqn, psi)

            end do
          end do

        end if
      end if

      SAFE_DEALLOCATE_A(psi)

      ! The nuclear velocities will be changed by
      ! Delta v_z = ( Z*e*E_0 / M) = - ( Z*k*\hbar / M)
      ! where M and Z are the ionic mass and charge, respectively.
      if(ion_dynamics_ions_move(ions_dyn)  .and. kick%delta_strength /= M_ZERO) then
        if(kick%delta_strength_mode /= KICK_MAGNON_MODE) then
          do iatom = 1, ions%natoms
            ions%vel(:, iatom) = ions%vel(:, iatom) + &
              kick%delta_strength * kick%pol(1:ions%space%dim, kick%pol_dir) * &
              P_PROTON_CHARGE * species_zval(ions%atom(iatom)%species) / ions%mass(iatom)
          end do
        end if
      end if

      SAFE_DEALLOCATE_A(kick_function)
    end if delta_strength

    POP_SUB(kick_apply)
  end subroutine kick_apply

  pure integer function kick_get_type(kick) result(kick_type)
    type(kick_t),    intent(in) :: kick

    kick_type = kick%delta_strength_mode
 
  end function kick_get_type

end module kick_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
