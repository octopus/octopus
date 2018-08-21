!! Copyright (C) 2013 U. De Giovannini
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

module cmplxscl_oct_m
  use accel_oct_m
  use global_oct_m
  use parser_oct_m
  use messages_oct_m
  use varinfo_oct_m

  implicit none
  
  private

  public ::                        &
    cmplxscl_t,                    &
    cmplxscl_init,                 & 
    cmplxscl_copy,                 & 
    cmplxscl_end ,                 &
    cmplxscl_energy_ordering_score
       
  !> Complex scaling module    
  type cmplxscl_t
    logical            :: space       !< scale spatial coordinates
    logical            :: time        !< scale time coordinate
    FLOAT              :: theta       !< spatial coordinates scaling angle
    FLOAT              :: alphaR      !< time coordinate scaling angle for right states
    FLOAT              :: alphaL      !< time coordinate scaling angle for left states
    FLOAT              :: rotatespectrumangle !< angle with which to rotate eigenvalues to obtain desired order from ARPACK
    FLOAT              :: penalizationfactor  !< factor which penalizes imaginary parts of eigenvalues when ordering states
    FLOAT              :: localizationradius  !< integration radius used for distinguishing localized vs continuum states
    FLOAT              :: localizationthreshold !< states with delocalization larger than this value will be considered continuum
    integer            :: nlocalizedstates    !< number of states to be considered localized and thus preferentially occupied
  end type cmplxscl_t

  integer, parameter ::     &
    CMPLXSCL_NONE    = 0,   &
    CMPLXSCL_SPACE   = 2,   &
    CMPLXSCL_TIME    = 4

contains

  !-----------------------------------------------------------------
  subroutine cmplxscl_init(this)
    type(cmplxscl_t),    intent(out) :: this

    integer :: cmplxscl_flags
    
    PUSH_SUB(cmplxscl_init)

    !%Variable ComplexScaling
    !%Type flag
    !%Default none
    !%Section Hamiltonian::ComplexScaling
    !%Description
    !% (Experimental) Global complex scaling. The options 
    !% allow to scale space and time coordinates in the Hamiltonian.
    !% Cannot be used with OpenCL currently.
    !% You can specify more than one value by giving them as a sum, 
    !% for example:
    !% <tt>ComplexScaling = space + time</tt>
    !%Option none 0
    !% No scaling is applied. This is the default.
    !%Option space 2
    !% This is implements the global coordinate transformation <math>r \rightarrow r e^{i \theta}</math>.
    !% When <tt>TheoryLevel=DFT</tt>, density-functional resonance theory (DFRT) is employed.  
    !% In order to reveal resonances, <tt>ComplexScalingTheta</tt> bigger than zero should be set.
    !% D. L. Whitenack and A. Wasserman, <i>Phys. Rev. Lett.</i> <b>107</b>, 163002 (2011).
    !%Option time 4
    !% This is implements the coordinate transformation <math>t \rightarrow t e^{i \alpha_r}</math> for 
    !% right states and <math>t \rightarrow t e^{i \alpha_l}</math> for left states.
    !% J. Bengtsson, E. Lindroth, and S. Selst&oslash;, <i>Phys. Rev. A</i> <b>85</b>, 013419 (2012).
    !%End
    call parse_variable('ComplexScaling', CMPLXSCL_NONE, cmplxscl_flags)
    if(.not.varinfo_valid_option('ComplexScaling', cmplxscl_flags, is_flag = .true.)) then
      call messages_input_error('ComplexScaling')
    end if

    this%space = bitand(cmplxscl_flags, CMPLXSCL_SPACE) /= 0
    this%time  = bitand(cmplxscl_flags, CMPLXSCL_TIME)  /= 0

    if(this%space) then
      !%Variable ComplexScalingTheta
      !%Type float 
      !%Default 0.3
      !%Section Hamiltonian::ComplexScaling
      !%Description
      !% The spatial coordinate complex scaling angle <math>\theta</math>, in radians.
      !% Allowed values must be in the range <math>0 <= \theta < \pi/4</math>.
      !%End
      call parse_variable('ComplexScalingTheta', CNST(0.3), this%theta)
      if(this%theta < M_ZERO .or. this%theta > M_PI/CNST(4.0)) call messages_input_error('ComplexScalingTheta')
    else
      this%theta = M_ZERO
    end if

    !%Variable ComplexScalingRotateSpectrum
    !%Type float
    !%Default 0.0
    !%Section Hamiltonian::ComplexScaling
    !%Description
    !% The order of occupation of eigenvalues is by
    !% real part unless otherwise specified.  This parameter 
    !% is used to rotate the whole spectrum before assigning
    !% occupations, thus customizing the sorting scheme.
    !% The spectrum is rotated back afterwards.
    !%End
    call parse_variable('ComplexScalingRotateSpectrum', M_ZERO, this%rotatespectrumangle)

    !%Variable ComplexScalingPenalizationFactor
    !%Type float
    !%Default 2
    !%Section Hamiltonian::ComplexScaling
    !%Description
    !% Eigenvalues <math>\varepsilon</math> will be ordered by
    !% <math>\Re(\varepsilon) + {\rm penalizationfactor} (\Im(\varepsilon))^2</math>.
    !%End
    call parse_variable('ComplexScalingPenalizationFactor', M_TWO, this%penalizationfactor)



    !%Variable ComplexScalingLocalizedStates
    !%Type integer
    !%Default 0
    !%Section Hamiltonian::ComplexScaling
    !%Description
    !% It will be assumed that this many localized states are present in the calculation.  Localized
    !% states will be occupied as ordered by real part of the energy.  Any remaining states will then
    !% be ordered by localization as measured.  Localization is measured by integrating the square
    !% of each wavefunction within <tt>ComplexScalingLocalizationRadius</tt>.
    !%End    
    call parse_variable('ComplexScalingLocalizedStates', 0, this%nlocalizedstates)

    !%Variable ComplexScalingLocalizationThreshold
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian::ComplexScaling
    !%Description
    !% If the part of a state beyond <tt>ComplexScalingLocalizationRadius</tt> has a norm greater than this value,
    !% the state will be considered a continuum state.
    !%End
    call parse_variable('ComplexScalingLocalizationThreshold', M_ONE, this%localizationthreshold)

    !%Variable ComplexScalingAlpha
    !%Type float 
    !%Default 2*<tt>ComplexScalingTheta</tt>
    !%Section Hamiltonian::ComplexScaling
    !%Description
    !% The time coordinate complex scaling angle <math>\alpha_r</math> used to evolve 
    !% right states.  
    !%End
    if(this%time .and. this%space) then
      call parse_variable('ComplexScalingAlpha', M_TWO*this%theta, this%alphaR)
    else
      call parse_variable('ComplexScalingAlpha', M_ZERO, this%alphaR)
    end if

    !%Variable ComplexScalingAlphaLeft
    !%Type float 
    !%Default <tt>ComplexScalingAlpha</tt>
    !%Section Hamiltonian::ComplexScaling
    !%Description
    !% The time coordinate complex scaling angle <math>\alpha_l</math> used to evolve 
    !% left states.  
    !%End
    call parse_variable('ComplexScalingAlphaLeft', this%alphaR, this%alphaL)

    !%Variable ComplexScalingLocalizationRadius
    !%Type float
    !%Default 0.0
    !%Section Hamiltonian::ComplexScaling
    !%Description
    !% Radius used to determine whether a state is a continuum state (in bohr).
    !%End
    call parse_variable('ComplexScalingLocalizationRadius', M_ZERO, this%localizationradius)

    if (this%space .or. this%time) then
      call messages_print_stress(stdout, "Complex Scaling")
      call messages_experimental('Complex Scaling')
    end if
    
    if (this%space) then
      write(message(1), '(a)') 'Space complex scaling transformation: r -> r * e^(i * theta) '
      write(message(2), '(a,f7.4)') 'Complex scaling angle theta =  ', this%theta
      call messages_info(2)
    end if
      
    if (this%time) then
      write(message(1), '(a)') 'Time complex scaling transformation:  t -> t * e^(i * alpha[R,L]) '
      write(message(2), '(a,f7.4)') 'Complex scaling angle alphaR = ', this%alphaR
      write(message(3), '(a,f7.4)') 'Complex scaling angle alphaL = ', this%alphaL
      call messages_info(3)
    end if

    if (this%space .or. this%time) then
      call messages_print_stress(stdout)
    end if


    POP_SUB(cmplxscl_init)

  end subroutine cmplxscl_init


  !-----------------------------------------------------------------
  subroutine cmplxscl_end(this)
    type(cmplxscl_t),     intent(inout) :: this

    PUSH_SUB(cmplxscl_end)


    POP_SUB(cmplxscl_end)
  end subroutine cmplxscl_end


  !-----------------------------------------------------------------
  subroutine cmplxscl_copy(in, out)
    type(cmplxscl_t),     intent(in) :: in
    type(cmplxscl_t),     intent(out) :: out

    PUSH_SUB(cmplxscl_copy)

    out%space  = in%space
    out%time   = in%time
    out%theta  = in%theta
    out%alphaL = in%alphaL
    out%alphaR = in%alphaR


    POP_SUB(cmplxscl_copy)
  end subroutine cmplxscl_copy

  !-----------------------------------------------------------------
  function cmplxscl_energy_ordering_score(energy, penalizationfactor) result(score)
    CMPLX, intent(in) :: energy
    FLOAT, intent(in) :: penalizationfactor
    CMPLX             :: score
    
    PUSH_SUB(complex_energy_ordering_score)

    score = real(energy) + penalizationfactor * aimag(energy)**2
    if (aimag(energy).gt.M_ZERO) then
      ! Pretty arbitrary.  Must work on something more systematic
      score = score + CNST(8.0) * aimag(energy)
    end if
    
    POP_SUB(complex_energy_ordering_score)
    
  end function cmplxscl_energy_ordering_score

end module cmplxscl_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
