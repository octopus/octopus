---
Title: "Maxwell input file"
weight: 4
---

# Maxwell Input File

## Input file variable description

### Calculation mode and parallelization strategy  

At the beginning of the input file, the basic option variable "CalculationMode" selects the run mode of Octopus and has always to be set.  
In case of a parallel run, there are some variables to set the proper parallelization options.  
Additionally, some relevant physical constants or parameters can be defined.

{{< code-block >}}
#include_input doc/html/includes/01_input_file_calculation_mode.inp
{{< /code-block >}}


### Multisystem setup

The Maxwell system is now implemented in Octopus' multisystem framework, which allows to calculate several different systems, which can interact with each other.

Currently implemented system types are:

*  electronic: An electronic system. (NOT IMPLEMENTED)
*  maxwell: A maxwell system.
*  classical_particle: A classical particle. Used for testing purposes only.
*  charged_particle: A charged classical particle.
*  multisystem: A system containing other systems.


### Definition of the systems

In this tutorial, we will use a pure Maxwell system:

```sh
%Systems
  'Maxwell' | maxwell
%
```

Subsequently, variables relating to a named system are prefixed by the system name.

### Maxwell box variables and parameters

Similar to the usual matter box, the Maxwell box can be determined equivalently only using the prefix Maxwell. We note, that both system boxes can be selected completely independently, but they have the same center. The Maxwell box is only defined as a parallelepiped box.

!!! Example
    ```
    # ----- Maxwell box variables ----------------------------------------------
    
    # Variable of the Maxwell grid size
    lsize_mx = 12.0
    
    # Variable of the Maxwell grid spacing
    dx_mx = 9.5
    
    # Number of dimensions of the Maxwell grid
    Maxwell.Dimensions = 3
    
    # Box shape of the Maxwell simulation box
    Maxwell.BoxShape = parallelepiped
    
    # Box size of the Maxwell simulation box in each dimension
    
    %Maxwell.[Lsize]
    #
    
    
    
    lsize_mx | lsize_mx | lsize_mx
    %
    
    # Grid spacing of the Maxwell simulation box in each dimension
    
    %Maxwell.Spacing
    dx_mx | dx_mx | dx_mx
    %
    ```


### Maxwell options

The Maxwell options determines the Maxwell propagation scheme. First, the Maxwell propagation operator options consist of the type of "MaxwellHamiltonianOperator", which can be set as a vacuum Maxwell-Hamiltonian (faraday_ampere) or a Maxwell-Hamiltonian in linear medium (faraday_ampere_medium). The Maxwell-Hamiltonian Operator himself can be applied as a finite difference operator (fd) or as a phase factor in momentum space with a fast Fourier transform (fft). 

Currently, the Maxwell system does not define its own propagator, but uses the system defauls propagator, defined in TDSystemPropagator. So far, only the exponential midpoint is implemented.

The Maxwell Boundary conditions can be chosen for each direction differently with "MaxwellBoundaryConditions". It can be simulated zero and absorbing boundary conditions, as well as plane waves and constant fields boundary conditions. Additionally to the boundary condition, the code can include absorbing boundaries. This means that all outgoing waves are absorbed while all incoming signals still arises at the boundaries and propagate into the simulation box. The absorbing boundaries can be achieved by a mask function of by a perfectly matchen layer calculation with additional parameters.

	# ----- Maxwell calculation variables -----------------------

	# Type of Maxwell Hamiltonian operator

	# - faraday_ampere (3x3 for vacuum calculation)
	# - faraday_ampere_medium (6x6 for medium calculation)

	MaxwellHamiltonianOperator = faraday_ampere

	# Maxwell operator method
	# - MaxwellTDOperatorMethod op_fd (finite difference operation)
	# - MaxwellTDOperatorMethod op_fft (fft into the fourier space, phase multiplication, fft back into real space)

	MaxwellTDOperatorMethod = op_fd


	# Options for the Maxwell derivative stencil and the exponential expanison of the time propagator

	# - MaxwellDerivativesStencil = stencil_star
	# - MaxwellDerivativesStencil = stencil_starplus
	# - MaxwellDerivativesStencil = stencil_cube

	MaxwellDerivativesStencil = stencil_starplus
	MaxwellDerivativesOrder = 4
	MaxwellTDExpOrder = 4


	# Maxwell boundary conditions for each dimension
	# - MaxwellBoundaryConditions = zero        (Boundaries are set to zero.)
	# - MaxwellBoundaryConditions = constant    (Boundaries are set to a constant.)
	# - MaxwellBoundaryConditions = mirror_pec  (Perfect electric conductor.)
	# - MaxwellBoundaryConditions = mirror_pmc  (Perfect magnetic conductor.)
	# - MaxwellBoundaryConditions = plane_waves (Boundaries feed in plane waves.)

	%MaxwellBoundaryConditions
	plane_waves | plane_waves | plane_waves
	%


	# Maxwell absorbing boundaries options
	# - MaxwellAbsorbingBoundaries = not_absorbing (zero boundary condition)
	# - MaxwellAbsorbingBoundaries = mask (Maxwell field is muliplied by a mask function)
	# - MaxwellAbsorbingBoundaries = cpml (Maxwell PEML with convolution method)

	%MaxwellAbsorbingBoundaries
	not_absorbing | not_absorbing | not_absorbing
	%

	# Absorbing boundary width for the Maxwell mask function (if using masked boundaries)
	MaxwellABWidth = 5.0


	# Parameters to tune the Maxwell PML
	# - numerical tests show best performance for MaxwellABPMLPower between 2 and 3
	# - MaxwellABPMLReflectionError should be rather "small"

	MaxwellABPMLPower = 2.0
	MaxwellABPMLReflectionError = 1e-16


*[PML]: Perfectly Matched Layer
*[fft]: Fast Fourtier Transform


### Output options

The output option "OutputFormat" is valid for both systems and choses the type of output format. It is possible to chose multiple formats.

	# ----- General output variables ----------------------------------------------------------------

	# Output format
	# - Output format = plane_{axis} (with {axis}=x,y,z: plotting plane with {axis}=O)
	# - Output format = vtk (full 3D box values in vtk format)
	# - Output format = xyz (output format for coordinates of ions)
	# - Output format = axis_x (1D plot along x axix with y=O and z=O)
	OutputFormat= plane_ x + plane_ y + plane_ 2 + vtk + xyz + axis_x




### Maxwell output options

The Maxwell output options for the Maxwell field variables are equivalent to those of the matter output only with the prefix Maxell. The interval option is independent on the matter one.

	# ----- Maxwell output variables ----------------------------------------------------------------

	# Maxwell output variables are written every MaxwellOutputInterval time steps into the output_ iter folder
	# - MaxwellOutput = maxwell_electric_field (electric field output)
	# - MaxwellOutput = maxwell_magnetic_field (magnetic field output)
	# - MaxwellOutput = maxwell_trans_electric_field (transverse electric field output)
	# - MaxwellOutput = maxwell_trans_magnetic_field (transverse magnetic field output)
	# - MaxwellOutput = maxwell_long_electric_field (longitudinal electric field output)
	# - MaxwellOutput = maxwell_long_magnetic_field (longitudinal magnetic field output)
	# - MaxwellOutput = maxwell_div_electric_field (divergence of the electric field output)
	# - MaxwellOutput = maxwell_div_magnetic_field (divergence of the magnetic field output)
	# - MaxwellOutput = maxwell_vector_potential (vector potential output)
	# - MaxwellOutput = maxwell_poynting_vector (poynting vector output)
	# - MaxwellOutput = maxwell_energy_density (electromagnetic energy density output)
	# - MaxwellOutput = maxwell_current (electromagnetic current density on the Maxwell grid)
	# - MaxwellOutput = maxwell_external_current (external current density on the Maxwell grid)
	# - MaxwellOutput = maxwell_electric_dipole_potential (electric dipole potential output)
	# - MaxwellOutput = maxwell_electric_quadrupole_potential (electric quadrupole potential output)
	# - MaxwellOutput = maxwell_charge_density (electromagnetic charge density via div(E))
	# - MaxwellOutput = maxwell_charge_density_diff (difference between div(E) and KS charge density)

	MaxwellOutput = maxwell_electric_field + maxwell_magnetic_field

	# Output interval steps for the Maxwell variables. After MaxwellOutputInterval steps,
	# the Maxwell variables on the grid are written into the output_iter folder
	MaxwellOutputInterval = 1

	# Output of the scalar Maxwell variables for each time step, written into the td.general folder

	# - MaxwellTDOutput = maxwell_energy (electromagnetic and Maxwell-matter energies)
	# - MaxwellTDOutput = maxwell_spectrum (variables to calculate the spectrum via the Maxwell system)
	# - MaxwellTDOutput = maxwell_fields (Mx fields at coordinates given by MaxwellFieldsCoordinate)
	# - MaxwellTDOutput = maxwell_mean_poynting (mean poynting vector in each direction)
	# - MaxwellTDOutput = maxwell_poynting_surface (mean poynting vector of boundary surfaces)
	# - MaxwellTDOutput = maxwell_e_field_surface (mean electric field vector of boundary surfaces)
	# - MaxwellTDOutput = maxwell_b_field_surface (mean magnetic field vector of boundary surfaces)

	MaxwellTDOutput = maxwell_energy + maxwell_fields

	# Coordinates of the grid points, which corresponding electromagnetic field values are
	# written into td.general/maxwell_fields
	%MaxwellFieldsCoordinate
	0.00 | 0.00 | 0.00
	%


#### Time step variables  

The "TDTimeStep" option defines the timestep for one Maxwell propagation step in case of only a free Maxwell propagation with one Maxwell system or it describes the time step for one propagation step of the Matter system.  
In general, the matter system is more stable for larger time steps, whereas the Maxwell system has always to fulfill the Courant condition. To handle both in a feasible way, the Maxwell system can apply several time steps between one matter time step. The "MaxwellTDIntervalSteps" gives the number of Maxwell intermediate steps between two matter steps.  
To speed up the code it is useful to add another approximation that the current density assumes to be constant as the arithmetic mean between two matter time steps. This approximation can be applied with the "MaxwellTDETRSApprox" option.  
The two remaining options "TDMaxwellTDRelaxationSteps" and "TDMaxwellKSRelaxationSteps" can be used if the starting groundstate is not stable enough, e.g. in case of absorbing boundaries for the matter wavefunction. Since the groundstate cacluations are always done with zero boundary condition, switching to absorbing boundaries in a time-depending mode leads to small dynamics until the system relaxes to its groundstate with absorbing boundaries. Sometimes this relaxation causes strong currents that lead to instable runs. Using "TDMaxwellTDRelaxationSteps" the simulation avoids that the relaxation dynamics effects the Maxwell fields. After that number of time steps, the matter to Maxwell coupling is switched on. The second number "TDMaxwellKSRelaxationSteps" describes the number when the external fields like laser pulse or external current are switched on. In other words the full requested simulation starts after TDMaxwellKSRelaxationSteps effectrive zero time step.

	# ----- Time step variables ---------------------------------------------------------------------

	# Time step of the propagation
	# In case of coupled Maxwell-Kohn Sham propagation, this variable is the matter timestep.
	# The Maxwell system time step is given by TDTimeStep/MaxwellTDIntervalSteps
	# so that the Maxwell system propagates several intermediate steps between consecutive matter steps.
	#
	# TDTimeStep should be equal or smaller than the Courant criterion, which is here
	# S_Courant = 1 / (sqrt(c^2/dx_mx^2 + c^2/dx_mx^2 + c^2/dx_mx^2))

	TDTimeStep = 0.002

	# Total number of Matter steps to run
	TDMaxSteps = 200

	# Maxwell intermediate time steps only plays a role in the maxwell_ks mode. It
	# indicates the number of Maxwell time propagations between consecutive matter steps
	MaxwellTDIntervalSteps = 1

	# The Maxwell TD relaxation steps are the number of steps in the maxwell_ks mode,
	# where the run propagates only in a standard TD mode to relax the starting
	# groundstate if necessary.

	TDMaxwellTDRelaxationSteps = 0

	# The Maxwell ks relaxation steps denote the number of steps in the maxwell_ks mode,
	# where the run propagates in a maxwell_ks mode without any external field to relax
	# the starting groundstate if necessary.

	TDMaxwellKSRelaxationSteps = 0

	# Approximation method for the Maxwell ETRS propagator
	# MaxwellTDETRSApprox = no (no approximation)
	# MaxwellTDETRSApprox = const_steps (assumption that the current density is approximately constant
	# between two matter time steps and using the arithmetic mean of the current density at the latest # matter time step and the next matter time step)
	MaxwellTDETRSApprox = no

	# Internal current test run, the test function is hard coded in current.F90
	CurrentPropagationTest = no


#### Maxwell field variables  

The initial Maxwell field can be manually defined using the "UserDefinedInitialMaxwellStates" block.
The first column represents the input type, here a predefined function instead of parsing a function. The second to the fourth column define the maximum electric field amplitude vector of the plane wave. The fifth column names the following defined Maxwell envelope function, and the last code denotes that the incident wave is a plane wave. Multiple different plane waves can be defined, and their superposition gives the final incident wave result at the boundary.  
The MaxwellFunctions block reads the additional required informations for the plane wave pulse. The first column refers to the envelope function name in the "UserDefinedInitialMaxwellStates" block. The second column gives the type of spatial Maxwell function. The next three columns represent the wavevector, followed by the three columns for the shift vector in space to move the pulse in space. The last column defines the width of the corresponding envelope function. In case of plane waves the corresponding magnetic field can be calculated by the electric field amplitude and the wavevector.  
Inside the simulation box, the initial electromagnetic field is set up by the "UserDefinedInitialMaxwellStates" block. The first column denotes the vector component of the field, and the second column is the input type. The third column choses the Maxwell field and the last one reads the spatial dependent formula.  
In case of "MaxwellIncidentWavesInBox = yes", the code uses the initial analytical plane wave values also inside the simulation box. Hence, there is no Maxwell field propagation with the algorithm, the corresponding values are set analytically.

	# ----- Maxwell field variables -----------------------------------------------------------------

	# Plane waves input block
	# - Column 1: input method (in this example, a a predefined Maxwell function)
	# - Column 2: electric field amplitude in the x- direction
	# - Column 3: electric field amplitude in the y- direction
	# - Column 4: electric field amplitude in the z- direction
	# - Column 5: Maxwell function name
	# - Column 6: wave type

	%UserDefinedMaxwellIncidentWaves
	plane_waves_mx_function | Exl | Eyl | E21 | "plane_waves_func_1" | plane_wave
	plane_waves_mx_function | Ex2 | Ey2 | E22 | "plane_waves_func_2" | plane_wave
	%

	# Predefined Maxwell function
	# - Column 1: Maxwell function name that corresponds to UserDefinedMaxwellIncidentWaves
	# - Column 2: envelope function of the plane wave pulse
	# - Column 3: wavevector component in x-direction
	# - Column 4: wavevector component in y-direction
	# - Column 5: wavevector component in z-direction
	# - Column 6: pulse shift in x-direction
	# - Column 7: pulse shift in y-direction
	# - Column 8: pulse shift in z-direction
	# - Column 9: pulse width
	# kxl, kx2 = wavevectors in x-direction for EM fields 1 and 2
	# ky1,ky2 = wavevectors in y-direction for EM fields 1 and 2
	# k21,k22 = wavevectors in z-direction for EM fields 1 and 2
	# psx1,psx2 = spatial shift of pulse 1 and 2 in x-direction
	# psy1,psy2 = spatial shift of pulse 1 and 2 in y-direction
	# pszl, p522 = spatial shift of pulse 1 and 2 in z-direction
	# pw1, pw2 = pulse width of pulse 1 and 2
	%MaxwellFunctions
	"plane_waves_func_1" | mxf_cosinoidal_wave | kx1 | ky1 | kz1 | psx1 | psy1 | psz1 | pw1
	"plane_waves_func_2" | mxf_cosinoidal_wave | kx2 | ky2 | kzZ | psx2 | psy2 | psz2 | pw2

	# cosinoidal pulse
	# - Column 1: field component
	# - Column 2: input type (formula or file)
	# - Column 3: what kind of Maxwell field
	# - Column 4: analytic formula of file input for the field component
	# e.g. Ex2(x,y,z) gives here the x component of the 2nd pulse's electric field 
	# e.g. By1(x,y,z) gives here the y component of the 1st pulse's magnetic field

	%UserDefinedInitialMaxwellStates
	1 | formula | electric_field | " Ex1(x,y,z) + Ex2(x,y,z) "
	2 | formula | electric_field | " Ey1(x,y,z) + Ey2(x,y,z) "
	3 | formula | electric_field | " Ez1(x,y,z) + Ez2(x,y,z) "
	1 | formula | magnetic_field | " Bx1(x,y,z) + Bx2(x,y,z) "
	2 | formula | magnetic_field | " By1(x,y,z) + By2(x,y,z) "
	3 | formula | magnetic_field | " Bz1(x,y,z) + Bz2(x,y,z) "
	%

	# If incident plane waves should be used to calculate the analytical Maxwell field inside the simulation box (no propagation, just evaluation of the plane waves inside the simulation box) 

	MaxwellIncidentWavesInBox = no

#### External current  

A external current can be switched that adds an external current contribution to the internal current. Such external current is calculated analytically by the "UserDefinedMaxwellExternalCurrent" block. The first column represents the input type. The next three columns describe the spatial distribution of the external current density. The fifth column reads the temporal frequency of the current pulse, and the last column names the envelope of the current pulse.  
The "TDFunctions" block defines the time-depending function for the current pulse. The first column corresponds to the name of the envelope function in the previous block. The second column reads the function type, the third the amplitude factor, the fourth the pulse width and the last column the time shift.

	# ----- External current options ------------------------

	# Switch on external current density
	MaxwellExternalCurrent = yes

	# External current parameters
	# - Column 1: input option for the td function which is here a predefined function
	# - Column 2: spatial distribution of the current in x- direction
	# - Column 3: spatial distribution of the current in y-direction
	# - Column 4: spatial distribution of the current in z-direction
	# - Column 5: frequency of the temoral current density pulse
	# - Column 6: name of the temporal current pulse envelope function
	# t0 = time shift of the current density signal
	# tw = width of the current density signal in time
	# jx(x,y,z) = spatial current density distribution in x-direction
	# jy(x,y,z) = spatial current density distribution in y-direction
	# jz(x,y,z) = spatial current density distribution in z-direction
	# omega = frequency of the temoral current density pulse
	# env_func = name of the temporal current pulse envelope function
	%UserDefinedMaxwellExternalCurrent
	external_current_td_function | "jx(x.y,z)" | ”jy(x,y,z)" | "jz(x,y,z)" | omega | "env_func"
	%

	# Function in time
	# - Column 1: name of the envelope function
	# - Column 2: function type
	# - Column 3: amplitude
	# - Column 4: envelope function width
	# - Column 5: shift of the function f(t-t0)

	%TDFunctions
	"env_func” | tdf_gaussian | 1.0 | tw | to
	%


#### Instantaneous field increasing  

The instantaneous field increasing modus adds a certain field value, which is constant in space to the current propagated field value. Such added constant field value is calculated by a time-depending function. The first six columns in "UserDefinedConstantSpatialMaxwellField" read the constant electromagnetic field amplitude. The last column refers to the time-depending function.  
The "TDFunctions" block defines the time-depending function.  

	# ----- Spatial constant field propagation ------------------------------------------------------

	Ez           = 0.000010000000
	By           = 0.000010000000
	pulse_width  = 500.0
	pulse_shift  = 270.0
	pulse_slope  = 100.0

	# Column 1: constant electric field component in x-direction
	# Column 2: constant electric field component in y-direction
	# Column 3: constant electric field component in z-direction
	# Column 4: constant magnetic field component in x-direction
	# Column 5: constant magnetic field component in y-direction
	# Column 6: constant magnetic field component in z-direction
	# Column 7: name of the td function

	%UserDefinedConstantSpatialMaxwellField
	0 | 0 | Ez | 0 | By | 0 | "time_function"
	%

	# Column 1: name of the td function
	# Column 2: function
	# Column 3: amplitude
	# Column 4: pulse 510
	# Column 5: pulse wid
	# Column 6: pulse shi

	%TDFunctions
	"time_function" | tdf_logistic | 1.0 | pulse_slope | pulse_width | pulse_shift
	%


### Linear Medium

Linear Media can be included in the calculation either through a simple box, or through a file describing more complex geometries.

#### Linear Medium Boxes

Linear medium boxes can be placed inside the Maxwell simulation box. They are defined by their centers, their sizes in each dimension, and their electromagnetic medium parameters. The first three columns read the center coordinate of the corresponding box, the following three columns define the box length in each direction. Columns seven and eight read the electric permittivity and magnetic permeability factors. Columns nine and ten give the electric and magnetic loss, and the last column choses the type of box edges. They can be hard edged or smoothed by a continuous function.

	# ----- Medium Box ------------------------------------------------------------------------------

	# Medium box parameters like size, epsilon, mu, electric loss and magnetic loss
	# Column 1: Length of the box in x-direction
	# Column 2: Length of the box in y-direction
	# Column 3: Length of the box in z-direction
	# Column 4: Center coordinate of the box in x-direction
	# Column 5: Center coordinate of the box in y-direction
	# Column 6: Center coordinate of the box in z-direction
	# Column 7: Electric permittivity factor (epsilon_r)
	# Column 8: Magnetic permeability factor (mu_r)
	# Column 9: Electric loss (sigma_e)
	# Column 10: Magnetic loss (sigma_m)
	# Column 11: Edge option ("edged" for sharp edges, "smooth" for using a smoothing function)

	% MaxwellMediumBox
	0.0 | 0.0 | 0.0 | 10.0 | 10.0 | 10.0 | 2.0 | 2.0 | 0.0 | 0.0 | edged
	%

	# Maxwell operation calculated via Riemann-Silberstein vector or separately by
	# the electric and magnetic field with standard Maxwell's equations to avoid errors
	# by calculating descrete gradients for sharped edged boxes


#### Complex geometries:

Medium surface file, followed by permittivity factor, electric conductivity
and magnetic conductivity, and finally type of numerical approximation used
for the derivatives at the edges.


	# ---- Linear Medium from file -----------------------------------------

	%LinearMediumFromFile
       medium_surface_file | epsilon_factor | mu_factor | sigma_e | sigma_m | edged/smooth
    %
 

For  linear  media the calculation of the Maxwell Operator acting on the RS
state  can be done directly using the Riemann-Silberstein representation or
by calculating the curl of the electric and magnetic fields.

	# ---- Medium Calculation options -------------------------------------
	
	# MaxwellMediumCalculation: Options
	# - MaxwellMediumCalculation = RS (Medium calculation via Riemann-Silberstein)
	# - MaxwellMediumCalculation = EM (Medium calculation via curl of electric field and magnetic field)
	
	MaxwellMediumCalculation = {{mymacro('RS')}}



Test for {{mymacro('test')}}.

--8<-- "includes/abberviations.md"