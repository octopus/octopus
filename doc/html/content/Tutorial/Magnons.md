---
title: "Magnons"
tags: ["Tutorial", "Advanced", "Bulk", "Magnons in real time"]
#series: "Tutorial"
---


The objective of this tutorial is to give a basic idea of how it is possible to compute magnons and transverse spin susceptibilities from a real-time calculation using {{< octopus >}}.

Important remark: this tutorial requires to be able to run calculations which have a non-negligible numerical cost and cannot be run locally.

{{< notice note >}}
This tutorial is still under construction.
{{< /notice >}}


###  Magnon kick in supercell  

As a first approach, we will investigate how to obtain transverse spin susceptibilities using supercells adapted to host a spin spiral of a specific momentum.

####  Ground-state input file  

As a prototypical example for investigating spin-waves, we will consider bulk iron, which is a ferromagnetic material. We start here by computing its ground-state by constructing a supercell corresponding to the cubic cell of iron, doubled along one direction.
For this, we use the following input file:

{{< code-block >}}
  {{< variable "CalculationMode" >}} = gs
  {{< variable "PeriodicDimensions" >}} = 3 
  {{< variable "BoxShape" >}} = parallelepiped
  {{< variable "ExperimentalFeatures" >}} = yes
  {{< variable "PseudopotentialSet" >}}=pseudodojo_lda
  a = 2.867*Angstrom
  %{{< variable "LatticeParameters" >}}
    a | a | 2*a
    90| 90| 90
  %
  %{{< variable "ReducedCoordinates" >}}
   "Fe" | 0.0 | 0.0 | 0.0
   "Fe" | 1/2 | 1/2 | 1/4
   "Fe" | 0.0 | 0.0 | 1/2
   "Fe" | 1/2 | 1/2 | 3/4
  %
  {{< variable "Spacing" >}} = 0.35
  %{{< variable "KPointsGrid" >}}
  4 | 4 | 2
  %
  {{< variable "Smearing" >}} = 0.1*eV
  {{< variable "SmearingFunction" >}} = fermi_dirac
  {{< variable "LCAOStart" >}} = lcao_states
  {{< variable "SpinComponents" >}} = spinors
  {{< variable "GuessMagnetDensity" >}} = user_defined
  %{{< variable "AtomsMagnetDirection" >}}
   0.0 | 0.0 | 4.0
   0.0 | 0.0 | 4.0
   0.0 | 0.0 | 4.0
   0.0 | 0.0 | 4.0
  %
  {{< variable "EigenSolver" >}} = rmmdiis
  {{< variable "ConvRelDens" >}} = 1e-7
  {{< variable "ExtraStates" >}} = 20
{{< /code-block >}}

Running this input file will give the ground state of bulk iron in a supercell. 
It is important to note here that we are using Pauli spinors to represent the wavefunctions. This is needed because we will later on kick the system with a spiral in spin space, that cannot be represented by collinear spins.

Note that the present input file is not converged in the number of k-points and should not be considered for production runs.
Similarly, the smearing of the occupation is set here to a high temperature of 0.1 eV, which is also not realistic for practical applications. This is done here because the number of k-point is too low in our example to properly sample the energy window close to the Fermi energy. 


####  Time-dependent run  

Now that we constructed the supercell, we can investigate a spin spiral. The momemtum q  corresponding to this cell has a value of q=(0,0,2 pi/a/2), where a is lattice parameter of iron, which we took as 2.867 angstrom in our example.

Let us now look at how to specify a magnon kick in {{< octopus >}}.
There are three different points that need to be specified:
* The strength of the kick. This determines how strongly we perturb the system.
* The momentum of the kick. This is the momentum of the spin spiral we are imposing to the system's spins.
* The easy axis of the material. This is usually the z direction, but the code allows for defining an arbitrary direction. This direction is used to determine the transverse magnetization from the total magnetization computed in Cartesian coordinates. 

First of all, we need to set that we will use a magnon kick:

{{< code-block >}}
 {{< variable TDDeltaStrengthMode >}} = kick_magnon
{{< /code-block >}}

Then we specify the strength of the perturbation:

{{< code-block >}}
 {{< variable TDDeltaStrength >}} = 0.01
{{< /code-block >}}

The strength of the perturbation should be converged such that the results are guarantied to be valid within linear response. The means that changing the strength of the kick should lead to the same transverse spin susceptibility than the original calculation, as the extracted susceptibility (see next section) does not depend on the kick strength within linear response.
The momentum of spin spiral is set using the block

{{< code-block >}}
  %{{< variable TDMomentumTransfer >}}
   0 | 0 | 2*pi/a/2
  %
{{< /code-block >}}

Finally, the easy axis of the material is defined by

{{< code-block >}}
  %{{< variable TDEasyAxis >}}
    0 | 0 | 1
  %
{{< /code-block >}}

The full input file should read as

{{< code-block >}}
  {{< variable "CalculationMode" >}} = gs
  {{< variable "PeriodicDimensions" >}} = 3 
  {{< variable "BoxShape" >}} = parallelepiped
  {{< variable "ExperimentalFeatures" >}} = yes
  {{< variable "PseudopotentialSet" >}}=pseudodojo_lda
  a = 2.867*Angstrom
  %{{< variable "LatticeParameters" >}}
    a | a | 2*a
    90| 90| 90
  %
  %{{< variable "ReducedCoordinates" >}}
   "Fe" | 0.0 | 0.0 | 0.0
   "Fe" | 1/2 | 1/2 | 1/4
   "Fe" | 0.0 | 0.0 | 1/2
   "Fe" | 1/2 | 1/2 | 3/4
  %
  {{< variable "Spacing" >}} = 0.35
  %{{< variable "KPointsGrid" >}}
  4 | 4 | 2
  %
  {{< variable "RestartFixedOccupations" >}} = yes
  {{< variable "SpinComponents" >}} = spinors

  {{< variable "ExtraStates" >}} = 8

  {{< variable "RestartWriteInterval" >}} = 5000
  {{< variable "TDDeltaStrength" >}} = 0.01
  {{< variable "TDDeltaStrengthMode" >}} = kick_magnon
  %{{< variable "TDMomentumTransfer" >}}
   0 | 0 | 2*pi/a/2 
  %
  %{{< variable "TDEasyAxis" >}}
   0 | 0 | 1
  %
  {{< variable "TDTimeStep" >}} = 0.075
  {{< variable "TDPropagator" >}} = aetrs
  {{< variable "TDExponentialMethod" >}} = lanczos
  {{< variable "TDExpOrder" >}} = 16
  {{< variable "TDPropagationTime" >}} = 1800
  {{< variable "TDOutput" >}} = total_magnetization + energy
{{< /code-block >}}

The parameters for time-dependent runs are already explained in other tutorials and are not further detailed here. 
Here the time-propagation is performed only up to 1800 atomic units. This frequency resolution for this time propagation is probably too large for most applications, but we use this value here to maintain the tutorial feasible within a reasonable amount of time.

Looking at the total energy (td.general/energy) of the system, one realizes that it starts to increase, which is not correct. This is due to the poor convergence parameters used for the calculation. In practical calculation, it is strongly advised to check that the energy is conserved throughout the complete simulation.

####  Computing the transverse spin susceptibility  

In order to extract the transverse spin susceptibilities from the total magnetization, one needs to use the utility oct-spin_susceptibility.
This utility produces files names td.general/spin_susceptibility_qXXX, where XXX is the index of the q-vector. There is typically only one file produced, except in the multi-q kick mode, see below for some example.

###  Magnon kick using the generalized Bloch theorem  

If the Hamiltonian does not include any off-diagonal terms in spin space, it is possible to use the so-called generalized Bloch theorem (GBT).
Thanks to the GBT, it is possible to investigate spin spirals using only the primitive cell of the system. This comes at the cost of modifying the periodic boundary conditions by some twisted boundary conditions in which a different phase is applied to each components of the Pauli spinors.

In order to investigate it, we first need to prepare the ground-state of the system in its primitive cell.
For this, we modify the above input file to define only the primitive cell

{{< code-block >}}
  {{< variable "CalculationMode" >}} = gs
  {{< variable "PeriodicDimensions" >}} = 3 
  {{< variable "BoxShape" >}} = parallelepiped
  {{< variable "ExperimentalFeatures" >}} = yes
  {{< variable "PseudopotentialSet" >}}=pseudodojo_lda
  a = 2.867*Angstrom
  %{{< variable "LatticeParameters" >}}
    a | a | a
  %
  %{{< variable "LatticeVectors" >}}
   -0.5 | 0.5 | 0.5
    0.5 |-0.5 | 0.5
    0.5 | 0.5 |-0.5
  %
  %{{< variable "ReducedCoordinates" >}}
   "Fe" | 0.0 | 0.0 | 0.0
  %
  {{< variable "Spacing" >}} = 0.35
  %{{< variable "KPointsGrid" >}}
  4 | 4 | 4
  %
  {{< variable "Smearing" >}} = 0.1*eV
  {{< variable "SmearingFunction" >}} = fermi_dirac
  {{< variable "LCAOStart" >}} = lcao_states
  {{< variable "SpinComponents" >}} = spinors
  {{< variable "GuessMagnetDensity" >}} = user_defined
  %{{< variable "AtomsMagnetDirection" >}}
   0.0 | 0.0 | 4.0
  %
  {{< variable "EigenSolver" >}} = rmmdiis
  {{< variable "ConvRelDens" >}} = 1e-7
  {{< variable "ExtraStates" >}} = 10
{{< /code-block >}}

Once the ground state is converged, we need to run the time-dependent calculation.
This is done using the same variable as used before, only adding the line
{{< code-block >}}
 {{< variable "SpiralBoundaryCondition" >}} = yes 
{{< /code-block >}}

This variable set the use of the GBT and allow to do spin-wave calculations in primitive cell.

