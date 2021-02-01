---
title: "Triplet excitations"
tags: ["Advanced", "Time-dependent", "Casida", "Molecule", "Pseudopotentials", "DFT", "Triplet Excitations", "oct-propagation_spectrum", "oct-casida_spectrum"]
weight: 5
#series: "Tutorial"
---


In this tutorial, we will calculate triplet excitations for methane with time-propagation and Casida methods.

##  Time-propagation  

###  Ground-state  
We begin with a spin-polarized calculation of the ground-state, as before but with {{< code-inline >}}{{< variable "SpinComponents" >}} = spin_polarized{{< /code-inline >}} specified now.

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = eV_angstrom
 
 {{< variable "Radius" >}} = 6.5*angstrom
 {{< variable "Spacing" >}} = 0.24*angstrom
 
 CH = 1.097*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
 
 {{< variable "SpinComponents" >}} = spin_polarized
{{< /code-block >}}

You can verify that the results are identical in detail to the non-spin-polarized calculation since this is a non-magnetic system.

###  Time-propagation  

Next, we perform the time-propagation using the following input file:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = td
 {{< variable "UnitsOutput" >}} = eV_angstrom
 
 {{< variable "Radius" >}} = 6.5*angstrom
 {{< variable "Spacing" >}} = 0.24*angstrom
 
 CH = 1.097*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
 
 {{< variable "SpinComponents" >}} = spin_polarized
 
 {{< variable "TDPropagator" >}} = aetrs
 {{< variable "TDTimeStep" >}} = 0.004/eV
 {{< variable "TDMaxSteps" >}} = 2500  - ~ 10.0/{{< variable "TDTimeStep" >}}
 
 {{< variable "TDDeltaStrength" >}} = 0.01/angstrom
 {{< variable "TDPolarizationDirection" >}} = 1
 {{< variable "TDDeltaStrengthMode" >}} = kick_spin
{{< /code-block >}}


Besides the {{< variable "SpinComponents" >}} variable, the main difference is the type of perturbation that is applied to the system. By setting {{< code-inline >}}
{{< variable "TDDeltaStrengthMode" >}} = kick_spin{{< /code-inline >}}, the kick will have opposite sign for up and down states. Whereas the ordinary kick ({{< code kick_density >}}) yields the response to a homogeneous electric field, ''i.e.'' the electric dipole response, this kick yields the response to a homogeneous magnetic field, ''i.e.'' the magnetic dipole response. Note however that only the spin degree of freedom is coupling to the field; a different calculation would be required to obtain the orbital part of the response. Only singlet excited states contribute to the spectrum with {{< code kick_density >}}, and only triplet excited states contribute with {{< code kick_spin >}}. We will see below how to use symmetry to obtain both at once with {{< code kick_spin_and_density >}}.

###  Spectrum  

When the propagation completes, run the {{< file "oct-propagation_spectrum" >}} utility to obtain the spectrum. This is how the {{< file "cross_section_vector" >}} should look like.

{{< code-block >}}
# nspin         2
# kick mode    1
# kick strength    0.005291772086
# pol(1)           1.000000000000    0.000000000000    0.000000000000
# pol(2)           0.000000000000    1.000000000000    0.000000000000
# pol(3)           0.000000000000    0.000000000000    1.000000000000
# direction    1
# Equiv. axes  0
# wprime           0.000000000000    0.000000000000    1.000000000000
# kick time        0.000000000000
#%
# Number of time steps =     2500
# PropagationSpectrumDampMode   =    2
# PropagationSpectrumDampFactor =     4.0817
# PropagationSpectrumStartTime  =     0.0000
# PropagationSpectrumEndTime    =    10.0000
# PropagationSpectrumMaxEnergy  =    20.0000
# PropagationSpectrumEnergyStep =     0.0100
#%
#       Energy        sigma(1, nspin=1)   sigma(2, nspin=1)   sigma(3, nspin=1)   sigma(1, nspin=2)   sigma(2, nspin=2)   sigma(3, nspin=2)   StrengthFunction(1) StrengthFunction(2)
#        [eV]               [A^2]               [A^2]               [A^2]               [A^2]               [A^2]               [A^2]               [1/eV]              [1/eV]       
      0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00
      0.10000000E-01     -0.24458945E-08      0.16678045E-11      0.19693772E-11      0.53059881E-08      0.20990576E-11      0.66109455E-12     -0.22283826E-08      0.48341299E-08
      0.20000000E-01     -0.97344740E-08      0.66632996E-11      0.78681086E-11      0.21157433E-07      0.83862450E-11      0.26412279E-11     -0.88687932E-08      0.19275916E-07
{{< /code-block >}}

{{< figure src="/images/Singlet_triplet_spectrum_CH4.png" width="500px" caption="Comparison of absorption spectrum of methane calculated with time-propagation for singlets and triplets." >}}

You can see that there are now separate columns for cross-section and strength function for each spin. The physically meaningful strength function for the magnetic excitation is given by {{< code "StrengthFunction(1)" >}} - {{< code "StrengthFunction(2)" >}} (since the kick was opposite for the two spins). (If we had obtained {{< file "cross_section_tensor" >}}, then the trace in the second column would be the appropriate cross-section to consider.) We can plot and compare to the singlet results obtained before. You can see how this looks on the right. The first triplet transition is found at 9.05 eV, slightly lower energy than the lowest singlet transition. 

If you are interested, you can also repeat the calculation for {{< code-inline >}}{{< variable "TDDeltaStrengthMode" >}} = kick_density{{< /code-inline >}} (the default) and confirm that the result is the same as for the non-spin-polarized calculation.

###  Using symmetries of non-magnetic systems  

As said before, methane is a non-magnetic system, that is, the up and down densities are the same and the magnetization density is zero everywhere:

$$
 \\rho^{\\uparrow}(\\mathbf r) = \\rho^{\\downarrow}(\\mathbf r)
$$

Note that it is not enough that the total magnetic moment of the system is zero as the previous condition does not hold for anti-ferromagnetic systems. The symmetry in the spin-densities can actually be exploited in order to obtain both the singlet and triplet spectra with a single calculation. This is done by using a special perturbation that is only applied to the spin up states.[^footnote-1] 
To use this perturbation, we need to set {{< code-inline >}}{{< variable "TDDeltaStrengthMode" >}} = kick_spin_and_density{{</ code-inline >}}. If you repeat the time-propagation with this kick, you should obtain a different {{< file "cross_section_vector" >}} file containing both the singlet and triplet spectra. The singlets are given by {{< code "StrengthFunction(1)" >}} + {{< code "StrengthFunction(2)" >}}, while the triplets are given by {{< code "StrengthFunction(1)" >}} - {{< code "StrengthFunction(2)" >}}.

##  Casida equation  

The calculation of triplets with the {{< code casida>}} mode for spin-polarized systems is currently not implement in {{< octopus >}}. Nevertheless, just like for the time-propagation, we can take advantage that for non-magnetic systems the two spin are equivalent. In this case it allow us to calculate triplets without the need for a spin-polarized run. The effective kernels in these cases are:

$f_{\rm Hxc}^{\rm singlet} \left[ \rho \right] = f^{\uparrow}_{\rm Hxc} \left[ \rho ^{\uparrow} \right] + f^{\uparrow}_{\rm Hxc} \left[ \rho^{\downarrow} \right] = f_{\rm H} \left[ \rho \right] + f^{\uparrow}_{\rm xc} \left[ \rho ^{\uparrow} \right] + f^{\uparrow}_{\rm xc} \left[ \rho^{\downarrow} \right]$

$f_{\rm Hxc}^{\rm triplet} \left[ \rho \right] = f^{\uparrow}_{\rm Hxc} \left[ \rho ^{\uparrow} \right] - f^{\uparrow}_{\rm Hxc} \left[ \rho^{\downarrow} \right] = f^{\uparrow}_{\rm xc} \left[ \rho ^{\uparrow} \right] - f^{\uparrow}_{\rm xc} \left[ \rho^{\downarrow} \right]$

Therefore, we start by doing a ground-state and unoccupied states runs exactly as was done in the {{< tutorial "Response/Optical spectra from Casida" "Optical spectra from Casida tutorial" >}}. Then, do a Casida run with the following input file:

{{< figure src="/images/Singlet_triplet_spectrum_Casida_CH4.png" width="500px" caption="Comparison of absorption spectrum of methane calculated with the Casida equation for singlets and triplets." >}}

{{< code-block >}}
 {{< variable "CalculationMode" >}} = casida
 {{< variable "UnitsOutput" >}} = eV_angstrom
 
 {{< variable "Radius" >}} = 4.*angstrom
 {{< variable "Spacing" >}} = 0.18*angstrom
 
 CH = 1.097*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0 
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
 
 {{< variable "ExtraStates" >}} = 10
 
 {{< variable "CasidaCalcTriplet" >}} = yes
{{< /code-block >}}

The only difference with respect to the calculation of the singlets, is the {{< variable "CasidaCalcTriplet" >}} input variable.

Once the Casida calculation is finished, run {{< file "oct-casida_spectrum" >}}, plot the results, and compare to the singlet calculation. On the right you can see how this plot should look like. How do the singlet and triplet energy levels compare? Can you explain a general relation between them? How does the run-time compare between singlet and triplet, and why?

##  Comparison  

As for the singlet spectrum, we can compare the time-propagation and Casida results. What is the main difference, and what is the reason for it?

{{< figure src="/images/Triplet_casida_td_CH4.png" width="500px" caption="Comparison of triplet absorption spectrum of methane calculated with time-propagation and with the Casida equation." >}}


[^footnote-1]: {{< article title="On the use of Neumann's principle for the calculation of the polarizability tensor of nanostructures" authors="M.J.T. Oliveira, A. Castro, M.A.L. Marques, and A. Rubio" journal="J. Nanoscience and Nanotechnology" volume="8" pages="1-7" year="2008" doi="10.1166/jnn.2008.142" arxiv="0710.2624v1" >}}


{{< tutorial-foot series="response" prev="Optical spectra from Sternheimer" next="Use of symmetries in optical spectra from time-propagation" >}}
