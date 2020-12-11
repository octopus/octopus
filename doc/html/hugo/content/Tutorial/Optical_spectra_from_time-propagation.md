---
title: "Optical spectra from time-propagation"
tags: ["Beginner", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-propagation_spectrum"]
series: "Tutorial"
---


In this tutorial we will learn how to obtain the absorption spectrum of a molecule from the explicit solution of the time-dependent Kohn-Sham equations. We choose as a test case methane (CH<sub>4</sub>).

##  Ground state  

Before starting out time-dependent simulations, we need to obtain the ground state of the system. For this we use basically the same {{< file "inp" >}} file as in the [Total energy convergence](../Total energy convergence) tutorial:

```text
 {{< Variable2 "CalculationMode" >}} = gs
 {{< Variable2 "UnitsOutput" >}} = eV_angstrom
 
 {{< Variable2 "Radius" >}} = 3.5*angstrom
 {{< Variable2 "Spacing" >}} = 0.18*angstrom
 
 CH = 1.097*angstrom
 %{{< Variable2 "Coordinates" >}}
   "C" |           0 |          0 |           0
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
```

After running {{< octopus >}}, we will have the Kohn-Sham wave-functions of the ground-state in the directory {{< file "restart/gs" >}}. As we are going to propagate these wave-functions, they have to be well converged. It is important not only to converge the energy (that is relatively easy to converge), but the density. The default <tt>{{< Variable2 "ConvRelDens" >}} = 1e-05</tt> is usually enough, though.

##  Time-dependent run  

To calculate absorption, we excite the system with an infinitesimal electric-field pulse, and then propagate the time-dependent Kohn-Sham equations for a certain time ''T''. The spectrum can then be evaluated from the time-dependent dipole moment.

####  Input  

This is how the input file should look for the time propagation. It is similar to the one from the [Time-dependent propagation](../Time-dependent propagation) tutorial.

```text
 {{< Variable2 "CalculationMode" >}} = td
 {{< Variable2 "FromScratch" >}} = yes
 {{< Variable2 "UnitsOutput" >}} = eV_angstrom
 
 {{< Variable2 "Radius" >}} = 3.5*angstrom
 {{< Variable2 "Spacing" >}} = 0.18*angstrom
 
 CH = 1.097*angstrom
 %{{< Variable2 "Coordinates" >}}
  "C" |           0 |          0 |           0
  "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
  "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
  "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
  "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
       
 {{< Variable2 "TDPropagator" >}} = aetrs
 {{< Variable2 "TDTimeStep" >}} = 0.0023/eV
 {{< Variable2 "TDMaxSteps" >}} = 4350  - ~ 10.0/TDTimeStep
 
 {{< Variable2 "TDDeltaStrength" >}} = 0.01/angstrom
 {{< Variable2 "TDPolarizationDirection" >}} = 1
```

Besides changing the {{< Variable2 "CalculationMode" >}} to <tt>td</tt>, we have added {{< Variable2 "FromScratch" >}} = yes. This will be useful if you decide to run the propagation for other polarization directions (see bellow). For the time-evolution we use again the Approximately Enforced Time-Reversal Symmetry (aetrs) propagator. The time-step is chosen such that the propagation remains numerically stable. You should have learned how to set it up in the [Time-dependent propagation](../Time-dependent propagation) tutorial. Finally, we set the number of time steps with the variable {{< Variable2 "TDMaxSteps" >}}. To have a maximum propagation time of 10 $\hbar/{\rm eV}$ we will need around 4350 iterations.

We have also introduced two new input variables to define our perturbation:

* {{< Variable2 "TDDeltaStrength" >}}: this is the strength of the perturbation. This number should be small to keep the response linear, but should be sufficiently large to avoid numerical problems.

* {{< Variable2 "TDPolarizationDirection" >}}: this variable sets the polarization of our perturbation to be on the first axis (''x'').

Note that you will be calculating the singlet dipole spectrum. You can also obtain the triplet by using {{< Variable2 "TDDeltaStrengthMode" >}}, and other multipole responses by using {{< Variable2 "TDKickFunction" >}}. For details on triplet calculations see the [Triplet excitations](../Triplet excitations) tutorial.

You can now start {{< octopus >}} and go for a quick coffee (this should take a few minutes depending on your machine). Propagations are slow, but the good news is that they scale very well with the size of the system. This means that even if methane is very slow, a molecule with 200 atoms can still be calculated without big problems.

####  Output  

The output should be very similar to the one from the [Time-dependent propagation](../Time-dependent propagation) tutorial. The main difference is the information about the perturbation:
```text

Info: Applying delta kick: k =    0.005292
Info: kick function: dipole.
Info: Delta kick mode: Density mode
```
</pre>

You should also get a {{< file "td.general/multipoles" >}} file, which contains the necessary information to calculate the spectrum. The beginning of this file should look like this:
```text

--------------------------------------------------------------------------------
- HEADER
- nspin         1
- lmax          1
- kick mode     0
- kick strength    0.005291772086
- pol(1)           1.000000000000    0.000000000000    0.000000000000
- pol(2)           0.000000000000    1.000000000000    0.000000000000
- pol(3)           0.000000000000    0.000000000000    1.000000000000
- direction     1
- Equiv. axes   0
- wprime           0.000000000000    0.000000000000    1.000000000000
- kick time        0.000000000000
- Iter             t          Electronic charge(1)       <x>(1)              <y>(1)              <z>(1)       
-[Iter n.]     [hbar/eV]           Electrons              [A]                 [A]                 [A]         
--------------------------------------------------------------------------------
       0  0.000000000000e+00  8.000000000000e+00 -5.105069841261e-11 -2.458728679159e-11 -6.191777548908e-11
       1  2.300000000000e-03  7.999999999898e+00 -1.441204189916e-03 -2.454483489431e-11 -6.178536414960e-11
       2  4.600000000000e-03  7.999999999794e+00 -2.849596623682e-03 -2.441311503813e-11 -6.139129111463e-11
       3  6.900000000000e-03  7.999999999692e+00 -4.212595751069e-03 -2.420204397245e-11 -6.075589410104e-11
```
</pre>
Note how the dipole along the ''x'' direction (forth column) changes in response to the perturbation. 

##  Optical spectra  

In order to obtain the spectrum for a general system one would need to perform three time-dependent runs, each for a perturbation along a different Cartesian direction (''x'', ''y'', and ''z''). In practice this is what you would need to do:
- Set the direction of the perturbation along ''x'' ({{< Variable2 "TDPolarizationDirection" >}} = 1),
- Run the time-propagation,
- Rename the {{< file "td.general/multipoles" >}} file to {{< file "td.general/multipoles.1" >}},
- Repeat the above step for directions ''y'' ({{< Variable2 "TDPolarizationDirection" >}} = 2) and ''z'' ({{< Variable2 "TDPolarizationDirection" >}} = 3) to obtain files {{< file "td.general/multipoles.2" >}} and {{< file "td.general/multipoles.3" >}}.

Nevertheless, the $T_d$ symmetry of methane means that the response is identical for all directions and the absorption spectrum for ''x''-polarization will in fact be equivalent to the spectrum averaged over the three directions. You can perform the calculations for the ''y'' and ''z'' directions if you need to convince yourself that they are indeed equivalent, but you this will not be necessary to complete this tutorial.

Note that {{< octopus >}} can actually use the knowledge of the symmetries of the system when calculating the spectrum. However, this is fairly complicated to understand for the purposes of this introductory tutorial, so it will be covered in the [Use of symmetries in optical spectra from time-propagation](../Use of symmetries in optical spectra from time-propagation) tutorial.

###  The {{< file "oct-propagation_spectrum" >}} utility  

{{< octopus >}} provides an utility called {{< file "oct-propagation_spectrum" >}} to process the {{< file "multipoles" >}} files and obtain the spectrum. 

####  Input  
This utility requires little information from the input file, as most of what it needs is provided in the header of the {{< file "multipoles" >}} files. If you want you can reuse the same input file as for the time-propagation run, but the following input file will also work:

```text
 {{< Variable2 "UnitsOutput" >}} = eV_angstrom
```

Now run the utility. Don't worry about the warnings generated, we know what we are doing!

####  Output  

This is what you should get:
```text

************************** Spectrum Options **************************
Input: [PropagationSpectrumType = AbsorptionSpectrum]
Input: [SpectrumMethod = fourier]
Input: [PropagationSpectrumDampMode = polynomial]
Input: [PropagationSpectrumTransform = sine]
Input: [PropagationSpectrumStartTime = 0.000 hbar/eV]
Input: [PropagationSpectrumEndTime = -.3675E-01 hbar/eV]
Input: [PropagationSpectrumEnergyStep = 0.1000E-01 eV]
Input: [PropagationSpectrumMaxEnergy = 20.00 eV]
Input: [PropagationSpectrumDampFactor = -27.21 hbar/eV^-1]
Input: [PropagationSpectrumSigmaDiagonalization = no]
**********************************************************************

File "multipoles" found. This will be the only file to be processed.
(If more than one file is to be used, the files should be called
"multipoles.1", "multipoles.2", etc.)


** Warning:
**   The file "multipoles" tells me that the system has no usable symmetry.
**   However, I am only using this file; cannot calculate the full tensor.
**   A file "XXXX_vector" will be generated instead.


Octopus emitted 1 warning.
```
</pre>

You will notice that the file {{< file "cross_section_vector" >}} is created.  If you have all the required information (either a symmetric molecule and one multipole file, or a less symmetric molecule and multiple multipole files), {{< file "oct-propagation_spectrum" >}} will generate the whole polarizability tensor, {{< file "cross_section_tensor" >}}.  If not enough multipole files are available, it can only generate the response to the particular polarization you chose in you time-dependent run, {{< file "cross_section_vector" >}}. 

###  Cross-section vector   

Let us look first at {{< file "cross_section_vector" >}}:
```text

- nspin         1
- kick mode    0
- kick strength    0.005291772086
- pol(1)           1.000000000000    0.000000000000    0.000000000000
- pol(2)           0.000000000000    1.000000000000    0.000000000000
- pol(3)           0.000000000000    0.000000000000    1.000000000000
- direction    1
- Equiv. axes  0
- wprime           0.000000000000    0.000000000000    1.000000000000
- kick time        0.000000000000
-%
```
</pre>

The beginning of the file just repeats some information concerning the run.

```text

- Number of time steps =     4350
- PropagationSpectrumDampMode   =    2
- PropagationSpectrumDampFactor =   -27.2114
- PropagationSpectrumStartTime  =     0.0000
- PropagationSpectrumEndTime    =    10.0050
- PropagationSpectrumMaxEnergy  =    20.0000
- PropagationSpectrumEnergyStep =     0.0100
-%
```
</pre>

Now comes a summary of the variables used to calculate the spectra. Note that you can change all these settings just by adding these variables to the input file before running {{< file "oct-propagation_spectrum" >}}.  Of special importance are perhaps {{< Variable2 "PropagationSpectrumMaxEnergy" >}} that sets the maximum energy that will be calculated, and {{< Variable2 "PropagationSpectrumEnergyStep" >}} that determines how many points your spectrum will contain. To have smoother curves, you should reduce this last variable.

```text

- Electronic sum rule       =         3.682690
- Static polarizability (from sum rule) =         2.064454 A
-%
```
</pre>

Now comes some information from the sum rules. The first is just the ''f''-sum rule, which should yield the number of active electrons in the calculations. We have 8 valence electrons (4 from carbon and 4 from hydrogen), but the sum rule gives a number that is much smaller than 8! The reason is that we are just summing our spectrum up to 20 eV (see {{< Variable2 "PropagationSpectrumMaxEnergy" >}}), but for the sum rule to be fulfilled, we should go to infinity. Of course infinity is a bit too large, but increasing 20 eV to a somewhat larger number will improve dramatically the ''f''-sum rule. The second number is the static polarizability also calculated from the sum rule. Again, do not forget to converge this number with {{< Variable2 "PropagationSpectrumMaxEnergy" >}}.

```text

-       Energy        sigma(1, nspin=1)   sigma(2, nspin=1)   sigma(3, nspin=1)   StrengthFunction(1)
-        [eV]               [A^2]               [A^2]               [A^2]               [1/eV]
      0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00
      0.10000000E-01     -0.42061649E-08     -0.71263367E-14      0.21595354E-13     -0.38321132E-08
      0.20000000E-01     -0.16805683E-07     -0.28462613E-13      0.86270383E-13     -0.15311164E-07
...
```
</pre>

Finally comes the spectrum. The first column is the energy (frequency), the next three columns are a row of the cross-section tensor, and the last one is the strength function for this run.

The dynamic polarizability is related to optical absorption cross-section via $\sigma \left( \omega \right) = \frac{4 \pi \omega}{c} \mathrm{Im}\ \alpha \left( \omega \right) $ in atomic units, or more generally $4 \pi \omega \tilde{\alpha}\ \mathrm{Im}\ \alpha \left( \omega \right) $ (where $\tilde{\alpha}$ is the fine-structure constant) or $\frac{\omega e^2}{\epsilon_0 c} \mathrm{Im}\ \alpha \left( \omega \right) $. The cross-section is related to the strength function by $S \left( \omega \right) = \frac{mc}{2 \pi^2 \hbar^2} \sigma \left( \omega \right)$.

###  Cross-section tensor   

If all three directions were done, we would have four files: {{< file "cross_section_vector.1" >}}, {{< file "cross_section_vector.2" >}}, {{< file "cross_section_vector.3" >}}, and {{< file "cross_section_tensor" >}}. The latter would be similar to:

```text

- nspin         1
- kick mode    0
- kick strength    0.005291772086
- pol(1)           1.000000000000    0.000000000000    0.000000000000
- pol(2)           0.000000000000    1.000000000000    0.000000000000
- pol(3)           0.000000000000    0.000000000000    1.000000000000
- direction    3
- Equiv. axes  0
- wprime           0.000000000000    0.000000000000    1.000000000000
- kick time        0.000000000000
-       Energy         (1/3)*Tr[sigma]    Anisotropy[sigma]      sigma(1,1,1)        sigma(1,2,1)        sigma(1,3,1)        ...
-        [eV]               [A^2]               [A^2]               [A^2]               [A^2]               [A^2]            ...
      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00     ...
      0.10000000E-01     -0.42062861E-08      0.16626408E-12     -0.42061649E-08      0.18742668E-12      0.18909806E-12     ...
      0.20000000E-01     -0.16806167E-07      0.66452380E-12     -0.16805683E-07      0.74882279E-12      0.75548381E-12     ...
...
```
</pre>

{{< figure src="/Absorption_spectrum_CH4.png" width="500px" caption="Absorption spectrum of CH<sub>4</sub>" >}}

The columns are now the energy, the average absorption coefficient (Tr $\sigma/3$), the anisotropy, and 
then all the 9 components of the tensor. The third number is the spin component, just 1 here since it is unpolarized. The anisotropy is defined as

$$
 (\\Delta \\sigma)^2 = \\frac{1}{3} \\{ 3{\\rm Tr}(\\sigma^2) - \[{\\rm Tr}(\\sigma)\]^2 \\}
$$

The reason for this definition is that it is identically equal to:

$$
 (\\Delta \\sigma)^2 = (1/3) \[ (\\sigma\_1-\\sigma\_2)^2 + (\\sigma\_1-\\sigma\_3)^2 + (\\sigma\_2-\\sigma\_3)^2 \]\\,
$$

where $\{\sigma_1, \sigma_2, \sigma_3\}\,$ are the eigenvalues of $\sigma\,$. An "isotropic" tensor is characterized by having three equal eigenvalues, which leads to zero anisotropy. The more different that the eigenvalues are, the larger the anisotropy is.

If you now plot the absorption spectrum (column 5 vs 1 in {{< file "cross_section_vector" >}}, but you can use {{< file "cross_section_tensor" >}} for this exercise in case you did propagate in all three directions), you should obtain the plot shown on the right. Of course, you should now try to converge this spectrum with respect to the calculation parameters. In particular:

* Increasing the total propagation time will reduce the width of the peaks. In fact, the width of the peaks in this methods is absolutely artificial, and is inversely proportional to the total propagation time. Do not forget, however, that the area under the peaks has a physical meaning: it is the oscillator strength of the transition.
* A convergence study with respect to the spacing and box size might be necessary. This is covered in the next tutorial.

Some questions to think about:
* What is the equation for the peak width in terms of the propagation time? How does the observed width compare to your expectation?
* What is the highest energy excitation obtainable with the calculation parameters here? Which is the key one controlling the maximum energy?
* Why does the spectrum go below zero? Is this physical? What calculation parameters might change this?

{{Tutorial_foot|series=Optical response|prev=|next=Convergence of the optical spectra}}








---------------------------------------------
