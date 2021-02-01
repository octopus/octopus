---
title: "Optical spectra from Sternheimer"
tags: ["Beginner", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "Electromagnetic Response", "Sternheimer"]
Weight: 4
#series: "Tutorial"
---


We have just seen how to calculate optical spectra in the time domain with a finite perturbation, and in a frequency-domain, linear-response matrix formulation with the Casida equation. Now we will try a third approach, which is in the frequency domain and linear response but rather than using a pseudo-eigenvalue equation as in Casida, uses a self-consistent linear equation, the Sternheimer equation. This approach is also known as density-functional perturbation theory. It has superior scaling, is more efficient for dense spectra, and is more applicable to nonlinear response. One disadvantage is that one needs to proceed one frequency point at a time, rather than getting the whole spectrum at once. We will find we can obtain equivalent results with this approach for the optical spectra as for time propagation and Casida, by calculating the polarizability and taking the imaginary part.

##  Ground state  

Before doing linear response, we need to obtain the ground state of the system, for which we can use the same input file as for [Optical spectra from Casida](../Optical spectra from Casida), but we will use a tighter numerical tolerance, which helps the Sternheimer equation to be solved more rapidly. Unlike for Casida, no unoccupied states are required. If they are present, they won't be used anyway.

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
 {{< variable "ConvRelDens" >}} = 1e-7
{{< /code-block >}}

##  Linear response  

Add the lines below to the input file (replacing the {{< variable "CalculationMode" >}} line).

The frequencies of interest must be specified, and we choose them based on the what we have seen from the Casida spectrum. The block below specifies 5 frequencies spanning the range 0 to 8 eV (below the resonances) and 9 frequencies spanning the range 10 to 12 eV (where there are peaks). If we didn't know where to look, then looking at a coarse frequency grid and then sampling more points in the region that seems to have a peak (including looking for signs of resonances in the real part of the polarizability) would be a reasonable approach. We must add a small imaginary part ($\eta$) to the frequency in order to be able to obtain the imaginary part of the response, and to avoid divergence at resonances. The resonances are broadened into Lorentzians with this width. The larger the $\eta$, the easier the SCF convergence is, but the lower the resolution of the spectrum.

To help in the numerical solution, we turn off the preconditioner (which sometimes causes trouble here), and use a linear solver that is experimental but will give convergence much faster than the default one.

{{< code-block >}}
 {{< variable "CalculationMode" >}} = em_resp
 %{{< variable "EMFreqs" >}}
 5 | 0*eV | 8*eV
 9 | 10*eV | 12*eV
 %
 {{< variable "EMEta" >}} = 0.1*eV
 
 {{< variable "Preconditioner" >}} = no
 {{< variable "LinearSolver" >}} = qmr_dotp
 {{< variable "ExperimentalFeatures" >}} = yes
{{< /code-block >}}

##  Spectrum  

This Perl script can extract the needed info out of the files in the {{< file "em_resp" >}} directory:

```text
 #!/bin/perl
 use Math::Trig;
 $c = 137.035999679;
 
 $ls = `ls -d freq*`;
 @list = ($ls =~/freq_([0-9.]*)/g);
 @freqs = sort {$a <=> $b} @list;
 
 print "# Energy (eV)    Re alpha         Im alpha    Cross-section (A^2)\n";
 format =
 @###.###      @####.#######   @####.########   @####.########
 $energy, $Re_alpha    ,  $Im_alpha,         $cross_section
 .
 
 foreach(@freqs)
 {
     if ($_ eq "0.0000")
     {
         $Im_alpha = 0;
     }
     else
     {
         $crossfile = `cat freq_$_/cross_section`;
         @crossbits = split(' ', $crossfile);
         $energy = $crossbits[26];
         $cross_section = $crossbits[27]; # isotropic average
         $Im_alpha = $c * $cross_section / (4 * pi * $energy);
     }
 
     $alphafile = `cat freq_$_/alpha`;
     @alphabits = split(' ', $alphafile);
     $Re_alpha = $alphabits[15];
 
     write;
 }
```

Our result is

```text
 # Energy (eV)    Re alpha         Im alpha    Cross-section (A^2)
    0.000          2.64122000       0.00000000       0.00000000
    2.000          2.69965700       0.00041822       0.00007670
    4.000          2.89886900       0.00101989       0.00037410
    6.000          3.35097500       0.00234919       0.00129254
    8.000          4.70072400       0.01004299       0.00736763
   10.000          3.58189900       0.04386600       0.04022566
   10.250          4.14888500       0.11959579       0.11241259
   10.500          4.86710100       0.05001950       0.04816193
   10.750          6.78277900       0.07330539       0.07226359
   11.000         10.49267200       0.26261161       0.26489990
   11.250          3.33778500       1.06435158       1.09802650
   11.500         -0.47501400       0.20728472       0.21859505
   11.750          4.28837500       0.28012324       0.30182986
   12.000         -1.06249200       0.21176755       0.23303215
```

##  See also  

{{< tutorial "Sternheimer linear response" "Sternheimer linear response" >}}

{{< tutorial-foot series="response" prev="Optical spectra from Casida" next="Triplet excitations" >}}

