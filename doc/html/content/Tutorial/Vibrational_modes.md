---
title: "Vibrational modes"
tags: ["Tutorial", "Advanced", "Vibrational Modes", "Geometry Optimization", "Molecule", "Pseudopotentials", "DFT", "Forces", "Vibrations"]
#series: "Tutorial"
---


This tutorial will show you how to obtain vibrational modes in {{< octopus >}} by using Sternheimer linear response. We will use the water molecule as example. To start we will need the geometry. Copy the following in a file called {{< file "h2o.xyz" >}}:

{{< code-block >}}
 3
 Water molecule in Angstroms 
 O  0.000 0.000 0.0
 H  0.757 0.586 0.0
 H -0.757 0.586 0.0
{{< /code-block >}}

Since to obtain proper vibrational modes we need to be as close as possible to the minimum of energy, first we will optimize the geometry under the parameters we have defined for the simulation. For that we use the following input file:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = go
 
 {{< variable "UnitsOutput" >}} =  eV_Angstrom
 {{< variable "XYZCoordinates" >}} = "h2o.xyz"
 
 {{< variable "Spacing" >}} = 0.16*angstrom
 {{< variable "Radius" >}} = 4.5*angstrom
 {{< variable "BoxShape" >}} = minimum
  
 {{< variable "FilterPotentials" >}} = filter_ts
 
 {{< variable "GOMethod" >}} = fire
 
 {{< variable "ExperimentalFeatures" >}} = yes
{{< /code-block >}}


As you see we will use a smaller than usual spacing and a larger spherical box. This is required since forces require more precision than other quantities. Of course for production runs these values should be properly converged. Also to increase the precision of the forces and the total energy we will use the variable {{< variable "FilterPotentials" >}}, which specifies that a filter should be applied to the ionic potentials that removes the Fourier components that are higher than what we can represent on our grid.

Now run {{< octopus >}}, once the calculation is finished you should get the optimized geometry in {{< file "min.xyz" >}}. Modify the value of {{< variable "XYZCoordinates" >}} to point to this file.

To perform a calculation of the vibrational modes we need a ground-state calculation with a very well converged forces. So change the {{< variable "CalculationMode" >}} to {{< code "gs" >}} and run {{< octopus >}} again, by adding {{< variable "ConvForce" >}} = 1.e-7 ( the criterion to stop the self-consistency process will be given by the change in the force). Once done, check that the forces in the {{< file "static/info" >}} file are reasonably small (<0.05 eV/A).

Now we are ready for the actual calculation of vibrational modes. To do it, change {{< variable "CalculationMode" >}} to {{< code "vib_modes" >}} and run {{< octopus >}}. This calculation will take a while, since $3N_{atoms}$ response calculations are required.

Finally, the results will be written into the {{< code "vib_modes" >}} folder.

- Check the {{< file "normal_frequencies_lr" >}} file. 

{{< code-block >}}
     1    3714.66865804
     2    3612.19180564
     3    1541.89763162
     4     282.95538165
     5     189.98486286
     6     155.62791393
     7    -197.74013582
     8    -246.65205687
     9    -258.63857751
{{< /code-block >}}

What kind of information does it give? Did we find the correct equilibrium geometry? 

Try now to visualize the vibrational eigenmodes. You can use XCrySDen to open the file {{< file "normal_modes_lr.axsf" >}}. In the "Display" menu set "Forces", and in "Modify" set the factor to 1. Try to classify the modes as translations, rotations, and true internal vibrations. Which ones have the "negative" (actually imaginary) frequencies?