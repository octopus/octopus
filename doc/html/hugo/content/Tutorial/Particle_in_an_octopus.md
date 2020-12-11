---
title: "Particle in an octopus"
tags: ["Beginner", "Ground State", "Model", "User-defined Species", "Independent Particles", "Visualization"]
series: "Tutorial"
---


{{< octopus >}} can actually run 2D systems where the shape of the simulation box is defined by what is white in an image file. Here is an example of a "particle in an octopus", in which we have a constant potential and an octopus-shaped quantum dot. To run it, you will need to have built the code with the [https://libgd.github.io optional library GDLIB].

##  Input  

For this example we will need two files:

####  {{< file "inp" >}}  

```text
 {{< Variable2 "CalculationMode" >}} = gs
 {{< Variable2 "FromScratch" >}} = yes
 {{< Variable2 "Dimensions" >}} = 2
 
 %{{< Variable2 "Species" >}}
 "null" | species_user_defined | potential_formula | "0" | valence | 1
 %
 
 %{{< Variable2 "Coordinates" >}}
 "null" | 0 | 0
 %
 
 {{< Variable2 "BoxShape" >}} = box_image
 {{< Variable2 "BoxShapeImage" >}} = "gdlib.png"
 
 ff = 20
 %{{< Variable2 "Lsize" >}}
  135 / ff | 95 / ff
 %
 
 {{< Variable2 "TheoryLevel" >}} = independent_particles
 {{< Variable2 "ConvEigenError" >}} = yes
```

####  {{< file "gdlib.png" >}}   

Make this file available in the run directory. You can download it by clicking on the image bellow. It is also available in the {{< file "PREFIX/share/octopus" >}} directory from your {{< octopus >}} installation.

[[Image:gdlib.png]]

##  Plotting  

You can obtain the wavefunction by adding this to the input file:

```text
 {{< Variable2 "Output" >}} = wfs
 {{< Variable2 "OutputFormat" >}} = plane_z
```

and rerunning. View it in {{< code "gnuplot" >}} with

```text
 plot 'static/wf-st0001.z=0' u 1:2:3 linetype palette
```

or

```text
 splot 'static/wf-st0001.z=0' u 1:2:(0):($3*500) with pm3d
```

Where does the wavefunction localize, and why?

##  Exercises  
* See how the total energy scales with the size of the system (controlled by the <tt>ff</tt> parameter in the input file). How does it compare to the formula for a particle in a box?
* Look at the wavefunctions of the unoccupied states.
* Think of a serious application that would use the ability to define the simulation box by an image!

{{Tutorial_foot|series=Model systems|prev=1D Helium|next=}}







---------------------------------------------
