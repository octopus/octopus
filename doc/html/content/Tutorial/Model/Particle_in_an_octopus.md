---
title: "Particle in an octopus"
tags: ["Beginner", "Ground State", "Model", "User-defined Species", "Independent Particles", "Visualization"]
weight: 4
#series: "Tutorial"
---


{{< octopus >}} can actually run 2D systems where the shape of the simulation box is defined by what is white in an image file. Here is an example of a "particle in an octopus", in which we have a constant potential and an octopus-shaped quantum dot. To run it, you will need to have built the code with the [optional library GDLIB](https://libgd.github.io).

##  Input  

For this example we will need two files:

####  {{< file "inp" >}}  

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "FromScratch" >}} = yes
 {{< variable "Dimensions" >}} = 2
 
 %{{< variable "Species" >}}
 "null" | species_user_defined | potential_formula | "0" | valence | 1
 %
 
 %{{< variable "Coordinates" >}}
 "null" | 0 | 0
 %
 
 {{< variable "BoxShape" >}} = box_image
 {{< variable "BoxShapeImage" >}} = "gdlib.png"
 
 ff = 20
 %{{< variable "Lsize" >}}
  135 / ff | 95 / ff
 %
 
 {{< variable "TheoryLevel" >}} = independent_particles
 {{< variable "ConvEigenError" >}} = yes
{{< /code-block >}}

####  {{< file "Gdlib.png" >}}   

Make this file available in the run directory. You can download it by clicking on the image bellow. It is also available in the {{< file "PREFIX/share/octopus" >}} directory from your {{< octopus >}} installation.

{{< figure src="/images/Gdlib.png" >}}

##  Plotting  

You can obtain the wavefunction by adding this to the input file:

{{< code-block >}}
 {{< variable "Output" >}} = wfs
 {{< variable "OutputFormat" >}} = plane_z
{{< /code-block >}}

and rerunning. View it in {{< code "gnuplot" >}} with

```bash
 plot 'static/wf-st0001.z=0' u 1:2:3 linetype palette
```

or

```bash
 splot 'static/wf-st0001.z=0' u 1:2:(0):($3*500) with pm3d
```

Where does the wavefunction localize, and why?

##  Exercises  
* See how the total energy scales with the size of the system (controlled by the {{< code ff >}} parameter in the input file). How does it compare to the formula for a particle in a box?
* Look at the wavefunctions of the unoccupied states.
* Think of a serious application that would use the ability to define the simulation box by an image!

{{< tutorial-foot series="Model" prev="1D Helium" next="" >}}

