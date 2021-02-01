---
title: "Centering a geometry"
tags: ["Basic", "Molecule", "oct-center-geom"]
#series: "Tutorial"
weight: 5
---


Before running an {{< octopus >}} calculation of a molecule, it is always a good idea to run the {{< manual "External utilities:oct-center-geom"  "oct-center-geom" >}} utility, which will translate the center of mass to the origin, and align the molecule, by default so its main axis is along the ''x''-axis. Doing this is often helpful for visualization purposes, and making clear the symmetry of the system, and also it will help to construct the simulation box efficiently in the code. The current implementation in {{< octopus >}} constructs a parallelepiped containing the simulation box and the origin, and it will be much larger than necessary if the system is not centered and aligned. For periodic systems, these considerations are not relevant (at least in the periodic directions).

##  Input  
For this example we need two files.

####  {{< file "inp" >}}  
We need only a very simple input file, specifying the coordinates file.

{{< code-block >}}
 {{< variable "XYZCoordinates" >}} = "tAB.xyz"
{{< /code-block >}}

####  {{< file "tAB.xyz" >}}  

We will use this coordinates file, for the molecule [http://en.wikipedia.org/wiki/Azobenzene ''trans''-azobenzene], which has the interesting property of being able to switch between ''trans'' and ''cis'' isomers by absorption of light.

{{< code-block >}}
 24
 
 C   -0.520939   -1.29036    -2.47763
 C   -0.580306    0.055376   -2.10547
 C    0.530312    0.670804   -1.51903
 C    1.71302    -0.064392   -1.30151
 C    1.76264    -1.41322    -1.67813
 C    0.649275   -2.02387    -2.26419
 H   -1.38157    -1.76465    -2.93131
 H   -1.48735     0.622185   -2.27142
 H    2.66553    -1.98883    -1.51628
 H    0.693958   -3.06606    -2.55287
 H    0.467127    1.71395    -1.23683
 N    2.85908     0.52877    -0.708609
 N    2.86708     1.72545    -0.355435
 C    4.01307     2.3187      0.237467
 C    3.96326     3.66755     0.614027
 C    5.0765      4.27839     1.20009
 C    6.2468      3.54504     1.41361
 C    6.30637     2.19928     1.04153
 C    5.19586     1.58368     0.455065
 H    5.25919     0.540517    0.17292
 H    3.06029     4.24301     0.4521
 H    5.03165     5.32058     1.48872
 H    7.10734     4.01949     1.86731
 H    7.21349     1.63261     1.20754
{{< /code-block >}}

##  Centering the geometry   

[[Manual:External utilities:oct-center-geom|{{< code "oct-center-geom" >}}]] is a serial utility, and it will be found in the {{< code "bin" >}} directory after installation of {{< octopus >}}.

When you run the utility, you should obtain a new coordinates file {{< code "adjusted.xyz" >}} for use in your calculations with {{< octopus >}}. The output is not especially interesting, except for perhaps the symmetries and moment of inertia tensor, which come at the end of the calculation:
 
{{< code-block >}}
***************************** Symmetries *****************************
Symmetry elements : (sigma)
Symmetry group    : Cs
**********************************************************************

Using main axis        1.000000       0.000000       0.000000
Input: [AxisType = inertia]
Center of mass [b] =        5.410288       2.130154      -1.005363
Moment of inertia tensor [amu*b^2]
      6542962.320615     -4064048.884840     -3486905.049193
     -4064048.884840      8371776.892023     -2789525.719087
     -3486905.049193     -2789525.719087     10730451.796958
Isotropic average      8548397.003199
Eigenvalues:            1199576.606657          11623018.898148          12822595.504791
Found primary   axis        0.707107       0.565686       0.424264
Found secondary axis       -0.624695       0.780869      -0.000000
{{< /code-block >}}

You can now visualize the original and new coordinates files with {{< code "xcrysden" >}} or your favorite visualization program (''e.g.'' Jmol, Avogadro, VMD, Vesta, etc.), and see what transformations the utility performed. You can see more options about how to align the system at {{< variable "AxisType" >}} and {{< variable "MainAxis" >}}.

{{< tutorial-foot series="basics" prev="Visualization" next="Periodic systems" >}}

