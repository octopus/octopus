---
title: "Kronig-Penney Model"
tags: ["Beginner", "Ground State", "Model", "Chain", "Band Structure", "User-defined Species", "Independent Particles"]
series: "Tutorial"
---


The Kronig-Penney model is a 1D system that demonstrates band gaps, which relate to the allowed energies for electrons in a material. In this tutorial we calculate the bandstructure for Kronig-Penney Model. The Kronig-Penney Model has a periodic potential of

$V(x) =
\begin{cases}
```text
      V_0 & -b<x<0 \\
      0 & 0<x<a
```
\end{cases}
$$

Where b is the width of each barrier, and a is the spacing between them. 
```text
 
```
##  Input  
The following input file will be used for the ground state calculation:

```text
 {{< Variable2 "CalculationMode" >}} = gs
 {{< Variable2 "ExtraStates" >}} = 4
 {{< Variable2 "PeriodicDimensions" >}} = 1
 {{< Variable2 "Dimensions" >}} = 1
 {{< Variable2 "TheoryLevel" >}} = independent_particles
 
  a = 5
  b = 1
  V = 3
 
 {{< Variable2 "Lsize" >}} = (a + b)/2
 
 {{< Variable2 "Spacing" >}} = .0075
 
 %{{< Variable2 "Species" >}}
  "A" | species_user_defined | potential_formula | "(x>-b)*V*(x<0)" | valence | 1
 %
 
 %{{< Variable2 "Coordinates" >}}
  "A" | 0 |
 %
 
 %{{< Variable2 "KPointsGrid" >}}
  11 |
 %
 
 %{{< Variable2 "KPointsPath" >}}
  11 |
 0.0 |
 0.5 |
 %
 
 {{< Variable2 "ConvEigenError" >}} = true
```

{{< figure src="/Kp_wavefunctions.png" width="500px" caption="Caption" >}}
The first two wavefunctions plotted alongside the potential.

##  Bandstructure   
To calculate the bandstructure simply change the {{< Variable2 "CalculationMode" >}} to unocc.

```text
 {{< Variable2 "CalculationMode" >}} = unocc
 {{< Variable2 "ExtraStates" >}} = 4
 {{< Variable2 "PeriodicDimensions" >}} = 1
 {{< Variable2 "Dimensions" >}} = 1
 {{< Variable2 "TheoryLevel" >}} = independent_particles
 
 a = 5
 b = 1
 V = 3
 
 {{< Variable2 "Lsize" >}} = (a + b)/2
 
 {{< Variable2 "Spacing" >}} = .0075
 
 %{{< Variable2 "Species" >}}
  "A" | species_user_defined | potential_formula | "(x>-b)*V*(x<0)" | valence | 1
 %
 
 %{{< Variable2 "Coordinates" >}}
  "A" | 0 |
 %
 
 %{{< Variable2 "KPointsGrid" >}}
  11 |
 %
 
 %{{< Variable2 "KPointsPath" >}}
  11 |
 0.0 |
 0.5 |
 %
 
 {{< Variable2 "ConvEigenError" >}} = true
```

{{< figure src="/Kp_bandstructure.png" width="500px" caption="The band structure for Kronig-Penney Model." >}}

To plot the bandstructure, we will use the same command from the [Periodic systems](../Periodic systems) (assuming you are using gnuplot).

```text
 plot for [col=5:5+9] 'static/bandstructure' u 1:(column(col)) w l notitle ls 1
```

Reference:
Sidebottom DL. Fundamentals of condensed matter and crystalline physics: an introduction for students of physics and materials science. New York: Cambridge University Press; 2012.








---------------------------------------------
