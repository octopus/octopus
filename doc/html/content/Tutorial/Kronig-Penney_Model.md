---
title: "Kronig-Penney Model"
tags: ["Beginner", "Ground State", "Model", "Chain", "Band Structure", "User-defined Species", "Independent Particles"]
#series: "Tutorial"
---


The Kronig-Penney model is a 1D system that demonstrates band gaps, which relate to the allowed energies for electrons in a material. In this tutorial we calculate the bandstructure for Kronig-Penney Model. The Kronig-Penney Model has a periodic potential of

$$
V(x) =
\begin{cases}
      V_0 & -b<x<0 \cr
      0 & 0<x<a
\end{cases}
$$

Where b is the width of each barrier, and a is the spacing between them. 
##  Input  
The following input file will be used for the ground state calculation:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "ExtraStates" >}} = 4
 {{< variable "PeriodicDimensions" >}} = 1
 {{< variable "Dimensions" >}} = 1
 {{< variable "TheoryLevel" >}} = independent_particles
 
  a = 5
  b = 1
  V = 3
 
 {{< variable "Lsize" >}} = (a + b)/2
 
 {{< variable "Spacing" >}} = .0075
 
 %{{< variable "Species" >}}
  "A" | species_user_defined | potential_formula | "(x>-b)*V*(x<0)" | valence | 1
 %
 
 %{{< variable "Coordinates" >}}
  "A" | 0 |
 %
 
 %{{< variable "KPointsGrid" >}}
  11 |
 %
 
 %{{< variable "KPointsPath" >}}
  11 |
 0.0 |
 0.5 |
 %
 
 {{< variable "ConvEigenError" >}} = true
{{< /code-block >}}

{{< figure src="/images/Kp_wavefunctions.png" width="500px" caption="The first two wavefunctions plotted alongside the potential." >}}


##  Bandstructure   

To calculate the bandstructure simply change the {{< variable "CalculationMode" >}} to unocc.

{{< code-block >}}
 {{< variable "CalculationMode" >}} = unocc
 {{< variable "ExtraStates" >}} = 4
 {{< variable "PeriodicDimensions" >}} = 1
 {{< variable "Dimensions" >}} = 1
 {{< variable "TheoryLevel" >}} = independent_particles
 
 a = 5
 b = 1
 V = 3
 
 {{< variable "Lsize" >}} = (a + b)/2
 
 {{< variable "Spacing" >}} = .0075
 
 %{{< variable "Species" >}}
  "A" | species_user_defined | potential_formula | "(x>-b)*V*(x<0)" | valence | 1
 %
 
 %{{< variable "Coordinates" >}}
  "A" | 0 |
 %
 
 %{{< variable "KPointsGrid" >}}
  11 |
 %
 
 %{{< variable "KPointsPath" >}}
  11 |
 0.0 |
 0.5 |
 %
 
 {{< variable "ConvEigenError" >}} = true
{{< /code-block >}}

{{< figure src="/images/Kp_bandstructure.png" width="500px" caption="The band structure for Kronig-Penney Model." >}}

To plot the bandstructure, we will use the same command from the {{< tutorial "Basics:Periodic systems" "Periodic systems" >}} (assuming you are using gnuplot).

```bash
 plot for [col=5:5+9] 'static/bandstructure' u 1:(column(col)) w l notitle ls 1
```

Reference:
Sidebottom DL. Fundamentals of condensed matter and crystalline physics: an introduction for students of physics and materials science. New York: Cambridge University Press; 2012.

