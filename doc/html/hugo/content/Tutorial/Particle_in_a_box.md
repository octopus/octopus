---
title: "Particle in a box"
tags: ["Beginner", "Ground State", "Model", "User-defined Species", "Independent Particles"]
series: "Tutorial"
---


In the previous tutorial, we considered applying a user-defined potential. What if we wanted to do the classic quantum problem of a particle in a box, ''i.e.'' an infinite square well?

$V(x) =
\begin{cases}
0, & \tfrac{L}{2} < x <\tfrac{L}{2}\\
\infty, & \text{otherwise}
\end{cases}
$$

There is no meaningful way to set an infinite value of the potential for a numerical calculation. However, we can instead use the boundary conditions to set up this problem. In the locations where the potential is infinite, the wavefunctions must be zero. Therefore, it is equivalent to solve for an electron in the potential above in all space, or to solve for an electron just in the domain $x \in [-\tfrac{L}{2}, \tfrac{L}{2}]$ with zero boundary conditions on the edges. In the following input file, we can accomplish this by setting the "radius" to $\tfrac{L}{2}$, for the default box shape of "sphere" which means a line in 1D.

##  Input  
As usual, we will need to write an input file describing the system we want to calculate:

```text
 {{< Variable2 "FromScratch" >}} = yes
 {{< Variable2 "CalculationMode" >}} = gs
 
 {{< Variable2 "Dimensions" >}} = 1 
 {{< Variable2 "TheoryLevel" >}} = independent_particles
 
 {{< Variable2 "Radius" >}} = 5
 {{< Variable2 "Spacing" >}} = 0.01
 
 %{{< Variable2 "Species" >}}
   "null" | species_user_defined | potential_formula | "0" | valence | 1
 %
 
 %{{< Variable2 "Coordinates" >}}
   "null" | 0
 %
 
 {{< Variable2 "LCAOStart" >}} = lcao_states
```

Run this input file and look at the ground-state energy and the eigenvalue of the single state.

##  Exercises  

* Calculate unoccupied states and check that they obey the expected relation.
* Vary the box size and check that the energy has the correct dependence.
* Plot the wavefunctions and compare to your expectation.
* Set up a calculation of a ''finite'' square well and compare results to the infinite one as a function of potential step. (Hint: along with the usual arithmetic operations, you may also use logical operators to create a piecewise-defined expression. See {{< Manual "Input file" "Input file" >}}).
* Try a 2D or 3D infinite square well.

{{Tutorial_foot|series=Model systems|prev=1D Harmonic Oscillator|next=1D Helium}}






---------------------------------------------
