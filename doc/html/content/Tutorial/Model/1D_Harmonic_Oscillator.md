---
title: "1D Harmonic Oscillator"
tags:: ["Beginner", "Ground State", "Unoccupied", "Model", "User-defined Species", "Independent Particles"]
weight: 1
#series: "Tutorial"
---


As a first example we use the standard textbook harmonic oscillator in one dimension and fill it with two non-interacting electrons. 

##  Input  

The first thing to do is to tell {{< octopus >}} what we want it to do. Write the following lines and save the file as {{< file "inp" >}}.

{{< code-block >}}
 {{< variable "FromScratch" >}} = yes
 {{< variable "CalculationMode" >}} = gs
 
 {{< variable "Dimensions" >}} = 1 
 {{< variable "TheoryLevel" >}} = independent_particles
 
 {{< variable "Radius" >}} = 10
 {{< variable "Spacing" >}} = 0.1
 
 %{{< variable "Species" >}} 
   "HO" | species_user_defined | potential_formula | "0.5*x^2" | valence | 2
 %
 
 %{{< variable "Coordinates" >}}
   "HO" | 0
 %
{{< /code-block >}}

Most of these input variables should already be familiar. Here is a more detailed explanation for some of the values:
* {{< code-inline >}}{{< variable "Dimensions" >}} = 1{{< /code-inline >}}: This means that the code should run for 1-dimensional systems. Other options are 2, or 3 (the default). You can actually run in 4D too if you have compiled with the configure flag <tt>--max-dim=4</tt>.

* {{< code-inline >}}{{< variable "TheoryLevel" >}} = independent_particles{{< /code-inline >}}: We tell {{< octopus >}} to treat the electrons as non-interacting.

* {{< code-inline >}}{{< variable "Radius" >}} = 10.0{{< /code-inline >}}: The radius of the 1D "sphere," ''i.e.'' a line; therefore domain extends from -10 to +10 bohr.

* {{< code-inline >}}%{{< variable "Species" >}}{{< /code-inline >}}: The species name is "HO", then the potential formula is given, and finally the number of valence electrons. See {{< manual "Input file" "Input file" >}} for a description of what kind of expressions can be given for the potential formula.

* {{< code-inline >}}%{{< variable "Coordinates" >}}{{< /code-inline >}}: add the external potential defined for species "HO" in the position (0,0,0). The coordinates used in the potential formula are relative to this point.

##  Output  
Now one can execute this file by running {{< octopus >}}. Here are a few things worthy of a closer look in the standard output.

First one finds the listing of the species:

{{< code-block >}}
 ****************************** {{< variable "Species" >}} *******************************
 Species "HO" is an user-defined potential.
    Potential = 0.5*x^2
 Number of orbitals:      5
 **********************************************************************
{{< /code-block >}}

The potential is $V(x) = 0.5 x^2$, and 5 Hermite-polynomial orbitals are available for LCAO (the number is based on the valence). The theory level is as we requested:

{{< code-block >}}
 **************************** Theory Level ****************************
 Input: [{{< variable "TheoryLevel" >}} = independent_particles]
 **********************************************************************
{{< /code-block >}}

The electrons are treated as "non-interacting", which means that the Hartree and exchange-correlation terms are not included. This is usually appropriate for model systems, in particular because the standard XC approximations we use for actual electrons are not correct for "effective electrons" with a different mass.

{{< code-block >}}
 Input: [{{< variable "MixField" >}} = none] (what to mix during SCF cycles)
{{< /code-block >}}

Since we are using independent particles (and only one electron) there is no need to use a mixing scheme to accelerate the SCF convergence. 

{{< code-block >}}
 Input: [{{< variable "LCAOStart" >}} = lcao_none]
 Info: Unnormalized total charge =      2.000000
 Info: Renormalized total charge =      2.000000
 Info: Setting up Hamiltonian.
 Orthogonalizing wavefunctions.
{{< /code-block >}}

Starting from scratch means that {{< octopus >}} generates a starting density from the sum of atomic densities. This is then renormalized to integrate to the total number of electrons present in the system. For atomic systems, the default is to find a first estimate for the wave-functions using a {{< manual "Calculations/Ground_State#LCAO" "linear combination of atomic orbitals" >}} (LCAO) technique, using the atomic wavefunctions from the pseudopotentials. However, we do not necessarily have such corresponding wavefunctions for a user-defined potential, so LCAO is turned off by default here. The starting orbitals will then be random but orthogonal.

Now the self-consistent cycle starts and you should see the following output: 

{{< code-block >}}
Info: Starting SCF iteration.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER #    1 ************************
 etot  =  1.24484521E+00 abs_ev   =  9.49E+01 rel_ev   =  7.62E+01
 ediff =        1.24E+00 abs_dens =  3.27E+00 rel_dens =  1.63E+00
Matrix vector products:     27
Converged eigenvectors:      0

#  State  Eigenvalue [H]  Occupation    Error
      1        0.622423    2.000000   (6.3E+00)

Elapsed time for SCF step     1:          0.00
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0
{{< /code-block >}}

and after a few iterations it converges:

{{< code-block >}}
*********************** SCF CYCLE ITER #   22 ************************
 etot  =  1.00000001E+00 abs_ev   =  1.95E-08 rel_ev   =  1.95E-08
 ediff =       -1.95E-08 abs_dens =  1.61E-05 rel_dens =  8.07E-06
Matrix vector products:     27
Converged eigenvectors:      0

#  State  Eigenvalue [H]  Occupation    Error
      1        0.500000    2.000000   (1.4E-03)

Elapsed time for SCF step    22:          0.00
**********************************************************************


             Info: Writing states. 2018/08/08 at 14:13:39


        Info: Finished writing states. 2018/08/08 at 14:13:39

Info: SCF converged in   22 iterations
{{< /code-block >}}

Note that, as the electrons are non-interacting, there is actually no self-consistency needed. Nevertheless, {{< octopus >}} takes several SCF iterations to converge. This is because it takes more than one SCF iteration for the eigensolver to converge the wave-functions and eigenvalues. To improve this we can try having a better initial guess for the wave-functions by turning the LCAO on. You can do this by setting <tt>{{< variable "LCAOStart" >}} = lcao_states</tt> -- compare how many iterations and matrix-vector products (in total) are required now. Why? (Hint: what are Hermite polynomials?)

##  Exercises  

The first thing to do is to look at the {{< file "static/info" >}} file. Look at the eigenvalue and total energy. Compare to your expectation from the analytic solution to this problem!

We can also play a little bit with the input file and add some other features. For example, what about the other eigenstates of this system? To obtain them we can calculate some unoccupied states. This can be done by changing the <tt>CalculationMode</tt> to

{{< code-block >}}
 {{< variable "CalculationMode" >}} = unocc
{{< /code-block >}}

in the input file. In this mode {{< octopus >}} will not only give you the occupied states (which contribute to the density) but also the unoccupied ones. Set the number of states with an extra line

{{< code-block >}}
 {{< variable "ExtraStates" >}} = 10
{{< /code-block >}}

that will calculate 10 empty states. A thing to note here is that {{< octopus >}} will need the density for this calculation. (Actually for non-interacting electrons the density is not needed for anything, but since {{< octopus >}} is designed for interacting electrons it will try to read the density anyways.) So if you have already performed a static calculation (<tt>CalculationMode = gs</tt>) it will just use this result. 

Compared to the ground-state calculation we also have to change the convergence criterion. The unoccupied states do not contribute to the density so they might not (and actually will not) converge properly if we use the density for the convergence check. Therefore, {{< octopus >}} now checks whether all calculated states converge separately by just looking at the biggest error in the wave-functions. From the calculation you get another file in the {{< file "static" >}} directory called {{< file "eigenvalues" >}}. The {{< file "info" >}} file will only contain the information about the ground-state; all eigenvalues and occupation numbers will be in the {{< file "eigenvalues" >}} file.

If we also want to plot, say, the wave-function, at the end of the calculation, we have to tell {{< octopus >}} to give us this wave-function and how it should do this. We just include

{{< code-block >}}
 {{< variable "Output" >}} = wfs
 {{< variable "OutputFormat" >}} = axis_x
{{< /code-block >}}

The first line tells {{< octopus >}} to give out the wave-functions and the second line says it should do so along the x-axis. We can also select the wave-functions we would like as output, for example the first and second, the fourth and the sixth. (If you don't specify anything {{< octopus >}} will give them all.)

{{< code-block >}}
 {{< variable "OutputWfsNumber" >}} = "1-2,4,6"
{{< /code-block >}}

{{< octopus >}} will store the wave-functions in the same folder {{< file "static" >}} where the {{< file "info" >}} file is, under a meaningful name. They are stored as pairs of the ''x''-coordinate and the value of the wave-function at that position ''x''. One can easily plot them with gnuplot or a similar program.

It is also possible to extract a couple of other things from {{< octopus >}} like the density or the Kohn-Sham potential.

{{< tutorial-foot series="Model"  prev="" next="Particle in a box" >}}







---------------------------------------------
