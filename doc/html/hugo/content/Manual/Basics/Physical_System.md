---
title: "Physical System"
series: "Manual"
---


The first thing that {{< octopus >}} has to know is the physical system you want to treat. To do this you have specify a group of species and their positions. 

##  Dimensions  

{{< octopus >}} can work in a space with 1, 2 or 3 dimensions. You can select the dimension of your system with the {{< Variable2 "Dimensions" >}} variable.

##  Species  

An {{< octopus >}} species is very generic and can be a nucleus (represented by pseudopotentials or by the full Coulomb potential), a jellium sphere or even a user-defined potential. The information regarding the atomic species goes into the {{< Variable2 "Species" >}} block.

###  Pseudopotentials  
The many-electron Schroedinger equation can be greatly simpliﬁed if electrons are
divided in two groups: valence electrons and inner core electrons. The electrons
in the inner shells are strongly bound and do not play a signiﬁcant role in the
chemical binding of atoms, thus forming with the nucleus an inert core. Binding
properties are almost completely due to the valence electrons, especially in metals
and semiconductors. This separation implies that inner electrons can be ignored,
reducing the atom to an inert ionic core that interacts with the valence electrons.
This suggests the use of an eﬀective interaction, a pseudopotential, that gives an
approximation to the potential felt by the valence electrons due to the nucleus and
the core electrons. This can signiﬁcantly reduce the number of electrons that have
to be dealt with. Moreover, the pseudo wave functions of these valence electrons
are much smoother in the core region than the true valence wave functions, thus
reducing the computational burden of the calculations.

Modern pseudopotentials are obtained by inverting the free atom Schroedinger equation for a given reference electronic conﬁguration, and forcing the pseudo wave
functions to coincide with the true valence wave functions beyond a certain cutoﬀ
distance. The pseudo wave functions are also forced to have the same norm as the
true valence wave functions, and the energy pseudo eigenvalues are matched to the
true valence eigenvalues. Diﬀerent methods of obtaining a pseudo eigenfunction
that satisﬁes all these requirements lead to diﬀerent non-local, angular momentum
dependent pseudopotentials. Some widely used pseudopotentials are the Troullier
and Martins potentials, the Hamann potentials, the Vanderbilt potentials and the Hartwigsen-Goedecker-Hutter potentials. The default potentials
used by {{< octopus >}} are of the Troullier and Martins type, although you can also opt
for the HGH potentials.

{{< octopus >}} comes with a package of pseudopotentials and the parameters needed to use them. If you want to have a look you can find them under {{< inst_file "share/octopus/pseudopotentials" >}}. These pseudopotentials serve to define many species that you might wish to use in your coordinates block, e.g. a helium atom, "He".  If you are happy to use these predefined pseudopotentials, you do not need to write a species block.

However it is also possible to define new species to use in your coordinates block by adding a {{< Variable2 "Species" >}} block.  You can check the documentation of that variable for the specific syntax.  With this block you may specify the format of your file that contains a pseudopotential and parameters such as the atomic number, or you may define an algebraic expression for the potential with the user-defined potential. A user defined potential should be finite everywhere in the region where your calculation runs.

If you want to search other repositories on the web or create your own pseudopotentials, check the [[Pseudopotentials|pseudopotentials]] page. Save the pseudopotential file in the same directory as the inp file and specify its format with the species block.

###  All-Electron Nucleus  

The potential of this species is the full Coulomb potential 

$$
V\\left( \\vec{r}\\right)=\\frac1{\\left|\\vec{R}\_i-\\vec{r}\\right|}\\ .
$$

The main problem to represent this potential is the discontinuity over $\vec{R}_i$. To overcome this problem we do the following: 

* First we assume that atoms are located over the closest grid point.
* Then we calculate the charge density associated with the nucleus: a delta distribution with the value of the charge at this point and zero elsewhere.
* Now we solve the Poisson equation for this density.

In this way we get a potential that is the best representation of the Coulomb potential for our grid (we will discuss about grids later) and is continuous in $\vec{R}_i$ (the value is the average of the potential over the volume associated with the grid point).

The main problem is that the requirement of having atoms over grid points is quite strong: it is only possible for a few systems with simple geometries and you can't move the atoms. Also the Coulomb potential is very hard, which means you will need a very small spacing, and as you have to consider both core and valence electrons, this species is only suitable for atoms or very small molecules.

###  User Defined  

It is also possible to define an external, user-defined potential in the input file. All {{< Manual "Input_file-Mathematical_expressions " "functions" >}} accepted by the parser can be used. Besides that, one can use the symbols $x$, $y$, $z$, and $r$. In this way it is trivial to calculate model systems, like harmonic oscillators, quantum dots, etc.

##  Coordinates  

For each instance of a species (even for user-defined potentials), you have to specify its position inside the simulation box. To do this you can use the {{< Variable2 "Coordinates" >}} block which describes the positions inside of the input file or one of the {{< Variable2 "XYZCoordinates" >}} or {{< Variable2 "PDBCoordinates" >}} variables, that specify an external file, in {{< name "xyz" >}} or {{< name "PDB" >}} format respectively, from where the coordinates will be read.

Before using a geometry with {{< octopus >}} we recommend that you center it. For this you can use the {{< Manual "External utilities:oct-center-geom" "oct-center-geom" >}} utility.

##  Velocities  

If you are going to do ion dynamics you may want to have an initial velocity for the particles. You have several choices for doing this:

* Don't put anything in the input file; particles will have zero initial velocity.
* Give them a random velocity according to a temperature (in degrees Kelvin) given by the {{< Variable2 "RandomVelocityTemp" >}} variable.
* Explicitly give the initial velocity for each particle, either through the {{< Variable2 "Velocities" >}} block or from a pseudo-xyz file detailed by the variable  {{< Variable2 "XYZVelocities" >}}.

##  Number of Electrons  

Each species adds enough electrons to make the system neutral. If you want to add or remove electrons you can specify the total charge of your system with the {{< Variable2 "ExcessCharge" >}} variable (a negative charge implies to add electrons).

{{< Manual_foot prev="Manual:Units" next="Manual:Hamiltonian" >}}
---------------------------------------------
