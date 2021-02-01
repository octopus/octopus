## Tutorial for Maxwell propagation and its coupling with TDDFT in Octopus

### Authors: Rene Jestädt, Franco Bonafé and Heiko Appel

#### Max Planck Institute for the Structure and Dynamics of Matter - Theory Department

<br />

This tutorial is based on the following publication:

[Light-matter interactions within the Ehrenfest–Maxwell–Pauli–Kohn–Sham framework: fundamentals, implementation, and nano-optical applications](https://doi.org/10.1080/00018732.2019.1695875) René Jestädt, Michael Ruggenthaler, Micael J. T. Oliveira, Angel Rubio, and Heiko Appel

---

[Slides of Heiko Appel's talk](https://theory.mpsd.mpg.de/talks/heiko.appel/2020-01-21-Uni-Jena)

---

Using the Riemann-Silberstein representation for the Maxwell's equations, the corresponding time-evolution of Maxwell fields is implemented as quantum mechanical like propagation into the TDFT code Octopus.  

The program can be run in a free Maxwell propagation mode where the electromagnetic field is calculated inside a simulation box. This box can include a linear medium distribution. An initial electromagnetic field can be set up into the box as well as a external current density or incoming plane waves at the box boundaries. The simulation box boundaries can be selected as zero boundaries, perfectly electric conductor (PEC mirror), perfectly magnetic conductor (PMC) or absorbing boundaries with different methods.
  
Using also all features of a free Maxwell propagation, Octopus can solve fully coupled Maxwell-Kohn-Sham systems by propagating both systems on two separated grids. In this case, the Kohn-Sham system gives an internal current density that influences the Maxwell field propagation, and in turn the electromagnetic field is part of the Kohn-Sham Hamiltonian. Both coupled systems are propagated self-consistently, but the code can also be used for only forward or only backward coupling.
  
In the following we introduce the different kind of input options for the free Maxwell propagation and the fully coupled Maxwell-Kohn-Sham propagation, and different examples to show the various features.


### Simulation box and relevant input file variables  

* [Simulation box](simulationbox.md)
* [Input file variables](maxwellinputfile.md)




### Examples of electromagnetic field propagation without coupling to TDDFT

* [Cosinoidal plane wave in vacuum](./run01.html)
* [Interference of two cosinoidal plane waves](./run02.html)
* [Cosinoidal plane wave hitting a linear medium box with and without absorpbing boundaries](./run03.html)
* [Gaussian-shaped spatial external current distribution density passed by a Gaussian temporal current pulse](./run04.html)



### Coupled Maxwell-TDDFT propagation

* [Benzene ground state and dynamics coupled to an external EM pulse](./benzene_mx_matt.html)

<br />


To check what will be the input syntax once the Maxwell propagation branch is merged into the current version of Octopus and how to obtain and compile the code, click on [this link](./multisystem.html).

 
