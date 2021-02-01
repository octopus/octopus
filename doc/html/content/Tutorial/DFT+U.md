---
title: "DFT+U"
tags: ["Tutorial", "Advanced", "Ground State", "Bulk", "Pseudopotentials", "DFT+U"]
#series: "Tutorial"
---


The objective of this tutorial is to give a basic idea of how DFT+U in {{< octopus >}} works.

###  Input  

As a prototypical example for DFT+U, we will consider bulk NiO in its anti-ferromagnetic configuration. We will neglect the small lattice distortion and consider only its cubic cell.

{{< code-block >}}
  {{< variable "CalculationMode" >}} = gs
  {{< variable "PeriodicDimensions" >}} = 3 
  {{< variable "BoxShape" >}} = parallelepiped
  {{< variable "ExperimentalFeatures" >}} = yes
  {{< variable "PseudopotentialSet" >}}=hscv_pbe
  a = 7.8809
  %{{< variable "LatticeParameters" >}}
    a | a | a
  %
  %{{< variable "LatticeVectors" >}}
   0.0 | 1/2 | 1/2
   1/2 | 0.0 | 1/2
   1.0 | 1.0 | 0.0
  %
  
  %{{< variable "Species" >}}
  "Ni" | species_pseudo | hubbard_l | 2 | hubbard_u | 5.0*eV
  %
  
  {{< variable "DFTULevel" >}} = dft_u_empirical
  
  %{{< variable "ReducedCoordinates" >}}
   "Ni" | 0.0 | 0.0 | 0.0
   "Ni" | 0.0 | 0.0 | 0.5
   "O"  | 0.5 | 0.5 | 0.25
   "O"  | 0.5 | 0.5 | 0.75
  %
  {{< variable "Spacing" >}} = 0.5
  %{{< variable "KPointsGrid" >}}
  2 | 2 | 2
  %
  {{< variable "ParDomains" >}} = no
  {{< variable "ParKPoints" >}} = auto
  
  {{< variable "SpinComponents" >}} = polarized
  {{< variable "GuessMagnetDensity" >}} = user_defined
  %{{< variable "AtomsMagnetDirection" >}}
   8.0
  -8.0
   0.0
   0.0
  %
  
  {{< variable "OutputLDA_U" >}} = occ_matrices
{{< /code-block >}}

As we are interested by the antiferromagnetic order, the primitive cell is doubled along the last lattice vector.
To help the convergence, an initial guess should be added, by adding to the input file the variables {{< variable "GuessMagnetDensity" >}} and {{< variable "AtomsMagnetDirection" >}}.

In order to perform a calculation with a U of 5eV on the 3d orbitals (corresponding to the quantum number l=2),
we define a block Species, where hubbard_l specifies the orbitals (l=0 for s orbitals, l=1 for p orbitals, ...) and hubbard_u is used to set the value of the effective Hubbard U.

In order to activate the DFT+U part, one finally needs to specify the level of DFT+U used . This is done using the variable {{< variable "DFTULevel" >}}. At the moment there is three possible options for this variable, which correspond to no +U correction (dft_u_none), an empirical correction (dft_u_empirical) or the ab initio U correction based on the ACBN0 functional[^footnote-1] (dft_u_acbn0).

Some specific outputs can then be added, such as the density matrix of the selected localized subspaces. 

###  Output  

##  References  







---------------------------------------------
[^footnote-1]: {{< article title="Reformulation of $\mathrm{DFT}+U$ as a Pseudohybrid Hubbard Density Functional for Accelerated Materials Discovery" authors="Agapito, Luis A. and Curtarolo, Stefano and Buongiorno Nardelli, Marco" journal="Phys. Rev. X" volume="5" pages="011006" year="2015" doi="10.1103/PhysRevX.5.011006" >}}

