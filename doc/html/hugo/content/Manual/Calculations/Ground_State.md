---
title: "Ground State"
series: "Manual"
weight: 1
description: "Converge the ground state of a system"
---


The ground-state electronic density in a Kohn-Sham (KS)-based electronic-structure code such as {{< octopus >}} is obtained after a self-consistent process that attempts to solve the KS equations. 

##  Kohn-Sham Ground State  

In essence, the problem is the following: at a given iteration step, one departs from an approximate solution – some KS eigenfunctions $\psi^{inp}_j$, eigenvalues $\epsilon^{inp}_j$ and density $\rho^{inp}$, which determines a KS “input” Hamiltonian. By diagonalizing this Hamiltonian, one obtains the corresponding “output” eigenfunctions, eigenvalues, and density. This density (or, alternatively, the corresponding Kohn-Sham potential) is then used to build a new input Hamiltonian, that will be diagonalized in the next iteration step. This cycle is considered to be closed, and the solution achieved, when the input and output are similar enough that some convergence criterion is fulfilled. In our case, we have allowed for four different criteria, to be defined below. The self-consistent procedure will stop either when the first of the convergence criterions is fulfilled, or when a maximum number of iterations has been performed.

###  Mixing  

The output density (or potential) of a given iteration is not used directly to construct the Kohn-Sham potential for the following iteration. Instead, it is "mixed" with some densities (or potentials) of previous iteration steps. The manner in which this mixing is produced is determined by the variables {{< Variable2 "MixingScheme" >}}, {{< Variable2 "Mixing" >}}, {{< Variable2 "MixField" >}} and {{< Variable2 "MixNumberSteps" >}}.

###  Convergence  
After each iteration {{< octopus >}} checks whether some convergence criterion is met.  One criterion is that the error in the electron density be smaller than some threshold.  Of course, the true electron density is not known, so this "error" is really the change in the density since the last iteration:
$\epsilon = \int {\rm d}^3r |\rho^{out}(\mathbf r) -\rho^{inp}(\mathbf r)|.$ We call this criterion {{< Variable2 "ConvAbsDens" >}}.

However, since the density is proportional to the number of electrons $N$, this absolute criterion is not very transferable between different system sizes. Therefore the default criterion to use is {{< Variable2 "ConvRelDens" >}}, in which the relative density, $ {\rho \over N}$, is used rather than the absolute density.  $\epsilon = {1\over N} \int {\rm d}^3r |\rho^{out}(\mathbf r) -\rho^{inp}(\mathbf r)|.$ By default we use a value of <tt>1e-5</tt>.

The other convergence variables {{< Variable2 "ConvAbsDens" >}}, {{< Variable2 "ConvAbsEv" >}}, and {{< Variable2 "ConvRelEv" >}} are set to <tt>0</tt>, indicating that they will not be used as criteria for convergence.  

These other available criteria use errors defined as follows:

{{< Variable2 "ConvAbsEv" >}}:  The change in each eigenvalue is found and the sum of these changes must be smaller than the threshold.
$ \epsilon = \sum_{j=1}^{N_{occ}} \vert \epsilon_j^{out}-\epsilon_j^{inp}\vert. $

{{< Variable2 "ConvRelEv" >}}: The eigenvalues are scaled by the number of electrons, otherwise as above.
$\epsilon = {1 \over N} \sum_{j=1}^{N_{occ}} \vert \epsilon_j^{out}-\epsilon_j^{inp}\vert.$

To use them, set the relevant variable equal to a number indicating your desired error.  Set the other three convergence variables to zero.

###  Eigensolver  

In each iteration of the self-consistency problem described above, {{< octopus >}} must diagonalize the Kohn-Sham input Hamiltonian to obtain the output eigenfunctions and eigenvalues.  This diagonalization is done by the eigensolver, also an iterative procedure.  There are several options for the iterative scheme used to diagonalize the Hamiltonian, which are specified in the documentation for the variable 
{{< Variable2 "Eigensolver" >}}.  You may specify the threshhold for considering this iterative diagonalization finished with the variable {{< Variable2 "EigensolverTolerance" >}}. The variable {{< Variable2 "EigensolverMaxIter" >}} sets a maximum number of steps in the diagonalization, so that if it is reached, the diagonalization is considered finished even if some of the eigenvectors are not fully converged.  

During each self-consistent field cycle iteration {{< octopus >}} reports the eigenvalues it has obtained by diagonalizing the Hamiltonian, and how many of those eigenvectors are fully converged.

```text
                                                                              
*********************** SCF CYCLE ITER -    3 ************************
 etot =  4.52092631E+00 abs_ev   =  7.91E+02 rel_ev   =  4.28E+01
                        abs_dens =  5.16E-01 rel_dens =  1.43E-02
Matrix vector products:   3669
Converged eigenvectors:      6
Eigenvalues [H]
 -st  Spin   Eigenvalue     Occupation       Error
   1   --    -1.288198       2.000000      (7.2E-07)
   2   --    -0.830676       2.000000      (1.0E-06)
   3   --    -0.826885       2.000000      (8.8E-07)
   4   --    -0.808297       2.000000      (6.2E-07)
...
```
</pre>

It is not too important whether the eigenvectors are not converged in the SCF steps, only whether they are converged at the end.

###  LCAO  

Since the solution of the ground-state problem is done iteratively, we need an initial guess; a set of initial Kohn-Sham orbitals. If we are doing the calculation with pseudopotentials (as opposed to model potentials defined by the user), we can use the pseudo-orbitals that are used to generate the pseudopotential. By default, the program will fill the initial guess states with pseudo-orbitals. Whether or not this is done is determined by the variable {{< Variable2 "LCAOStart" >}}. The guess density is the sum of the atomic densities.

Note, however, that those pseudo-orbitals are not passed directly to the iterative cycle. Instead, the code performs an initial diagonalization with the Hamiltonian generated by the guess density. Therefore, the SCF cycle is started with the linear combination of those atomic orbitals, with the coefficients that result of that diagonalization (LCAO stands for linear combination of atomic orbitals). In other words, the first step of the SCF cycle is performed inside the LCAO subspace, whereas the following steps are performed in the full space.

This diagonalization will typically be done with a number of pseudo-orbitals that is larger than the number that will be used later in the KS SCF cycle. Once we diagonalize that LCAO matrix, we take the lowest lying eigenstates to proceed with the calculation. There is some default number of pseudo-orbitals that will be used, but one can change it making use of variable {{< Variable2 "LCAODimension" >}}.

If you set {{< Variable2 "SCFinLCAO" >}}, the LCAO calculation will be performed self-consistently. Or, in other words, the whole SCF cycle will be done inside the LCAO subspace.


##  Unoccupied states  

This is {{< Variable2 "CalculationMode" >}} {{< code "<nowiki>=</nowiki>unocc" >}}. The purpose of this run mode is to calculate higher lying Kohn-Sham orbitals. For that purpose, it reads the restart information from a converged previous ground-state calculation, and builds the corresponding Hamiltonian. Then, it calculates the unoccupied eigenvalues and eigenfunctions. The number of unoccupied orbitals calculated is given by the {{< Variable2 "ExtraStates" >}}.

{{< manual_foot prev="Manual:Troubleshooting" next="Manual:Time-Dependent" >}}
---------------------------------------------
