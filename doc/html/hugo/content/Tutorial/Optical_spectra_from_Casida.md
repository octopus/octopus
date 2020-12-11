---
title: "Optical spectra from Casida"
tags: ["Beginner", "Unoccupied", "Casida", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-casida_spectrum"]
series: "Tutorial"
---


In this tutorial we will again calculate the absorption spectrum of methane, but this time using Casida's equations.


##  Ground-state  

Once again our first step will be the calculation of the ground state. We will use the following input file:

```text
 {{< Variable2 "CalculationMode" >}} = gs
 {{< Variable2 "UnitsOutput" >}} = eV_angstrom
 
 {{< Variable2 "Radius" >}} = 6.5*angstrom
 {{< Variable2 "Spacing" >}} = 0.24*angstrom
 
 CH = 1.097*angstrom
 %{{< Variable2 "Coordinates" >}}
   "C" |           0 |          0 |           0 
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
```

Note that we are using the values for the spacing and radius that were found in the [Convergence of the optical spectra](../Convergence_of_the_optical_spectra) tutorial to converge the absorption spectrum.

##  Unoccupied States  

The Casida equation is a (pseudo-)eigenvalue equation written in the basis of particle-hole states. This means that we need both the occupied states -- computed in the ground-state calculation -- as well as the unoccupied states, that we will now obtain, via a non-self-consistent calculation using the density computed in <tt>gs</tt>. The input file we will use is

```text
 {{< Variable2 "CalculationMode" >}} = unocc
 {{< Variable2 "UnitsOutput" >}} = eV_angstrom
 
 {{< Variable2 "Radius" >}} = 6.5*angstrom
 {{< Variable2 "Spacing" >}} = 0.24*angstrom
 
 CH = 1.097*angstrom
 %{{< Variable2 "Coordinates" >}}
   "C" |           0 |          0 |           0 
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
 
 {{< Variable2 "ExtraStates" >}} = 10
```

Here we have changed the {{< Variable2 "CalculationMode" >}} to unocc and added 10 extra states by setting the {{< Variable2 "ExtraStates" >}} input variable.

By running {{< octopus >}}, you will obtain the first 10 unoccupied states (do not forget to run the ground-state calculation first). The solution of the unoccupied states is controlled by the variables in section {{< Variable2 "SCF::Eigensolver" >}}. You can take a look at the eigenvalues of the unoccupied states in the file {{< file "static/eigenvalues" >}}:

```text

All states converged.
Criterion =      0.100000E-05

Eigenvalues [eV]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --   -16.787587       2.000000      (7.6E-07)
   2   --    -9.464233       2.000000      (7.6E-07)
   3   --    -9.464233       2.000000      (7.6E-07)
   4   --    -9.464233       2.000000      (7.6E-07)
   5   --    -0.287342       0.000000      (9.8E-07)
   6   --     0.822599       0.000000      (8.0E-07)
   7   --     0.822599       0.000000      (8.6E-07)
   8   --     0.822599       0.000000      (6.5E-07)
   9   --     1.599977       0.000000      (7.9E-07)
  10   --     1.599977       0.000000      (9.2E-07)
  11   --     1.599977       0.000000      (7.3E-07)
  12   --     1.813297       0.000000      (7.5E-07)
  13   --     1.813297       0.000000      (7.6E-07)
  14   --     2.309432       0.000000      (8.8E-07)
```
</pre>

##  Casida calculation  

Now modify the {{< Variable2 "CalculationMode" >}} to <tt>casida</tt> and rerun {{< octopus >}}. Note that by default {{< octopus >}} will use all occupied and unoccupied states that it has available. 

Sometimes, it is useful not to use all states. For example, if you have a molecule with 200 atoms and 400 occupied states ranging from -50 to -2 eV, and you are interested in looking at excitations in the visible, you can try to use only the states that are within 10 eV from the Fermi energy. You could select the states to use with {{< Variable2 "CasidaKohnShamStates" >}} or {{< Variable2 "CasidaKSEnergyWindow" >}}.

A new directory will appear named {{< file "casida" >}}, where you can find the file {{< file "casida/casida" >}}:
```text

                E [eV]         <x> [A]         <y> [A]         <z> [A]             <f>
     1  9.27286861E+00  7.88719340E-02  3.02317792E-01 -1.41580187E-01  9.54564687E-02
     2  9.27291282E+00 -3.06114652E-01  5.96207735E-03 -1.57881218E-01  9.62734206E-02
     3  9.27292772E+00 -1.36862383E-01  1.62030534E-01  2.71506435E-01  9.63001401E-02
     4  1.02415035E+01 -3.18483852E-05  1.77939091E-02 -7.67661722E-04  2.84230878E-04
     5  1.02439430E+01 -6.36128746E-05  8.14736582E-03 -5.21336390E-04  5.97390614E-05
     6  1.02564557E+01 -1.79749506E-03 -1.09950087E-01  4.16624712E-03  1.08663407E-02
     7  1.02593552E+01  3.83016766E-03  3.82893874E-03  1.03764102E-01  9.69062218E-03
     8  1.02594070E+01 -1.04409592E-01  2.07398640E-03  3.78813752E-03  9.80169803E-03
...
```
</pre>

The \<x\>, \<y\> and \<z\> are the transition dipole moments:

$$
  \<x\> = \<\\Phi\_0|x|\\Phi\_I\>
  \\,;\\qquad
  \<y\> = \<\\Phi\_0|y|\\Phi\_I\>
  \\,;\\qquad
  \<z\> = \<\\Phi\_0|z|\\Phi\_I\>
$$

where $\Phi_0$ is the ground state and $\Phi_I$ is the given excitation. The
oscillator strength is given by:

$$
  f\_I = \\frac{2 m\_e}{3 \\hbar^2} \\omega\_I \\sum\_{n\\in x,y,z} |\<\\Phi\_0|n|\\Phi\_I\>|^2\\,
$$

as the average over the three directions. The optical absorption spectrum can be given as the "strength function",
which is

$$
  S(\\omega) = \\sum\_I f\_I \\delta(\\omega-\\omega\_I)\\,
$$

Note that the excitations are degenerate with degeneracy 3. This could already be expected from the $T_d$ symmetry of methane.

''Note that within the degenerate subspaces, there is some arbitrariness (possibly dependent on the details of your compilation and machine) in the linear combinations of transitions. Therefore, you should not be concerned if you do not have the same results for the components of the transition dipole moments (above) and analysis of the excitations (below). However, the energies and resulting spectra should agree.''

Further information concerning the excitations can be found in the directory {{< file "casida/casida_excitations" >}}. For example, the first excitation at 9.27 eV is analyzed in the file {{< file "casida/casida_excitations/00001" >}}:
```text

- Energy [eV] =    9.27287E+00
- <x> [A] =    7.88719E-02
- <y> [A] =    3.02318E-01
- <z> [A] =   -1.41580E-01
           1           5           1  -1.5104487570024249E-005
           2           5           1  0.46406175367379687     
           3           5           1  0.86441581738495976     
           4           5           1 -0.18515534208137396     
           1           6           1   3.2556357546167634E-003
...
```
</pre>

These files contain basically the eigenvector of the Casida equation. The first two columns are respectively the index of the occupied and the index of the unoccupied state, the third is the spin index (always 1 when spin-unpolarized), and the fourth is the coefficient of that state in the Casida eigenvector. This eigenvector is normalized to one, so in this case one can say that 74.7% (0.864<sup>2</sup>) of the excitation is from state 3 to state 5 (one of 3 HOMO orbitals->LUMO) with small contribution from some other transitions.

##  Absorption spectrum  

{{< figure src="/Absorption_spectrum_CH4_casida.png" width="500px" caption="Absorption spectrum of CH<sub>4</sub> calculated with time-propagation and with the Casida equation." >}}

To visualize the spectrum, we need to broaden these delta functions with the utility {{< file "oct-casida_spectrum" >}} (run in your working directory, not in {{< file "casida" >}}). It convolves the delta functions with Lorentzian functions. The operation of {{< file "oct-casida_spectrum" >}} is controlled by the variables {{< Variable2 "CasidaSpectrumBroadening" >}} (the width of this Lorentzian), {{< Variable2 "CasidaSpectrumEnergyStep" >}}, {{< Variable2 "CasidaSpectrumMinEnergy" >}}, and {{< Variable2 "CasidaSpectrumMaxEnergy" >}}. If you run  {{< file "oct-casida_spectrum" >}} you obtain the file {{< file "casida/spectrum.casida" >}}. It contains all columns of {{< file "casida" >}} broadened. If you are interested in the total absorption spectrum, then you should plot the first and fifth columns. You should obtain a picture like the one on the right.

Comparing the spectrum obtained with the time-propagation in the [Convergence of the optical spectra](../Convergence of the optical spectra) tutorial with the one obtained with the Casida approach using the same grid parameters, we can see that

* The peaks of the time-propagation are broader. This can be solved by either increasing the total propagation time, or by increasing the broadening in the Casida approach.
* The first two peaks are nearly the same. Probably also the third is OK, but the low resolution of the time-propagation does not allow to distinguish the two close peaks that compose it.
* For high energies the spectra differ a lot. The reason is that we only used 10 empty states in the Casida approach. In order to describe better this region of the spectrum we would need more. This is why one should always check the convergence of relevant peaks with respect to the number of empty states.

You probably noticed that the Casida calculation took much less time than the time-propagation. This is clearly true for small or medium-sized systems. However, the implementation of Casida in {{< octopus >}} has a much worse scaling with the size of the system than the time-propagation, so for larger systems the situation may be different. Note also that in Casida one needs a fair amount of unoccupied states which are fairly difficult to obtain.

If you are interested, you may also compare the Casida results against the Kohn-Sham eigenvalue differences in {{< file "casida/spectrum.eps_diff" >}} and the Petersilka approximation to Casida in {{< file "casida/spectrum.petersilka" >}}.

{{Tutorial_foot|series=Optical response|prev=Convergence of the optical spectra|next=Optical spectra from Sternheimer}}









---------------------------------------------
