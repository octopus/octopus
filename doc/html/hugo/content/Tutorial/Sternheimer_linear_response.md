---
title: "Sternheimer linear response"
tags: ["Tutorial", "Advanced", "Electromagnetic Response", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "Sternheimer"]
series: "Tutorial"
---


The Sternheimer approach to perturbation theory allows efficient calculations of linear and non-linear response properties.
[^footnote-1]

The basis of this method, just as in standard perturbation theory, is to calculate the variation of the wave-functions $\psi^{1}$ under a given perturbing potential. The advantage of the method is that the variations are obtained by solving the linear equation 

$$
(H^0-\\epsilon^0 + \\omega )|\\psi^{1}\>=-P\_{\\rm c} H^{1}|\\psi^{0}\>\\ ,
$$

that only depends on the occupied states instead of requiring an (infinite) sum over unoccupied states. In the case of (time-dependent) density functional theory the variation of the Hamiltonian includes a term that depends on the variation of the density, so this equation must be solved self-consistently.

To run a Sternheimer calculation with {{< octopus >}}, the only previous calculation you need is a ground-state calculation.
For this tutorial we will use a water molecule, with this basic input file for the ground state:

```text
 {{< Variable2 "CalculationMode" >}} = gs
 
 %{{< Variable2 "Coordinates" >}}
  'O'  |  0.000000  | -0.553586  |  0.000000
  'H'  |  1.429937  |  0.553586  |  0.000000
  'H'  | -1.429937  |  0.553586  |  0.000000
 %
 
 {{< Variable2 "Radius" >}} = 10
 {{< Variable2 "Spacing" >}} = 0.435
 {{< Variable2 "ConvRelDens" >}} = 1e-6
```

We use a tighter setting on SCF convergence (ConvRelDens) which will help the ability of the Sternheimer calculation to converge numerically, and we increase a bit the size of the box as response calculations tend to require more space around the molecule than ground-state calculations to be converged.
[^footnote-2]


After the ground-state calculation is finished, we change the run mode to {{< code "em_resp" >}}, to run a calculation of the electric-dipole response:

```text
 {{< Variable2 "CalculationMode" >}} = em_resp
```

Next, to specify the frequency of the response we use the {{< Variable2 "EMFreqs" >}} block; in this case we will use three values 0.00, 0.15 and 0.30 {{< units "{{< Hartree >" >}}}}:

```text
 %{{< Variable2 "EMFreqs" >}} 
 3 | 0.0 | 0.3
 %
```

and we will also specify a small imaginary part to the frequency of 0.1 {{< units "eV" >}}, which avoids divergence on resonance:

```text
 {{< Variable2 "EMEta" >}} = 0.1*eV
```

and finally we add a specification of the linear solver, which will greatly speed things up compared to the default:

```text
 {{< Variable2 "LinearSolver" >}} = qmr_dotp
 {{< Variable2 "ExperimentalFeatures" >}} = yes
```

In the run, you will see calculations for each frequency for the ''x'', ''y'', and ''z'' directions, showing SCF iterations, each having linear-solver iterations for the individual states' $\psi^{1}$, labelled by the k-point/spin (ik) and state (ist). The norm of $\psi^{1}$, the number of linear-solver iterations (iter), and the residual $\left|(H^0-\epsilon^0 + \omega )|\psi^{1}>+P_{\rm c} H^{1}|\psi^{0}>\right|$ are shown for each. First we see the static response:

```text

****************** Linear-Response Polarizabilities ******************
Wavefunctions type: Complex
Calculating response for   3 frequencies.
**********************************************************************

Info: Calculating response for the x-direction and frequency 0.0000.
Info: EM Resp. restart information will be written to 'restart/em_resp'.
Info: EM Resp. restart information will be read from 'restart/em_resp'.

** Warning:
**   Unable to read response wavefunctions from 'wfs_x_f1+': Initializing to zero.

Info: Finished reading information from 'restart/em_resp'.
--------------------------------------------
LR SCF Iteration:   1
Frequency:             0.000000 Eta :             0.003675
   ik  ist                norm   iters            residual
    1    1            0.216316      21        0.110754E-03
    1    2            1.573042      21        0.286815E-02
    1    3            1.620482      21        0.973599E-02
    1    4            1.177410      21        0.540540E-03

             Info: Writing states. 2016/01/14 at 19:50:27


        Info: Finished writing states. 2016/01/14 at 19:50:27

SCF Residual:     0.122899E+01 (abs),     0.153623E+00 (rel)
```
</pre>

Later will come the dynamical response. The negative state indices listed indicate response for $-\omega$. For each frequency, the code will try to use a saved response density from the closest previously calculated frequency.

```text

Info: Calculating response for the x-direction and frequency 0.1500.
Info: EM Resp. restart information will be written to 'restart/em_resp'.
Info: EM Resp. restart information will be read from 'restart/em_resp'.
Read response density 'rho_0.0000_1'.
Info: Finished reading information from 'restart/em_resp'.
--------------------------------------------
LR SCF Iteration:   1
Frequency:             0.150000 Eta :             0.003675
   ik  ist                norm   iters            residual
    1    1            0.181350       8        0.136183E-02
    1   -1            0.240803      19        0.794480E-04
    1    2            0.953447       8        0.508247E-02
    1   -2            1.708845      19        0.166997E-02
    1    3            0.864834       8        0.858907E-02
    1   -3            1.780595      19        0.984578E-02
    1    4            0.725681       8        0.529092E-02
    1   -4            1.421823      19        0.536713E-03
```
</pre>

At the end, you will have a directory called {{< file "em_resp" >}} containing a subdirectory for each frequency calculated, each in turn containing {{< file "eta" >}} (listing $\eta$ = 0.1 {{< units "eV" >}}), {{< file "alpha" >}} (containing the real part of the polarizability tensor), and {{< file "cross_section" >}} (containing the cross-section for absorption, based on the imaginary part of the polarizability).

For example, {{< file "em_resp/freq_0.0000/alpha" >}} says

```text
 - Polarizability tensor [b^3]
            10.238694           -0.000000           -0.000000
             0.000000           10.771834           -0.000000
            -0.000000           -0.000000            9.677212
 Isotropic average           10.229247
```

Exercise: compare results for polarizability or cross-section to a calculation from time-propagation or the Casida approach.

##  References  
<references/>

<span class=noprint><hr>
Back to [[Tutorials]]









---------------------------------------------
[^footnote-1]: {{< Article title="Time-dependent density functional theory scheme for efficient calculations of dynamic (hyper)polarizabilities" authors="Xavier Andrade, Silvana Botti, Miguel Marques and Angel Rubio" journal="J. Chem. Phys" volume="126" pages="184106" year="2007" doi="10.1063/1.2733666" >}}

[^footnote-2]: {{< Article title="Basis set effects on the hyperpolarizability of CHCl<sub>3</sub>: Gaussian-type orbitals, numerical basis sets and real-space grids" authors="F. D. Vila, D. A. Strubbe, Y. Takimoto, X. Andrade, A. Rubio, S. G. Louie, and J. J. Rehr" journal="J. Chem. Phys." volume="133" pages="034111" year="2010" doi="10.1063/1.3457362" >}}

