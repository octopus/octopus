---
title: "Linear Response"
series: "Manual"
weight: 3
description: "Polarizabilities"
---


{{< Octopus >}} can calculate dynamic polarizabilities and first-order hyperpolarizabilites in a linear-response scheme using the Sternheimer equation. It is also possible to calculate optical spectra with this technique, but it is slower than time-evolution.

#####  Ground state  

The first thing we will need for linear response is a [[Manual:Ground State|Ground State]] calculation. Unlike the Casida approach, when using the Sterheimer equation you needn't do a unoccupied-states calculation. To improve the convergence of the linear-response calculation, it is better to use tightly converged wavefunctions. For example, you can add these parameters to your gs calculation:

```text

EigenSolverFinalTolerance = 1e-10
ConvRelDens = 1e-9
```
</pre>

#####  Input  

The {{< Variable2 "CalculationMode" >}} for polarizability calculations is {{< code "em_resp" >}}. The main parameter you have to specify is the frequency of the perturbation, given by the {{< Variable2 "EMFreqs" >}} block. You can also add an imaginary part to the frequency by setting the variable {{< Variable2 "EMEta" >}}. Adding a small imaginary part is required if you want to get the imaginary part of the polarizability or to calculate polarizabilities near resonance; a reasonable value is {{< code "0.1 eV" >}}.

To get the hyperpolarizabilties, you also have to specify the variable {{< Variable2 "EMHyperpol" >}} with the three coefficients with respect to the base frequency; the three values must sum to zero.

#####  Output  

After running, for each frequency in the input file, {{< octopus >}} will generate a subdirectory under {{< file "em_resp/" >}}. In each subdirectory there is a file called {{< file "alpha" >}} that contains the real part of the polarizability tensor $\alpha_{ij}$ and the average polarizability

$$
\\bar{\\alpha}=\\frac13\\sum\_{i=1}^3\\alpha\_{ii}\\,
$$

The imaginary part $\eta$ is written to file {{< file "eta" >}}. If $\eta > 0$, there is also a file called {{< file "cross_section_tensor" >}} that contains the photo-absorption cross section tensor for that frequency, related to the imaginary part of the polarizability ($\sigma = \frac{4 \pi \omega}{c} \mathrm{Im} \alpha $).

The hyperpolarizability will be in a file called {{< file "beta" >}} at the base frequency, containing all the 27 components and some reduced quantities:

$$
\\beta\_{||\\,i} = \\frac15 \\sum\_{j=1}^3(\\beta\_{ijj}+\\beta\_{jij}+\\beta\_{jji})\\ .
$$

Optionally, Born charges can also be calculated.

####  Finite differences  

In this mode only static polarizability can be obtained. The calculation is done by taking the numerical derivative of the energy with respect to an external static and uniform electric field. To use this, run with {{< Variable2 "ResponseMethod" >}} {{< code "=finite_differences" >}}. {{< octopus >}} will run several ground-state energy calculations and then calculate the polarizability using a finite-differences formula for the derivative. The results will be in the {{< file "em_resp_fd" >}} directory. Hyperpolarizability and Born charges can also be calculated.

{{< manual_foot prev="Manual:Casida" next="Manual:Optimal Control" >}}
---------------------------------------------
