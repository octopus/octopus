---
title: "Use of symmetries in optical spectra from time-propagation"
tags: ["Advanced", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-propagation_spectrum"]
series: "Tutorial"
---


In this tutorial we will see how the spatial symmetries of a molecule can be used to reduce the number of time-propagations required to compute the absorption cross-section.

##  Introduction  

The dynamic polarizability (related to optical absorption cross-section via $\sigma = \frac{4 \pi \omega}{c} \mathrm{Im} \alpha $) is, in its most general form, a 3x3 tensor. The reason is that we can shine light on the system polarized in any of the three Cartesian axes, and for each of these three cases measure how the dipole of the molecule oscillates along the three Cartesian axes. This usually means that to obtain the full dynamic polarizability of the molecule we usually need to apply 3 different perturbations along $x, y, z\,$, by setting {{< Variable2 "TDPolarizationDirection" >}} to 1, 2, or 3.

However, if the molecule has some symmetries, it is in general possible to reduce the total number of calculations required to obtain the tensor from 3 to 2, or even 1.[^footnote-1]

To use this formalism in {{< octopus >}}, you need to supply some extra information. The most important thing that the code requires is the information about equivalent axes, that is, directions that are related through some symmetry transformation. Using these axes, we construct a reference frame and specify it with the {{< Variable2 "TDPolarization" >}} block. Note that these axes need not be orthogonal, but they must be linearly-independent. The {{< Variable2 "TDPolarizationEquivAxes" >}} tells the code how many equivalent axes there are. Ideally, the reference frame should be chosen to maximize the number of equivalent axes. When using three equivalent axes, an extra input variable, {{< Variable2 "TDPolarizationWprime" >}}, is also required.

Let us give a couple of examples, which should make all these concepts easier to understand.

##  Methane  

As seen in previous tutorials, the methane molecule has $T_d$ symmetry. This means that it is trivial to find three linearly-independent equivalent axes such that only one time-propagation is needed to obtain the whole tensor. As it happens, we can use the usual $x$, $y$, and $z$ directions, with all of them being equivalent (check the [https://en.wikipedia.org/wiki/Symmetry_operation symmetry operations] of the $T_d$  [https://en.wikipedia.org/wiki/Molecular_symmetry point group] if you are not convinced). Therefore we can perform just one propagation with the perturbation along the $x$ direction adding the following to the input file used previously: 

```text
 %{{< Variable2 "TDPolarization" >}}
  1 | 0 | 0
  0 | 1 | 0
  0 | 0 | 1
 %
 {{< Variable2 "TDPolarizationDirection" >}} = 1
 {{< Variable2 "TDPolarizationEquivAxes" >}} = 3
 %{{< Variable2 "TDPolarizationWprime" >}}
  0 | 0 | 1
 %
```

Note that we had omitted the blocks {{< Variable2 "TDPolarization" >}} and {{< Variable2 "TDPolarizationWprime" >}} in the previous tutorials, as these are their default 

Once the time-propagation is finished, you will find, as usual, a {{< file "td.general/multipoles" >}} file. This time the file contains all the necessary information for {{< octopus >}} to compute the full tensor, so running the  {{< file "oct-propagation_spectrum" >}} utility will produce two files: {{< file "cross_section_vector.1" >}} and {{< file "cross_section_tensor" >}}. The later should look like this (assuming you used the converged grid parameters from the
[Convergence of the optical spectra](../Convergence of the optical spectra) tutorial):  
```text

- nspin         1
- kick mode    0
- kick strength    0.005291772086
- pol(1)           1.000000000000    0.000000000000    0.000000000000
- pol(2)           0.000000000000    1.000000000000    0.000000000000
- pol(3)           0.000000000000    0.000000000000    1.000000000000
- direction    1
- Equiv. axes  3
- wprime           0.000000000000    0.000000000000    1.000000000000
- kick time        0.000000000000
-       Energy         (1/3)*Tr[sigma]    Anisotropy[sigma]      sigma(1,1,1)        sigma(1,2,1)        sigma(1,3,1)       ...
-        [eV]               [A^2]               [A^2]               [A^2]               [A^2]               [A^2]           ...
      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00    ...
      0.10000000E-01     -0.21016843E-08      0.17355148E-11     -0.21016843E-08      0.12207713E-11      0.88665555E-13    ...
      0.20000000E-01     -0.83911822E-08      0.69337823E-11     -0.83911822E-08      0.48772812E-11      0.35411614E-12    ...
...
```
</pre>

Try comparing the spectrum for each component of the $\sigma$ tensor.

##  Linear molecule  

Now let us look at a linear molecule. In this case, you might think that we need two calculations to obtain the whole tensor, one for the direction along the axis of the molecule, and another for the axis perpendicular to the molecule. The fact is that we need only one, in a specially chosen direction, so that our field has components both along the axis of the molecule and perpendicular to it. Let us assume that the axis of the molecule is oriented along the $x\,$ axis. Then we can use

```text
 %{{< Variable2 "TDPolarization" >}}
  1/sqrt(2) | -1/sqrt(2) | 0
  1/sqrt(2) |  1/sqrt(2) | 0
  1/sqrt(2) |  0         | 1/sqrt(2)
 %
 {{< Variable2 "TDPolarizationDirection" >}} = 1
 {{< Variable2 "TDPolarizationEquivAxes" >}} = 3
 %{{< Variable2 "TDPolarizationWprime" >}}
  1/sqrt(2) |  0         | 1/sqrt(2)
 %
```

You should try to convince yourself that the three axes are indeed equivalent and linearly independent. The first and second axes are connected through a simple reflection in the $xz$ plane, transforming the $y$ coordinate from $-1/\sqrt{2}$ into $1/\sqrt{2}$. {{< Variable2 "TDPolarizationWprime" >}} should be set to the result obtained by applying the inverse symmetry operation to the third axis. This actually leaves the third axis unchanged.

##  Planar molecule  

Finally, let us look at a general planar molecule (in the $xy$ plane). In principle we need only two calculations (that is reduced to one if more symmetries are present like, ''e.g.'', in benzene). In this case we chose one of the polarization axes on the plane and the other two rotated 45 degrees:

```text
 %{{< Variable2 "TDPolarization" >}}
  1/sqrt(2) | 0 | 1/sqrt(2)
  1/sqrt(2) | 0 |-1/sqrt(2)
  0         | 1 | 0
 %
 
 {{< Variable2 "TDPolarizationEquivAxes" >}} = 2
```

In this case, we need two runs, one for {{< Variable2 "TDPolarizationDirection" >}} equal to 1, and another for it equal to 3. Note that if there are less than 3 equivalent axes, {{< Variable2 "TDPolarizationWprime" >}} is irrelevant.

##  References  
<references/>

{{Tutorial_foot|series=Optical response|prev=Triplet Excitations|next=}}








---------------------------------------------
[^footnote-1]: {{< Article title="On the use of Neumann's principle for the calculation of the polarizability tensor of nanostructures" authors="M.J.T. Oliveira, A. Castro, M.A.L. Marques, and A. Rubio" journal="J. Nanoscience and Nanotechnology" volume="8" pages="1-7" year="2008" doi="10.1166/jnn.2008.142" arxiv="0710.2624v1" >}}

