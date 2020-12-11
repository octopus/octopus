---
title: "Symmetry"
series: "Manual"
---


There is not much use of symmetry in Octopus.

In finite systems, you will get an analysis just for your information, but which will not be used anywhere in the code. It is of this form (''e.g.'' for silane):

```text
 ***************************** Symmetries *****************************
 Symmetry elements : 4*(C3) 3*(C2) 3*(S4) 6*(sigma)
 Symmetry group    : Td
 **********************************************************************
```

Many symmetries will in fact be broken by the real-space mesh. Since it is always orthogonal, it will break three-fold rotational symmetry of benzene, for example.

In periodic systems, you will also get an analysis, ''e.g'' like this for bulk silicon in its 8-atom convention cell:

```text
 ***************************** Symmetries *****************************
 Space group No.227
  International: Fd -3 m
  International(long): Fd -3 m _1
  Schoenflies: Oh^7
  Multiplicity: 192
 Point group
  International: m -3 m
  Schoenflies: Oh
 Identity has a fractional translation     0.500000    0.500000    0.000000
 Identity has a fractional translation     0.500000    0.000000    0.500000
 Identity has a fractional translation     0.000000    0.500000    0.500000
 Disabling fractional translations. System appears to be a supercell.
 Info: The system has    24 symmetries that can be used.
 **********************************************************************
```

The analysis is done by the library [http://spglib.sourceforge.net/ spglib]. The comments on fractional translations and supercell are due to the use of the conventional cell, since the 2-atom primitive cell does not have orthogonal lattice vectors and therefore cannot be used in Octopus currently.

For large systems, the symmetry analysis might be very time-consuming; in rare cases, the symmetry analysis might crash and stop the calculation. In either situation, you can use the variable {{< Variable2 "SymmetriesCompute" >}} to turn off the symmetry analysis. You can manually identify a direction in which symmetry is broken with {{< Variable2 "SymmetryBreakDir" >}}, which is appropriate for the case that an external field is applied in a particular direction.

Symmetries are used in two ways for periodic systems: first, to reduce the set of k-points needed. This behavior is controlled by {{< Variable2 "KPointsUseSymmetries" >}}. For high-symmetry systems, this can make a dramatic reduction in the time required for a calculation. For example, in our silicon example, if we set

```text
 %KPointsGrid
  4   | 4   | 4
  0.5 | 0.5 | 0.5
 %
```

we will obtain not 4×4×4 = 64 k-points but only 4:

```text
    4 k-points generated from parameters :
  ---------------------------------------------------
     n =    4    4    4      s =  0.50  0.50  0.50
 
  index |    weight    |             coordinates              |
      1 |     0.125000 |    0.125000    0.125000    0.125000  |
      2 |     0.375000 |    0.125000    0.125000    0.375000  |
      3 |     0.375000 |    0.375000    0.375000    0.125000  |
      4 |     0.125000 |    0.375000    0.375000    0.375000  |
```

The density and other scalar quantities are straightforwardly computed from the density due to each k-point and the weight; the same has been done for some vector quantities such as the gradients used in GGA functionals and the forces.. Note, however, that the proper calculation of more complicated tensorial properties such as the response functions in the Sternheimer calculation modes have not been implemented with symmetry. You can also use the symmetry operations to symmetrize the density, via {{< Variable2 "SymmetrizeDensity" >}}. In general, you should beware of use of symmetries for partially periodic systems (''e.g.'' wire or sheet geometry) for which some problems have been found.

The use of symmetry for reducing the number of calculations needed for obtaining absorption spectra by time-propagation is discussed here:  [[Tutorial:Optical_Spectra_from_TD:Symmetries]]

{{< manual_foot prev="Manual:Output" next="Manual:Troubleshooting" >}}
---------------------------------------------
