---
title: "Large systems:   the Fullerene molecule"
tags: ["Tutorial", "Advanced", "Ground State", "Molecule", "Pseudopotentials", "DFT"]
series: "Tutorial"
---


In this tutorial, by solving the Fullerene molecule, we show you how to do simulations of systems with a large number of atoms with {{< octopus >}}. Most of the default parameters work fine for large systems, but there are a few things you might need to adjust.

###  The eigensolver  

The default eigensolver of {{< octopus >}} is simple and reliable, but for large systems (> 50 orbitals) it can be quite slow. For these systems it is better to use the RMMDIIS eigensolver. This is a very efficient diagonalizer, but it can have convergence problems. In {{< octopus >}} we follow the algorithm implementation of [http://dx.doi.org/10.1103/PhysRevB.54.11169 Kresse and Furthmüller] that provides a very robust algorithm if the following considerations are taken:

* A good initial approximation of the eigenvalues and eigenvectors is required. In {{< octopus >}}, this means you have to do an LCAO before using the RMMDIIS eigensolver (this is the default). You might as well use an alternative eigensolver for one or two iterations.

* You have to provide a reasonable number of {{< Variable2 "ExtraStates" >}}, it should be around 10% to 20% of the number of occupied states. You must include them by hand.

Let's see an example: this is an input file, {{< file "inp" >}}, to calculate the ground state of the fullerene molecule with the RMMDIIS eigensolver:

```text
 {{< Variable2 "XYZCoordinates" >}} = "C60.xyz"
 {{< Variable2 "UnitsOutput" >}} = eV_Angstrom
 {{< Variable2 "CalculationMode" >}} = gs
 {{< Variable2 "EigenSolver" >}} = rmmdiis
 {{< Variable2 "ExtraStates" >}} = 20
```

In this case, to use RMMDIIS the only thing we need to do is to change the {{< Variable2 "Eigensolver" >}} variable and to include 20 {{< Variable2 "ExtraStates" >}} (by default {{< octopus >}} does not include any unoccupied states).

This is the content of the coordinates file, {{< file "C60.xyz" >}}:

```text

60

 C    -1.415441     3.015027    -1.167954
 C     2.583693     2.292084    -0.721147
 C     0.721155    -2.583700     2.292090
 C     3.015024    -1.167941     1.415429
 C    -3.015019    -1.167950     1.415434
 C    -0.721134    -2.583697     2.292085
 C    -2.583698     2.292079    -0.721148
 C    -3.461397    -0.000006     0.694328
 C     1.415423     3.015023    -1.167951
 C     3.461399     0.000011     0.694319
 C    -1.415430    -3.015023    -1.167948
 C    -2.292086    -0.721156    -2.583690
 C    -0.000003     0.694321    -3.461398
 C     2.292087    -0.721146    -2.583706
 C     1.415436    -3.015019    -1.167948
 C    -2.292087     0.721141    -2.583692
 C     1.167941     1.415428    -3.015025
 C     3.015020    -1.167945    -1.415437
 C     0.694330    -3.461397     0.000006
 C    -2.583701    -2.292092    -0.721145
 C     0.721144     2.583697     2.292082
 C     1.167948     1.415436     3.015018
 C    -0.000002     0.694327     3.461397
 C    -1.167951     1.415435     3.015031
 C    -0.721152     2.583693     2.292084
 C    -0.694329     3.461397    -0.000006
 C     2.583701     2.292093     0.721144
 C     2.292088    -0.721140     2.583692
 C    -1.167941    -1.415428     3.015025
 C    -3.015020     1.167945     1.415437
 C     1.167955    -1.415429     3.015019
 C    -2.292077    -0.721143     2.583701
 C    -2.583697     2.292080     0.721148
 C     0.694318     3.461398    -0.000005
 C     3.015030     1.167958     1.415430
 C    -0.694317    -3.461398     0.000005
 C    -3.015030    -1.167958    -1.415430
 C    -1.167955     1.415429    -3.015019
 C     2.292078     0.721143    -2.583701
 C     2.583697    -2.292079    -0.721148
 C     1.167951    -1.415434    -3.015031
 C     0.721152    -2.583693    -2.292084
 C    -0.721144    -2.583696    -2.292082
 C    -1.167948    -1.415436    -3.015019
 C     0.000002    -0.694327    -3.461397
 C     0.721134     2.583697    -2.292085
 C     3.461397     0.000007    -0.694329
 C     1.415440    -3.015027     1.167954
 C    -2.583694    -2.292084     0.721147
 C    -3.015025     1.167942    -1.415429
 C    -2.292087     0.721146     2.583706
 C    -1.415436     3.015019     1.167948
 C     1.415430     3.015023     1.167947
 C     2.292086     0.721156     2.583690
 C     0.000003    -0.694321     3.461398
 C    -1.415423    -3.015023     1.167952
 C    -3.461399    -0.000011    -0.694320
 C    -0.721155     2.583700    -2.292090
 C     3.015019     1.167950    -1.415434
 C     2.583698    -2.292079     0.721148
```
</pre>

Create both files and run {{< octopus >}}. It should take a few minutes (depending on the speed of your machine). There are some peculiarities of the RMMDIIS eigensolver that you might notice:

* The first self-consistency iterations take less time than the rest of the iterations. To ensure a proper initial approximation, the eigesolver does a different diagonalization during the first steps. 

* With respect to other eigensolvers, it takes many more self-consistency iterations to converge the calculation. In compensation each RMMDIIS step is faster. This is because at each step a small and fixed number of eigensolver steps are done, so in practice the eigensolver and self-consistency convergence are done together.

* At some iterations the eigenvalues might not be ordered -- this is normal. When convergence is reached, the eigenvalues will be properly ordered.

###  Self-consistency convergence  

Once your calculation is finished, you might notice that the final residue in the eigenvectors is rather large (of the order of $10^{-4}$). This is for two reasons: first, it is possible to converge the density without properly converging the eigenvectors sometimes. Try a tighter tolerance and rerun (with restart), adding

```text
 {{< Variable2 "ConvRelDens" >}} = 1e-6
```

to your {{< file "inp" >}} file. Now your eigenvectors should be better converged.

Can you converge all the eigenvectors? How many states are actually properly converged? Suppose you wanted to calculate correct unoccupied states. How would you do it and be sure they were converged?

The second reason for poor convergence is that we are using a large (default) spacing -- check what it is in the output. The problem is just not as well-conditioned with a large spacing. So, you can try using a smaller spacing such as the one you found in the [methane tutorial](../Methane_molecule ).

<span class=noprint><hr>
Back to [[Tutorials]]







---------------------------------------------
