---
title: "Geometry optimization"
tags: ["Tutorial", "Beginner", "Geometry Optimization", "Molecule", "Pseudopotentials", "DFT", "Total Energy", "Forces"]
series: "Tutorial"
---


In this tutorial we will show how to perform geometry optimizations using {{< octopus >}}.

## Methane molecule 

###  Manual optimization  
We will start with the methane molecule. Before using any sophisticated algorithm for the geometry optimization, it is quite instructive to do it first "by hand". Of course, this is only practical for systems with very few degrees of freedom, but this is the case for methane, as we only need to optimize the CH bond-length. We will use the following files.

####  inp  

```text
 {{< Variable2 "CalculationMode" >}} = gs
 {{< Variable2 "UnitsOutput" >}} = eV_Angstrom
 
 {{< Variable2 "Radius" >}} = 3.5*angstrom
 {{< Variable2 "Spacing" >}} = 0.18*angstrom
 
 CH = 1.2*angstrom
 %{{< Variable2 "Coordinates" >}}
   "C" |           0 |          0 |           0 
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3) 
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3) 
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
```

This input file is similar to ones used in other tutorials. We have seen in the [Total energy convergence](../Total_energy_convergence) tutorial that for these values of the spacing and radius the total energy was converged to within 0.1 eV.

####  ch.sh  
You can run the code by hand, changing the bond-length for each calculation, but you can also use the following script:
```text

-!/bin/bash
echo "-CH distance   Total Energy" > ch.log
list="0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30"
for CH in $list
do
    sed -ie "/CH = /s|[0-9.]\+|$CH|" inp
    octopus >& out-$CH
    energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
    echo $CH $energy >> ch.log
    rm -rf restart
done
```
</pre>

{{< figure src="/Methane Energy vs CH.png" width="500px" caption="420px" >}}

####  Results  
Once you run the calculations, you should obtain values very similar to the following ones:

```text

-CH distance   Total Energy
0.90 -215.22444186
0.95 -217.00785398
1.00 -218.09421507
1.05 -218.64692827
1.10 -218.78607337
1.15 -218.61468911
1.20 -218.20448384
1.25 -217.61525239
1.30 -216.89630205
```
</pre>

You can see on the right how this curve looks like when plotted. By simple inspection we realize that the equilibrium CH distance that we obtain is around 1.1 Å. This should be compared to the experimental value of 1.094 Å. From this data can you calculate the frequency of the vibrational breathing mode of methane? 

###  FIRE  
Since changing the bond-lengths by hand is not exactly the simplest method to perform the geometry optimization, we will now use automated algorithms for doing precisely that. Several of these algorithms are available in {{< octopus >}}. The default and recommend one is the the FIRE method[^footnote-1]
. 
{{< octopus >}} can also use the [http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html GSL optimization routines], so all the methods provided by this library are available.

####  Input  

First we will run a simple geometry optimization using the following input file:

```text
 {{< Variable2 "FromScratch" >}} = yes
 {{< Variable2 "CalculationMode" >}} = go
 {{< Variable2 "UnitsOutput" >}} = eV_Angstrom
 
 {{< Variable2 "BoxShape" >}} = sphere
 {{< Variable2 "Radius" >}} = 5.0*angstrom
 {{< Variable2 "Spacing" >}} = 0.18*angstrom
 
 CH = 1.2*angstrom
 %{{< Variable2 "Coordinates" >}}
   "C" |           0 |          0 |           0 
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3) 
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3) 
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
```

In this case we are using default values for all the input variables related to the geometry optimization. Compared to the input file above, we have changed the {{< Variable2 "CalculationMode" >}} to <tt>go</tt>. Because the default minimum box of spheres centered on the atoms has some caveats when performing a geometry optimization (see bellow for more details), we also changed the box shape to a single sphere with a 5 Å radius. Like always, convergence of the final results should be checked with respect to the mesh parameters. In this case you can trust us that a radius of 5 Å is sufficient.

The geometry in the input file will be used as a starting point. Note that although a single variable CH was used to specify the structure, the coordinates will be optimized independently.

####  Output  
When running {{< octopus >}} in geometry optimization mode the output will be a bit different from the normal ground-state mode. Since each minimization step requires several SCF cycles, only one line of information for each SCF iteration will printed: 

```text

...
Info: Starting SCF iteration.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    1 : etot -2.28028257E+02 : abs_dens 2.03E+00 : force  1.51E+00 : etime     1.1s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    2 : etot -2.28357657E+02 : abs_dens 2.31E-01 : force  4.90E-01 : etime     1.0s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    3 : etot -2.20168027E+02 : abs_dens 2.44E-01 : force  1.94E-02 : etime     1.0s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    4 : etot -2.18151090E+02 : abs_dens 6.26E-02 : force  2.02E-01 : etime     0.9s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    5 : etot -2.18383137E+02 : abs_dens 5.06E-03 : force  2.01E-02 : etime     0.9s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    6 : etot -2.18212576E+02 : abs_dens 3.38E-03 : force  1.24E-02 : etime     0.8s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    7 : etot -2.18212135E+02 : abs_dens 3.30E-05 : force  1.81E-07 : etime     0.5s

             Info: Writing states. 2018/08/06 at 11:01:41


        Info: Finished writing states. 2018/08/06 at 11:01:41

Info: SCF converged in    7 iterations



++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++ MINIMIZATION ITER -:    1 ++++++++++++++++++++++
  Energy    = -218.2121352982 eV
  Max force =    2.5241058228 eV/A
  Max dr    =    0.0002435392 A
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

...
```
</pre>

At the end of each minimization step the code will print the total energy, the maximum force acting on the ions and the maximum displacement of the ions from the previous geometry to the next one. {{< octopus >}} will also print the current geometry to a file named {{< file "go.XXXX.xyz" >}} (XXXX is the minimization iteration number) in a directory named {{< file "geom" >}}. In the working directory there is also a file named {{< file "work-geom.xyz" >}} that accumulates the geometries of each SCF calculation. The results are also summarized in {{< file "geom/optimization.log" >}}. The most recent individual geometry is in {{< file "last.xyz" >}}.

After some minimization steps the code will exit. 

```text

...

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++ MINIMIZATION ITER -:   32 ++++++++++++++++++++++
  Energy    = -218.7910677857 eV
  Max force =    0.0513974737 eV/A
  Max dr    =    0.0052755032 A
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



Writing final coordinates to min.xyz
Info: Finished writing information to 'restart/gs'.

...
```
</pre>

If the minimization was successful a file named {{< file "min.xyz" >}} will be written in the working directory containing the minimized geometry. In this case it should look like this:
```text

   5
 units: A
      C                   -0.000000   -0.000000   -0.000000
      H                    0.633663    0.633663    0.633663
      H                   -0.633663   -0.633663    0.633663
      H                    0.633663   -0.633663   -0.633663
      H                   -0.633663    0.633663   -0.633663
```
</pre>
Has the symmetry been retained? What is the equilibrium bond length determined? Check the forces in {{< file "static/info" >}}. You can also visualize the optimization process as an animation by loading {{< file "work-geom.xyz" >}} as an XYZ file into XCrySDen.

##  Sodium trimer  

Let's now try something different: a small sodium cluster. This time we will use an alternate method related to conjugate gradients. 

####  Input  
Here is the input file to use

```text
 {{< Variable2 "FromScratch" >}} = yes
 {{< Variable2 "CalculationMode" >}} = go
 {{< Variable2 "UnitsOutput" >}} = eV_Angstrom
 
 {{< Variable2 "Radius" >}} = 9.0*angstrom
 {{< Variable2 "Spacing" >}} = 0.3*angstrom
 
 {{< Variable2 "XYZCoordinates" >}} = "inp.xyz"
 
 {{< Variable2 "GOMethod" >}} = cg_bfgs
```

and the corresponding {{< file "inp.xyz" >}} file
```text

  3

  Na     -1.65104552837    -1.51912734276    -1.77061890357
  Na     -0.17818042320    -0.81880768425    -2.09945089196
  Na      1.36469690888     0.914070162416    0.40877139068
```
</pre>
In this case the starting geometry is a random geometry.

####  Output  
After some time the code should stop successfully. The information about the last minimization iteration should look like this:
```text

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++ MINIMIZATION ITER -:    3 ++++++++++++++++++++++
  Energy    =  -16.9047778355 eV
  Max force =    0.0257973692 eV/A
  Max dr    =    0.1384047905 A
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
```
</pre>
```text
 
```
You might have noticed that in this case the code might perform several SCF cycles for each minimization iteration. This is normal, as the minimization algorithm used here requires the evaluation of the energy and forces at several geometries per minimization iteration. The minimized geometry can be found in the {{< file "min.xyz" >}} file:
```text

   3
 units: A
      Na                  -2.269359   -1.808026   -1.617806
      Na                   0.457604   -0.515049   -2.237456
      Na                   1.347047    0.899166    0.394155
```
</pre>

####  Caveats of using a minimum box  
For this example we are using the minimum box (<tt>{{< Variable2 "BoxShape" >}} = minimum</tt>). This box shape is the one used by default, as it is usually the most efficient one, but it requires some care when using it for a geometry optimization. The reason is that the box is generated at the beginning of the calculation and is not updated when the atoms move. This means that at the end of the optimization, the centers of the spheres composing the box will not coincide anymore with the position of the atoms. This is a problem if the atoms move enough such that their initial and final positions differ significantly. If this is the case, the best is to restart the calculation using the previously minimized geometry as the starting geometry, so that the box is constructed again consistently with the atoms. In this example, this can be done by copying the contents of the {{< file "min.xyz" >}} file {{< file "inp.xyz" >}} and running {{< octopus >}} again without changing input file.

This time the calculation should converge after only one iteration:
```text

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++ MINIMIZATION ITER -:    1 ++++++++++++++++++++++
  Energy    =  -16.8896899780 eV
  Max force =    0.0447462812 eV/A
  Max dr    =    0.0450282495 A
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
```
</pre>
You can now compare the old {{< file "min.xyz" >}} file with the new one:
```text

   3
 units: A
      Na                  -2.252443   -1.788415   -1.588626
      Na                   0.484193   -0.498523   -2.232081
      Na                   1.302019    0.862970    0.362209
```
</pre>
Although the differences are not too large, there are also not negligible.

####  Further improving the geometry  
Imagine that you now want to improve the geometry obtained and reduce the value of the forces acting on the ions. The stopping criterion is controlled by the {{< Variable2 "GOTolerance" >}} variable (and {{< Variable2 "GOMinimumMove" >}}). Its default value is 0.051 eV/Å (0.001 H/b), so let's set it to 0.01 eV/Å:

```text
 {{< Variable2 "FromScratch" >}} = yes
 {{< Variable2 "CalculationMode" >}} = go
 {{< Variable2 "UnitsOutput" >}} = eV_Angstrom
 
 {{< Variable2 "Radius" >}} = 9.0*angstrom
 {{< Variable2 "Spacing" >}} = 0.3*angstrom
 
 {{< Variable2 "XYZCoordinates" >}} = "inp.xyz"
 
 {{< Variable2 "GOMethod" >}} = cg_bfgs
 {{< Variable2 "GOTolerance" >}} = 0.01*eV/angstrom
```

Copy the contents of file {{< file "min.xyz" >}} to file {{< file "inp.xyz" >}} and rerun the code. This time, after some minimization steps {{< octopus >}} should return the following error message:

```text
 **************************** FATAL ERROR *****************************
 *** Fatal Error (description follows)
 *--------------------------------------------------------------------
 * Error occurred during the GSL minimization procedure:
 * iteration is not making progress towards solution
 **********************************************************************
```

This means the minimization algorithm was unable to improve the solutions up to the desired accuracy. Unfortunately there are not many things that can be done in this case. The most simple one is just to restart the calculation starting from the last geometry (use the last {{< file "go.XXXX.xyz" >}} file) and see if the code is able to improve further the geometry.

##  References  
<references/>

<span class=noprint><hr>
Back to [[Tutorials]]









---------------------------------------------
[^footnote-1]: {{< Article title="Structural relaxation made simple" authors="Erik Bitzek, Pekka Koskinen, Franz Gähler, Michael Moseler, and Peter Gumbsch" journal="Phys. Rev. Lett." volume="97" pages="170201" year="2006" url="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201" doi="10.1103/PhysRevLett.97.170201" link="APS" >}}

