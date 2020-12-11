---
title: "Periodic systems"
tags: ["Basic", "Ground State", "Unoccupied", "Bulk", "Pseudopotentials", "DFT", "Band Structure", "DOS"]
series: "Tutorial"
---


The extension of a ground-state calculation to a periodic system is quite straightforward in {{< octopus >}}. In this tutorial we will explain how to perform some basic calculation using bulk silicon as an example.

##  Input  
As always, we will start with a simple input file. In this case we will use a primitive cell of Si, composed of two atoms. 

```text
 {{< Variable2 "CalculationMode" >}} = gs
 
 {{< Variable2 "PeriodicDimensions" >}} = 3
 
 {{< Variable2 "Spacing" >}} = 0.5
 
 %{{< Variable2 "LatticeVectors" >}}
   0.0 | 0.5 | 0.5 
   0.5 | 0.0 | 0.5
   0.5 | 0.5 | 0.0
 %
 
 a = 10.18
 %{{< Variable2 "LatticeParameters" >}}
  a | a | a
 %
 
 %{{< Variable2 "ReducedCoordinates" >}}
  "Si" | 0.0 | 0.0 | 0.0 
  "Si" | 1/4 | 1/4 | 1/4 
 %
 
 nk = 2
 %{{< Variable2 "KPointsGrid" >}}
   nk |  nk |  nk
  0.5 | 0.5 | 0.5
  0.5 | 0.0 | 0.0
  0.0 | 0.5 | 0.0
  0.0 | 0.0 | 0.5
 %
 {{< Variable2 "KPointsUseSymmetries" >}} = yes
 
 {{< Variable2 "ExtraStates" >}} = 1
 {{< Variable2 "Output" >}} = dos
```

Lets see more in detail some of the input variables:

* <tt>{{< Variable2 "PeriodicDimensions" >}} = 3</tt>: this input variable must be set equal to the number of dimensions you want to consider as periodic. Since the system is 3D (<tt>{{< Variable2 "Dimensions" >}} = 3</tt> is the default), by setting this variable to 3 we impose periodic boundary conditions at all borders. This means that we have a fully periodic infinite crystal.

* <tt>{{< Variable2 "LatticeVectors" >}}</tt> and <tt>{{< Variable2 "LatticeParameters" >}}</tt>: these two blocks are used to define the primitive lattice vectors that determine the unit cell. <tt>{{< Variable2 "LatticeVectors" >}}</tt> defines the direction of the vectors, while {{< Variable2 "LatticeParameters" >}} defines their length.

* <tt>{{< Variable2 "ReducedCoordinates" >}}</tt>: the position of the atoms inside the unit cell, in reduced coordinates.

* <tt>{{< Variable2 "KPointsGrid" >}}</tt>: this specifies the ''k''-point grid to be used in the calculation. Here we employ a 2x2x2 Monkhorst-Pack grid with four shifts. The first line of the block defines the number of ''k''-points along each axis in the Brillouin zone. Since we want the same number of points along each direction, we have defined the auxiliary variable <tt>nk = 2</tt>. This will be useful later on to study the convergence with respect to the number of ''k''-points. The other four lines define the shifts, one per line, expressed in reduced coordinates of the Brillouin zone. Alternatively, one can also define the reciprocal-space mesh by explicitly setting the position and weight of each ''k''-point using the {{< Variable2 "KPoints" >}} or {{< Variable2 "KPointsReduced" >}} variables.

* <tt>{{< Variable2 "KPointsUseSymmetries" >}} = yes</tt>: this variable controls if symmetries are used or not. When symmetries are used, the code shrinks the Brillouin zone to its irreducible portion and the effective number of ''k''-points is adjusted. 

* <tt>{{< Variable2 "Output" >}} = dos</tt>: we ask the code to output the density of states.

Here we have taken the value of the grid spacing to be 0.5 bohr. Although we will use this value throughout this tutorial, remember that in a real-life calculation the convergence with respect to the grid spacing must be performed for all quantities of interest.

Note that for periodic systems the default value for the <tt>{{< Variable2 "BoxShape" >}}</tt> variable is <tt>parallelepiped</tt>, although in this case the name can be misleading, as the actual shape also depends on the lattice vectors. This is the only box shape currently available for periodic systems.

##  Output  
Now run {{< octopus >}} using the above input file. Here are some important things to note from the output.

```text

******************************** Grid ********************************
Simulation Box:
  Type = parallelepiped
  Lengths [b] = (   3.627,   3.627,   3.627)
  Octopus will run in 3 dimension(s).
  Octopus will treat the system as periodic in 3 dimension(s).

  Lattice Vectors [b]
    0.000000    5.130000    5.130000
    5.130000    0.000000    5.130000
    5.130000    5.130000    0.000000
  Cell volume =           270.0114 [b^3]
  Reciprocal-Lattice Vectors [b^-1]
   -0.612396    0.612396    0.612396
    0.612396   -0.612396    0.612396
    0.612396    0.612396   -0.612396
Main mesh:
  Spacing [b] = ( 0.484, 0.484, 0.484)    volume/point [b^3] =      0.08000
  - inner mesh =       3375
  - total mesh =      10639
  Grid Cutoff [H] =    21.095389    Grid Cutoff [Ry] =    42.190778
**********************************************************************
```
</pre>
Here {{< octopus >}} outputs some information about the cell in real and reciprocal space.

```text

***************************** Symmetries *****************************
Space group No.  227
International: Fd-3m
Schoenflies: Oh^7
  Index                Rotation matrix                      Fractional translations
    1 :     1   0   0     0   1   0     0   0   1      0.000000    0.000000    0.000000
    2 :     0  -1   1     0  -1   0     1  -1   0      0.000000    0.000000    0.000000
    3 :    -1   0   0    -1   0   1    -1   1   0      0.000000    0.000000    0.000000
    4 :     0   1  -1     1   0  -1     0   0  -1      0.000000    0.000000    0.000000
    5 :    -1   1   0    -1   0   0    -1   0   1      0.000000    0.000000    0.000000
    6 :     0   0   1     1   0   0     0   1   0      0.000000    0.000000    0.000000
    7 :     1  -1   0     0  -1   1     0  -1   0      0.000000    0.000000    0.000000
    8 :     0   0  -1     0   1  -1     1   0  -1      0.000000    0.000000    0.000000
    9 :     0  -1   0     1  -1   0     0  -1   1      0.000000    0.000000    0.000000
   10 :    -1   0   1    -1   1   0    -1   0   0      0.000000    0.000000    0.000000
   11 :     0   1   0     0   0   1     1   0   0      0.000000    0.000000    0.000000
   12 :     1   0  -1     0   0  -1     0   1  -1      0.000000    0.000000    0.000000
   13 :     0   0  -1     1   0  -1     0   1  -1      0.000000    0.000000    0.000000
   14 :    -1   1   0    -1   0   1    -1   0   0      0.000000    0.000000    0.000000
   15 :     1  -1   0     0  -1   0     0  -1   1      0.000000    0.000000    0.000000
   16 :     0   0   1     0   1   0     1   0   0      0.000000    0.000000    0.000000
   17 :     1   0  -1     0   1  -1     0   0  -1      0.000000    0.000000    0.000000
   18 :     0  -1   0     0  -1   1     1  -1   0      0.000000    0.000000    0.000000
   19 :     0   1   0     1   0   0     0   0   1      0.000000    0.000000    0.000000
   20 :    -1   0   1    -1   0   0    -1   1   0      0.000000    0.000000    0.000000
   21 :     0   1  -1     0   0  -1     1   0  -1      0.000000    0.000000    0.000000
   22 :     1   0   0     0   0   1     0   1   0      0.000000    0.000000    0.000000
   23 :    -1   0   0    -1   1   0    -1   0   1      0.000000    0.000000    0.000000
   24 :     0  -1   1     1  -1   0     0  -1   0      0.000000    0.000000    0.000000
Info: The system has    24 symmetries that can be used.
**********************************************************************
```
</pre>
This block tells us about the space-group and the symmetries found for the specified structure.

```text

     2 k-points generated from parameters :
 ---------------------------------------------------
    n =    2    2    2
 
    s1  =  0.50  0.50  0.50
 
    s2  =  0.50  0.00  0.00
 
    s3  =  0.00  0.50  0.00
 
    s4  =  0.00  0.00  0.50
  
 index |    weight    |             coordinates              |
     1 |     0.250000 |    0.250000    0.000000    0.000000  |
     2 |     0.750000 |   -0.250000    0.250000    0.250000  |
```
</pre>
Next we get the list of the ''k''-points in reduced coordinates and their weights. Since symmetries are used, only two ''k''-points are generated. If we had not used symmetries, we would have 32 ''k''-points instead.

The rest of the output is much like its non-periodic counterpart. After a few iterations the code should converge:
```text

*********************** SCF CYCLE ITER -    8 ************************
 etot  = -7.92843852E+00 abs_ev   =  4.39E-06 rel_ev   =  1.79E-05
 ediff =       -1.19E-06 abs_dens =  3.49E-05 rel_dens =  4.37E-06
Matrix vector products:     77
Converged eigenvectors:     10

-  State  KPoint  Eigenvalue [H]  Occupation    Error
      1       1       -0.257101    2.000000   (6.5E-07)
      1       2       -0.184777    2.000000   (9.3E-07)
      2       2       -0.079243    2.000000   (9.4E-07)
      2       1        0.010873    2.000000   (8.4E-07)
      3       2        0.023247    2.000000   (8.2E-07)
      4       2        0.073689    2.000000   (8.6E-07)
      3       1        0.128013    2.000000   (9.0E-07)
      4       1        0.128013    2.000000   (9.6E-07)
      5       2        0.209032    0.000000   (9.3E-07)
      5       1        0.232110    0.000000   (9.2E-07)

Density of states:

----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
------------------------------------------------------%---------------
%---------%--------------%------------%-%------%------%-----------%--%
                                                      ^


Elapsed time for SCF step     8:          0.10
**********************************************************************


             Info: Writing states. 2018/08/06 at 17:15:24


        Info: Finished writing states. 2018/08/06 at 17:15:24

Info: SCF converged in    8 iterations
```
</pre>

As usual, the {{< file "static/info" >}} file contains the most relevant information concerning the calculation. Since we asked the code to output the density of states, we also have a few new files in the {{< file "static" >}} directory:

* {{< file "dos-XXXX.dat" >}}: the band-resolved density of states (DOS);

* {{< file "total-dos.dat" >}}: the total DOS (summed over all bands);

* {{< file "total-dos-efermi.dat" >}}: the Fermi Energy in a format compatible with {{< file "total-dos.dat" >}}.

Of course you can tune the output type and format in the same way you do in a finite-system calculation.

##  Convergence in k-points  

Similar to the convergence in spacing, a convergence must be performed for the sampling of the Brillouin zone. To do this one must try different numbers of ''k''-points along each axis in the Brillouin zone. This can easily be done by changing the value of the <tt>nk</tt> auxiliary variable in the previous input file. You can obviously do this by hand, but this is something that can also be done with a script. Here is such a script, which we will name {{< file "kpts.sh" >}}:
```text

-!/bin/bash
echo "-nk Total Energy" > kpts.log
list="2 4 6 8"
for nkpt in $list
do
    sed -i "s/nk = [0-9]\+/nk = $nkpt/g" inp
    octopus >& out-$nkpt
    energy=`grep Total static/info  | head -2 | tail -1 | cut -d "=" -f 2`
    echo $nkpt $energy >> kpts.log
    rm -rf restart
done
```
</pre>
After running the script (<tt>source kpts.sh</tt>), you should obtain a file called {{< file "kpts.log" >}} that should look like this:
```text

-nk Total Energy
2 -7.92843852
4 -7.93452207
6 -7.93462084
8 -7.93462226
```
</pre>
As you can see, the total energy is converged to within 0.0001 hartree for <tt>nk = 6</tt>.
For an extended range, e.g. from 2 to 12, the script should be
```text

-!/bin/bash
echo "-nk Total Energy" > kpts.log
list="2 4 6 8 10 12"
for nkpt in $list
do
    sed -i  "s/nk = [0-9][0-9]*/nk = $nkpt/g" inp
    octopus >& out-$nkpt
    energy=`grep Total static/info  | head -2 | tail -1 | cut -d "=" -f 2`
    echo $nkpt $energy >> kpts.log
    rm -rf restart
done
```
</pre>
Then, the {{< file "kpts.log" >}} file should look like this:
```text

-nk Total Energy
2 -7.92843852
4 -7.93452207
6 -7.93462084
8 -7.93462226
10 -7.93462251
12 -7.93462217
```
</pre>

##  Band-structure  

We now proceed with the calculation of the band-structure of Si. In order to compute a band-structure, we must perform a non-self-consistent calculation, using the density of a previous ground-state calculation. So the first step is to obtain the initial ground-state. To do this, rerun the previous input file, but changing the number of k-points to <tt>nk=6</tt>. Next, modify the input file such that it looks like this:

```text
 {{< Variable2 "CalculationMode" >}} = unocc
 
 {{< Variable2 "PeriodicDimensions" >}} = 3
 
 {{< Variable2 "Spacing" >}} = 0.5
 
 %{{< Variable2 "LatticeVectors" >}}
   0.0 | 0.5 | 0.5 
   0.5 | 0.0 | 0.5
   0.5 | 0.5 | 0.0
 %
 
 a = 10.18
 %{{< Variable2 "LatticeParameters" >}}
  a | a | a
 %
 
 %{{< Variable2 "ReducedCoordinates" >}}
  "Si" | 0.0 | 0.0 | 0.0 
  "Si" | 1/4 | 1/4 | 1/4 
 %
 
 {{< Variable2 "ExtraStates" >}} = 10
 {{< Variable2 "ExtraStatesToConverge" >}} = 5
 
 %{{< Variable2 "KPointsPath" >}}
   10 |  10 |  15
  0.5 | 0.0 | 0.0  - L point
  0.0 | 0.0 | 0.0  - Gamma point
  0.0 | 0.5 | 0.5  - X point
  1.0 | 1.0 | 1.0  - Another Gamma point
 %
 {{< Variable2 "KPointsUseSymmetries" >}} = no
```

Here are the things we changed:

* <tt>{{< Variable2 "CalculationMode" >}} = unocc</tt>: we are now performing a non-self-consistent calculation, so we use the <tt>unoccupied</tt> calculation mode;

* <tt>{{< Variable2 "ExtraStates" >}} = 10</tt>: this is the number of unoccupied bands to calculate;

* <tt>{{< Variable2 "ExtraStatesToConverge" >}} = 5</tt>: the highest unoccupied states are very hard to converge, so we use this variable to specify how many unoccupied states are considered for the stopping criterion of the non-self-consistent run.

* <tt>{{< Variable2 "KPointsPath" >}}</tt>: this block is used to specify that we want to calculate the band structure along a certain path in the Brillouin zone. This replaces the {{< Variable2 "KPointsGrid" >}} block. The first row describes how many ''k''-points will be used to sample each segment. The next rows are the coordinates of the ''k''-points from which each segment starts and stops. In this particular example, we chose the following path: L-Gamma, Gamma-X, X-Gamma using a sampling of 10-10-15 ''k''-points.

* <tt>{{< Variable2 "KPointsUseSymmetries" >}} = no</tt>: we have turned off the use of symmetries.

After running {{< octopus >}} with this input file, you should obtain a file named {{< file "bandstructure" >}} inside the {{< file "static" >}} directory. This is how the first few lines of the file should look like:

```text

- coord. kx ky kz (red. coord.), bands:    14 [H]
     0.00000000    0.50000000    0.00000000    0.00000000    -0.19984750    -0.10463050     0.11059562     0.11059567     0.21320524     0.27708534     0.27708537     0.43517458     0.54466879     0.54466894     0.57089044     0.57496796     0.57496801     0.63339690
     0.02640129    0.45000000   -0.00000000   -0.00000000    -0.20506492    -0.09711271     0.11125581     0.11125586     0.21395854     0.27788925     0.27788927     0.43640709     0.51950978     0.51950989     0.55510464     0.60282729     0.60282738     0.64782337
     0.05280258    0.40000000   -0.00000000   -0.00000000    -0.21725636    -0.07804331     0.11323409     0.11323414     0.21620689     0.28007561     0.28007563     0.43967737     0.48696242     0.48696251     0.52612903     0.64322288     0.64322300     0.66832900
...
```
</pre>

{{< figure src="/images/Si_bandstructure.png" width="500px" caption="Band structure of bulk silicon. The zero of energy has been shifted to the maximum of the occupied bands." >}}

The first column is the coordinate of the ''k''-point along the path. The second, third, and fourth columns are the reduced coordinates of the ''k''-point. The following columns are the eigenvalues for the different bands. In this case there are 14 bands (4 occupied and 10 unoccupied).

On the right you can see the plot of the band structure. This plot shows the occupied bands (purple) and the first 5 unoccupied bands (green). If you are using gnuplot, you can obtain a similar plot with the following command:
```text

plot for [col=5:5+9] 'static/bandstructure' u 1:(column(col)) w l notitle ls 1
```
</pre>

Note that when using the {{< Variable2 "KPointsPath" >}} input variable, {{< octopus >}} will run in a special mode, and the restart information of the previous ground-state calculation will not be altered in any way. The code informs us about this just before starting the unoccupied states iterations:
```text

Info: The code will run in band structure mode.
     No restart information will be printed.
```
</pre>

{{Tutorial_foot|series=Octopus basics|prev=Centering a geometry|next=Time-dependent propagation}}









---------------------------------------------
