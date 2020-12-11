---
title: "Basic input options"
tags: ["Basic", "Ground State", "Molecule", "Pseudopotentials", "DFT", "Total Energy"]
series: "Tutorial"
---


Now we will move to a more complicated (and realistic) input file. We will obtain the ground state of the nitrogen atom. We will introduce several basic input variables and will give a more detailed description of the output for this example.

##  The input files  

This sample input file lets us obtain the ground state of the nitrogen atom, within the LDA approximation, in a closed-shell (unpolarized) configuration (as explained below, you need an auxiliary {{< file ".xyz" >}} input). Note that this is not the correct ground state of the nitrogen atom! However, it will permit us to describe some of the most important input variables:

```text
 {{< Variable2 "CalculationMode" >}} = gs
 {{< Variable2 "UnitsOutput" >}} = eV_Angstrom
 
 Nitrogen_mass = 14.0
 
 %{{< Variable2 "Species" >}}
  'N' | species_pseudo | set | standard | lmax | 1 | lloc | 0 | mass | Nitrogen_mass
 %
 
 {{< Variable2 "XYZCoordinates" >}} = 'N.xyz'
 
 {{< Variable2 "ExtraStates" >}} = 1
 %{{< Variable2 "Occupations" >}}
   2 | 1 | 1 | 1
 %
 
 {{< Variable2 "BoxShape" >}} = sphere
 {{< Variable2 "Radius" >}} = 5.0*angstrom
 {{< Variable2 "Spacing" >}} = 0.18*angstrom
```

We have introduced here several new variables:

* <tt>{{< Variable2 "UnitsOutput" >}} = eV_Angstrom</tt>: Two different unit systems may be used for output: the usual atomic units (which is the default, and the ones used internally in the code); and the system in which the Ångström is substituted for the atomic unit of length, and the electronvolt is substituted for the atomic unit of energy. You can find a more detailed description of units in {{< Octopus >}} in the {{< Manual "Units" "Units" >}} page of the [[manual]].

* The following entry in the input file is not a variable that {{< octopus >}} will read directly, but rather illustrates the possibility of writing "user-defined" values and expressions to simplify the input file. In this case, we define the nitrogen  mass (<tt>Nitrogen_mass = 14.0</tt>) (note that in this case, as an exception, the value is expected to be in the so-called "atomic mass units", rather than in "atomic units"). This definition may be used elsewhere in the input file.

* The {{< Variable2 "Species" >}} block should contain the list of species that are present in the system to be studied. In this case we have only one species: nitrogen. The first field is a string that defines the name of the species, "N" in this case. The second field defines the type of species, in this case  <tt> species_pseudo </tt>. Then a list of parameters follows. The parameters are specified by a first field with the parameter name and the field that follows with the value of the parameter. Some parameters are specific to a certain species while others are accepted by all species. In our example <tt> set </tt> instructs {{< octopus >}} to use a pseudopotential for nitrogen from the <tt> standard </tt> set. This happens to be a Troullier-Martins pseudopotential defined in the {{< file "N.psf" >}} file found in the directory {{< file "share/octopus/pseudopotentials/PSF" >}}. Then come maximum <tt> lmax </tt> - component of the pseudopotential to consider in the calculation, and the <tt> lloc </tt> - component to consider as local. Generally, you want to set the maximum ''l'' to the highest available in the pseudopotential and the local ''l'' equal to the maximum ''l''. Finally, the mass of the species can also be modified from the default values by setting <tt> mass </tt> parameter.

* <tt>{{< Variable2 "XYZCoordinates" >}} = 'N.xyz'</tt>: The geometry of the molecule (in this case, a single atom in the grid origin) is described in this case in a file with the well known <tt>XYZ</tt> format. The file for this outrageously simple case is given by:

```text
 1
 This is a comment line
 N 0 0 0
```

* <tt>{{< Variable2 "ExtraStates" >}} = 1</tt>: By default, {{< octopus >}} performs spin-unpolarized calculations (restricted closed-shell, in Hartree-Fock terminology). It then places two electrons in each orbital. The number of orbitals, or Kohn-Sham states, is then calculated by counting the number of valence electrons present in the system, and dividing by two. In this case, since we have five valence electrons, the code would use three orbitals. However, we know beforehand that the HOMO orbital has a three-fold degeneracy, and as a consequence we need to put each one of the three <i>p</i> electrons in a different orbital. We therefore need one more orbital, which we get with this line in the input file.

* <tt>%{{< Variable2 "Occupations" >}}</tt> block: Generally, the occupations of the Kohn-Sham orbitals are automatically decided by the code, filling the lowest-energy orbitals. However, if we have degeneracies in the LUMO as in this case, the user may want to accommodate the electrons in a certain predefined way. In this example, the obvious way to fill the orbitals of the nitrogen atom is to put two electrons in the first and deepest orbital (the <i>s</i> orbital), and then one electron on each of the second, third and fourth orbitals (the <i>p</i> orbitals, which should be degenerate).

* <tt>{{< Variable2 "BoxShape" >}} = sphere</tt>: This is the choice of the shape of the simulation box, which in this case is set to be a sphere (other possible choices are <tt>minimum</tt>, <tt>cylinder</tt>, or <tt>parallelepiped</tt>).

* <tt>{{< Variable2 "Radius" >}} = 5.0*angstrom</tt>: The radius of the sphere that defines the simulation box.

* <tt>{{< Variable2 "Spacing" >}} = 0.18*angstrom</tt>: As you should know, {{< octopus >}} works in a real-space regular cubic mesh. This variable defines the spacing between points, a key numerical parameter, in some ways equivalent to the energy cutoff in plane-wave calculations.

##  Output  

Once you have constructed the input file and created the {{< file "N.xyz" >}} file, you may unleash {{< octopus >}} on it. Lets now go over some of the sections of the output.

####  Species  
```text

****************************** Species *******************************
  Species 'N'
    type             : pseudopotential
    file             : '/opt/share/octopus/pseudopotentials/PSF/N.psf'
    file format      : PSF
    valence charge   : 5.0
    atomic number    :   7
    form on file     : semilocal
    orbital origin   : calculated
    lmax             : 1
    llocal           : 0
    projectors per l : 1
    total projectors : 1
    application form : kleinman-bylander
    orbitals         : 16
    bound orbitals   : 16

**********************************************************************
```
</pre>
Here the code searches for the needed pseudopotential files, and informs the user about its success or failure. In this case, only the {{< file "N.psf" >}} file is required. Once that file has been processed, some information about it is written to the output. One of the most important pieces of information to be found here is the valence charge, which tells us how many electrons from this species will be considered in the calculation.

####  Grid  
```text

******************************** Grid ********************************
Simulation Box:
  Type = sphere
  Radius  [A] =   5.000
  Octopus will run in 3 dimension(s).
  Octopus will treat the system as periodic in 0 dimension(s).
Main mesh:
  Spacing [A] = ( 0.180, 0.180, 0.180)    volume/point [A^3] =      0.00583
  - inner mesh =      89727
  - total mesh =     127183
  Grid Cutoff [eV] =  1160.586810    Grid Cutoff [Ry] =    85.301565
**********************************************************************
```
</pre>
This step is about the construction of the mesh. As requested in the input file, a sphere of radius 5 Å is used, which contains a cubic regular real-space grid with spacing 0.18 Å. This implies 89727 points (<tt>inner mesh =  89727</tt>). For the sake of comparison with plane-wave-based codes, this is more or less equivalent to a plane-wave calculation that imposes a density cutoff of 1160.595 eV = 42.6 Hartree (except that in this case there is no artificial periodic repetition of the system).

####  Mixing  
```text

Input: [MixField = potential] (what to mix during SCF cycles)
Input: [MixingScheme = broyden]
```
</pre>
During the self-consistent procedure one has to use a [[Manual:Ground_State-Mixing|mixing scheme]] to help convergence. One can mix either the density or the potential, and there are several mixing schemes available.

####  Eigensolver  
```text

**************************** Eigensolver *****************************
Input: [Eigensolver = cg]
Input: [Preconditioner = pre_filter]
Input: [PreconditionerFilterFactor = 0.5000]
Input: [SubspaceDiagonalization = standard]
**********************************************************************
```
</pre>
Here we see that the {{< Manual "Ground_State-Eigensolver" "eigensolver" >}} used will be simple conjugate gradients (cg), and a preconditioner is used to speed up its convergence.

####  LCAO  
After some output you should see something like:
```text

Info: Performing initial LCAO calculation with      4 orbitals.
Info: Getting Hamiltonian matrix elements.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

Eigenvalues [eV]
 -st  Spin   Eigenvalue      Occupation
   1   --   -17.398866       2.000000
   2   --    -6.414013       1.000000
   3   --    -6.414013       1.000000
   4   --    -6.414013       1.000000
```
</pre>
This is the first step of a ground-state calculation: obtaining a reasonably good starting density and Kohn-Sham orbitals to feed in the self-consistent (SCF) procedure. For this purpose, {{< octopus >}} performs an initial calculation restricted to the basis set of atomic orbitals ([[Manual:Ground_State-LCAO|Linear Combination of Atomic Orbitals]], LCAO). The resulting eigenvalues of this calculation are written to standard output.

####  Wavefunction kind  

```text
 Info: SCF using real wavefunctions.
```

Very often one can work with real wave-functions. This is particularly helpful as calculations with real wave-functions are much faster than with complex ones. However, if a magnetic field is present, if the system is periodic, or if spin-orbit coupling is present, complex wave-functions are mandatory. But don't worry: the program is able to figure out by itself what to use.

####  SCF  
```text

*********************** SCF CYCLE ITER -    1 ************************
 etot  = -2.54485949E+02 abs_ev   =  4.52E-01 rel_ev   =  8.29E-03
 ediff =        2.69E-01 abs_dens =  2.82E-01 rel_dens =  5.65E-02
Matrix vector products:    108
Converged eigenvectors:      0

-  State  Eigenvalue [eV]  Occupation    Error
      1      -17.468243    2.000000   (5.8E-04)
      2       -6.518418    1.000000   (1.3E-04)
      3       -6.518418    1.000000   (1.3E-04)
      4       -6.518418    1.000000   (1.3E-04)

Density of states:

----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
---------------------------------------------------------------------%
---------------------------------------------------------------------%
%--------------------------------------------------------------------%
                                                                     ^


Elapsed time for SCF step     1:          1.11
**********************************************************************
```
</pre>

Now the SCF cycle starts. For every step, {{< octopus >}} outputs several pieces of information:

* The values <tt>abs_dens</tt> and <tt>rel_dens</tt> are to monitor the absolute and relative convergence of the density, while <tt>rel_ev</tt> and <tt>abs_ev</tt> are two alternative measures of the convergence, based on measuring the difference between input and output eigenvalues. The SCF procedure, by default, is stopped when <tt>rel_dens</tt> is smaller than $10^{-5}$. This may be altered with the appropriate input variables (see in the manual the variables {{< Variable2 "ConvAbsDens" >}}, {{< Variable2 "ConvRelDens" >}}, {{< Variable2 "ConvAbsEv" >}} and {{< Variable2 "ConvRelEv" >}}).

* The line <tt>Matrix vector products:    108</tt> tells us that the Hamiltonian was applied 108 times. This gives us an idea of the computational cost.

* The line <tt>Converged eigenvectors:      0</tt> tells us that upon completion of the diagonalization procedure, none of the orbitals met the required precision criterion for the wavefunctions. In a following example, we will modify this criterion in the input file.

* The list of eigenvalues is then printed, along with their errors: how much they deviate from "exact" eigenvalues of the current Hamiltonian. This number is the so-called "residue".

You can now take a look at the file {{< file "static/info" >}} that will hold a summary of the calculation.

##  Restarting  

Any ground-state calculation may be restarted later (to refine it if it did not converge properly, or with any other purpose), provided that the contents of the <tt>restart</tt> directory are preserved. You can try this now, just by running {{< octopus >}} again. You will notice that {{< octopus >}} did not give any warning after the line

```text
 Info: Loading restart information.
```

This is useful if you change slightly the parameters of the simulation (for example the XC functional or the convergence criteria). If you change the grid parameters {{< octopus >}} will not be able to restart from the previous calculation. If you do not want {{< octopus >}} to try to restart a calculation, you can set the variable {{< Variable2 "FromScratch" >}}.

In case you ware wondering what the restart information looks like, you can have a look at the contents of the {{< file "restart" >}} directory. This is where the files needed to restart a calculation are stored. It may contain several sub-directories depending on the calculations previously performed. In this case, it just contains one:

```text
 % ls restart
 gs
 % ls restart/gs
 0000000001.obf  0000000003.obf  density      df_010101.obf  df_010103.obf  dv_010102.obf  f_old_0101.obf  lxyz.obf  mixing  states  vhxc.obf          wfns
 0000000002.obf  0000000004.obf  density.obf  df_010102.obf  dv_010101.obf  dv_010103.obf  grid            mesh      occs    vhxc    vin_old_0101.obf
```

{{< octopus >}} stores each individual state in a different binary (yet platform-independent) file. In this case, we only have four states (files {{< file "0000000001.obf" >}} to {{< file "0000000004.obf" >}}). Some other useful quantities, like the density, are also stored in binary form. The other files are text files that contain diverse control information. It is unlikely that you will ever have to work directly with these files, but you may take a look around if you are curious. 


<!-- {{< figure src="/images/Nitrogen_pi02.jpg" width="500px" caption="pi orbital of N" >}}
{{< figure src="/images/Nitrogen_pi03.jpg" width="500px" caption="pi orbital of N" >}}
{{< figure src="/images/Nitrogen_pi04.jpg" width="500px" caption="pi orbital of N" >}}

Also, please add the following three lines to the {{< file "inp" >}} file:

{{< Variable2 "Output" >}} = wfs
{{< Variable2 "OutputWfsNumber" >}} = "1-4"
{{< Variable2 "OutputHow" >}} = dx + axis_x + axis_y + axis_z

You will now notice that the convergence procedure is stopped in only one iteration; the reason is that the starting point was now the (already converged) output state of the previous run. This is why the standard output has the following lines:

```text
 ******************** Loading restart information *********************
 All the needed files were succesfully read.
 **********************************************************************
```

Now take a look at the {{< file "static" >}} directory. Besides the {{< file "info" >}} file there are a bunch of new files, called {{< file "wf-001-00x-1.dx" >}}, and {{< file "wf-001-00x-1.a&-x3d;0,b&-x3d;0" >}}, where <tt>x</tt> runs from 1 to 4 (the four KS states), and where <tt>a</tt> and <tt>b</tt> are either <tt>x</tt>, <tt>y</tt> or <tt>z</tt>. Instructing {{< octopus >}} to generate these files was the task of the input variables that you just added. These files contain various functions related to the system (in this case, wavefunctions) for their visualization.


* <tt>{{< Variable2 "Output" >}} = wfs</tt>: This variable asks {{< octopus >}} to print the wavefunctions. One may also wish to print densities, potentials, etc. In the manual you may find the corresponding variables (<tt>density</tt>, etc).

* <tt>{{< Variable2 "OutputWfsNumber" >}} = '1-4'</tt>: This variable specifies which wavefunctions to print: <tt>'1-4'</tt> asks for the KS orbitals running from one to four; <tt>'1,2,7-10'</tt> would ask for the first two, and the ones from seven to ten.

* <tt>{{< Variable2 "Output" >}} = dx + axis_x + axis_y + axis_z</tt> This variable tells the code to print the functions in the format native to the OpenDX program. This program permits "sophisticated" data visualization (isosurfaces, contour plots, etc). In the figure you may see isosurfaces of the obtained <tt>p_x</tt>, <tt>p_y</tt> and <tt>p_z</tt> Kohn-Sham orbitals of Nitrogen. You will also find the files {{< file "wf-001-00x.x&-x3d;0,y&-x3d;0" >}}, {{< file "wf-001-000x.x&-x3d;0,z&-x3d;0" >}} and {{< file "wf-001-000x.y&-x3d;0,z&-x3d;0" >}}. These files contain the values of the wavefunctions along the $z$, $y$ and $x$, respectively, in the format coordinate, real value and imaginary value. This way you can plot the functions with easier-to-use plotting programs.
-->

{{Tutorial_foot|series=Octopus basics|prev=Getting started|next=Total energy convergence}}







---------------------------------------------
