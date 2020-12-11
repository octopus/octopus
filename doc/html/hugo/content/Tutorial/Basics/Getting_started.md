---
title: "Getting started"
tags: ["Basic", "Ground State", "Molecule", "Pseudopotentials", "DFT", "Total Energy"]
series: "Tutorial"
---


The objective of this tutorial is to give a basic idea of how {{< octopus >}} works.

###  Generating the input file  

With a text editor, create a text file called {{< file "inp" >}} containing the following text:

```text
 {{< Variable2 "CalculationMode" >}} = gs
 
 %{{< Variable2 "Coordinates" >}}
  'H' | 0 | 0 | 0
 %
```

This is the simplest example of an {{< octopus >}} input file:

* {{< Variable2 "CalculationMode" >}}<code> = gs</code>: This variable defines the run mode -- please consult the manual for the full list of the possible run modes. In this case we set it to <tt>gs</tt>, which instructs the code to start a ground-state calculation.

*  <tt>%{{< Variable2 "Coordinates" >}}</tt>: The entry is not just the definition of a variable, but rather of a full set of them -- a "block" of variables. The beginning of a block is marked by the <tt>%identifier</tt> line, and ended by a <tt>%</tt> line. In this case the identifier is <tt>%{{< Variable2 "Coordinates" >}}</tt>, where we list the atoms or species in our calculation and its coordinates, one per line. In this case, we put a single hydrogen atom in the center of our simulation box. 

The reason this input file can be so simple is that {{< octopus >}} comes with default values for the simulation parameters, and a set of default pseudopotentials for several elements (for properly converged calculations you might need to adjust these parameters, though).

To get a general idea of the format of the {{< octopus >}} input file, go and read the page about the [[Manual:Input file|Input file]] in the manual.

The documentation for each input variable can be found in the [http://octopus-code.org/doc/{{< octopus_version >}}/html/vars.php variable reference] online, and can also be accessed via the [[Manual:External utilities:oct-help|{{< code "oct-help" >}}]] utility.

###  Running Octopus  

Once you have written your input file, run the {{< command "octopus" >}} command (using {{< command "mpirun" >}} and perhaps a job script if you are using the parallel version). If everything goes correctly, you should see several lines of output in the terminal (if you don't, there must be a problem with your installation). As this is probably the first time you run {{< octopus >}}, we will examine the most important parts of the output:

* First there is an octopus drawn in ASCII art, the copyright notice and some information about the octopus version you are using and the system where you are running:

```text

      <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                                ___
                             .-'   `'.
                            /         \
                            |         ;
                            |         |           ___.--,
                   _.._     |0) ~ (0) |    _.---'`__.-( (_.
            __.--'`_.. '.__.\    '--. \_.-' ,.--'`     `""`
           ( ,.--'`   ',__ /./;   ;, '.__.'`    __
           _`) )  .---.__.' / |   |\   \__..--""  """--.,_
          `---' .'.''-._.-'`_./  /\ '.  \ _.-~~~````~~~-._`-.__.'
                | |  .' _.-' |  |  \  \  '.               `~---`
                 \ \/ .'     \  \   '. '-._)
                  \/ /        \  \    `=.__`~-.
             jgs  / /\         `) )    / / `"".`\
            , _.-'.'\ \        / /    ( (     / /
             `--~`   ) )    .-'.'      '.'.  | (
                    (/`    ( (`          ) )  '-;
                     `      '-;         (-'

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA

    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

                           Running octopus

Version                : 8.1
Commit                 : 9f1e56678d2698471df70e38f7dd4e37fd3d9044
Build time             : Thu Jul 19 14:51:54 CEST 2018
Configuration options  : max-dim=4 sse2 avx
Optional libraries     : arpack berkeleygw etsf_io gdlib netcdf pspio sparskit nlopt
Architecture           : x86_64
C compiler             : gcc
C compiler flags       : -g -Wall -O2 -march=native
Fortran compiler       : gfortran (GCC version 6.4.0)
Fortran compiler flags : -g -Wall -ffree-line-length-none -O2 -march=native

             The octopus is swimming in fenugreek (Linux)


            Calculation started on 2018/07/19 at 15:04:27

```
</pre>
Note that it also gives you the revision number, the compiler, and the compiler flags used. You should always include this information when submitting a bug report!

* The type of calculation it was asked to perform:
```text

************************** Calculation Mode **************************
Input: [CalculationMode = gs]
**********************************************************************
```
</pre>

* The species and pseudopotentials it is using:
```text

****************************** Species *******************************
  Species 'H'
    type             : pseudopotential
    file             : '/opt/local/share/octopus/pseudopotentials/PSF/H.psf'
    file format      : PSF
    valence charge   : 1.0
    atomic number    :   1
    form on file     : semilocal
    orbital origin   : calculated
    lmax             : 0
    llocal           : 0
    projectors per l : 1
    total projectors : 0
    application form : local
    orbitals         : 16
    bound orbitals   :  1
 
**********************************************************************
```
</pre>

* After some other output, {{< octopus >}} prints information about the grid: as we didn't say anything in the input file, {{< octopus >}} used the parameters recommended for this pseupopotential:
```text

******************************** Grid ********************************
Simulation Box:
  Type = minimum
  Species =     H     Radius =   7.559 b
  Octopus will run in 3 dimension(s).
  Octopus will treat the system as periodic in 0 dimension(s).
Main mesh:
  Spacing [b] = ( 0.435, 0.435, 0.435)    volume/point [b^3] =  0.08210
  - inner mesh =      22119
  - total mesh =      37759
 Grid Cutoff [H] =    26.123439    Grid Cutoff [Ry] =    52.246878
**********************************************************************
```
</pre>

* The level of theory and, in the case of (TD)DFT, the approximation to the exchange-correlation term:
```text

**************************** Theory Level ****************************
Input: [TheoryLevel = dft]

Exchange-correlation:
  Exchange
    Slater exchange (LDA)
    [1] P. A. M. Dirac, Math. Proc. Cambridge Philos. Soc. 26, 376 (1930)
    [2] F. Bloch, Z. Phys. 57, 545 (1929)
  Correlation
    Perdew & Zunger (Modified) (LDA)
    [1] J. P. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981), modified to improve the matching between the low- and high-rs

Input: [SICCorrection = sic_none]
**********************************************************************
```
</pre>

* At this point, {{< octopus >}} tries to read the wave-functions from a previous calculation. As there are none, it will give a warning.
```text

** Warning:
**   Could not find 'restart/gs' directory for restart.
**   No restart information will be read.

** Warning:
**   Unable to read wavefunctions.
**   Starting from scratch!
```
</pre>

* Now {{< octopus >}} commences the calculation. To get a reasonable starting point for the DFT calculation, the initial wavefunctions are calculated as a [[Manual:Ground_State-LCAO|Linear Combination of Atomic Orbitals]] (LCAO).
```text

Info: Performing initial LCAO calculation with      1 orbitals.
Info: Getting Hamiltonian matrix elements. 
ETA: .......1......2.......3......4......5.......6......7.......8......9......0
 
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation
   1   --    -0.233314       1.000000
Info: Ground-state restart information will be written
```
</pre>
```text
 
```
* After the LCAO, the real DFT calculation starts. For each self-consistency step some information is printed. When SCF {{< Manual "Ground_State-Convergence" "converges" >}}, the calculation is done.
```text

*********************** SCF CYCLE ITER -    1 ************************
 etot  = -4.48042396E-01 abs_ev   =  1.09E-03 rel_ev   =  4.66E-03
 ediff =       -2.76E-03 abs_dens =  8.90E-03 rel_dens =  8.90E-03
Matrix vector products:     27
Converged eigenvectors:      0

-  State  Eigenvalue [H]  Occupation    Error
      1       -0.234405    1.000000   (6.0E-05)

Elapsed time for SCF step     1:          0.12
**********************************************************************
```
</pre>
...
```text

*********************** SCF CYCLE ITER -    5 ************************
 etot  = -4.46377047E-01 abs_ev   =  3.50E-06 rel_ev   =  1.50E-05
 ediff =       -4.24E-06 abs_dens =  4.33E-06 rel_dens =  4.33E-06
Matrix vector products:      7
Converged eigenvectors:      1

-  State  Eigenvalue [H]  Occupation    Error
      1       -0.233013    1.000000   (9.6E-07)

Elapsed time for SCF step     5:          0.04
**********************************************************************


             Info: Writing states. 2018/07/19 at 15:13:58


        Info: Finished writing states. 2018/07/19 at 15:13:58

Info: SCF converged in    5 iterations

Info: Finished writing information to 'restart/gs'.

             Calculation ended on 2018/07/19 at 15:13:58

                          Walltime:  01. 10s

Octopus emitted 2 warnings.
```
</pre>
```text
 
```
Just running the command {{< command "octopus" >}} will write the output directly to the terminal. To have a saved copy of the output, it is generally advisable to redirect the output into a file, and to capture the standard error stream as well, which can be done like this: {{< command "octopus &> log" >}}. That would create a file called {{< file "log" >}} containing all output including warnings and errors in their context.

###  Analyzing the results  

After finishing the calculation you will find a series of files in the directory you ran:

```text
 % ls
 '''exec'''  inp  '''restart'''  '''static'''
```

For the moment we will ignore the '''exec'''  and  '''restart''' directories and focus on the {{< file "static/info" >}} file, which contains the detailed results of the ground-state calculation. If you open that file, first you will see some parameters of the calculations (that we already got from the output) and then the calculated energies and eigenvalues in Hartrees:

```text
 Eigenvalues [H]
  -st  Spin   Eigenvalue      Occupation
    1   --    -0.233013       1.000000
 
 Energy [H]:
       Total       =        -0.44637705
       Free        =        -0.44637705
       -----------
       Ion-ion     =         0.00000000
       Eigenvalues =        -0.23301327
       Hartree     =         0.28415332
       Int[n*v_xc] =        -0.30429841
       Exchange    =        -0.19375604
       Correlation =        -0.03975282
       vanderWaals =         0.00000000
       Delta XC    =         0.00000000
       Entropy     =         1.38629436
       -TS         =        -0.00000000
       Kinetic     =         0.41780616
       External    =        -0.91483022
       Non-local   =         0.00000000
```


Since by default {{< octopus >}} does a spin-unpolarized density-functional-theory calculation with the local-density approximation, our results differ from the exact total energy of 0.5 H. Our exchange-correlation functional can be set by the variable {{< Variable2 "XCFunctional" >}}, using the set provided by the [http://www.tddft.org/programs/Libxc Libxc] library.

###  Extra  

If you want to improve the LDA results, you can try to repeat the calculation with spin-polarization:

```text
 {{< Variable2 "SpinComponents" >}} = spin_polarized
```

And if you want to obtain the exact Schödinger equation result (something possible only for very simple systems like this one) you have to remove the self-interaction error (a problem of the LDA). Since we only have one electron the simplest way to do it for this case is to use independent electrons:
```text
 
 {{< Variable2 "TheoryLevel" >}} = independent_particles
```

A more general way would be to include self-interaction correction.

{{Tutorial_foot|series=Octopus basics|prev=|next=Basic input options}}







---------------------------------------------
