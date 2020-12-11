---
title: "Hello world"
series: "Manual"
---


As a first example, we will take a sodium atom. With your favourite text editor, create the file {{< file "inp" >}}.

```text
 {{< Variable2 "CalculationMode" >}} = gs
 %{{< Variable2 "Coordinates" >}}
     'Na' | 0.0 | 0.0 | 0.0 
 %
```

This input file should be essentially self-explanatory. 

Note that when a species is not specified in the {{< Variable2 "Species" >}} block, octopus reads the information of pseudopotentials from the {{< file "defaults" >}} file (located under {{< inst_file "share/octopus/PP/" >}}. This file also contains default values for {{< Variable2 "Radius" >}} and {{< Variable2 "Spacing" >}}.

Then run octopus – for example, do 

{{< command_line "octopus > out" >}}

so that the output is stored in {{< file "out" >}} file. If everything goes OK, {{< file "out" >}} should look like:
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

Version                : 5.0.1
Revision               : 15042
Build time             : Tue Jan 12 11:04:57 EST 2016
Configuration options  : max-dim=3 mpi sse2
Optional libraries     : arpack berkeleygw etsf_io gdlib metis mpi2 netcdf newuoa parmetis parpack pfft scalapack sparskit
Architecture           : x86_64
C compiler             : /opt/local/bin/mpicc-mpich-mp (/usr/bin/clang)
C compiler flags       :  -pipe -O3 -arch x86_64
Fortran compiler       : /opt/local/bin/mpif90-mpich-mp (/opt/local/bin/gfortran-mp-5)
Fortran compiler flags : -pipe -O3

  The octopus is swimming in dhcp-18-189-27-45.dyn.MIT.EDU (Darwin)


            Calculation started on 2016/01/13 at 20:24:58


************************** Calculation Mode **************************
Input: [CalculationMode = gs]
**********************************************************************

Reading Coordinates from Coordinates block

****************************** Species *******************************
Reading pseudopotential from file:
      '/opt/local/share/octopus/PP/PSF/Na.psf'
      Calculating atomic pseudo-eigenfunctions for species Na....
Info: l =  0 component used as local potential.
Info: l =  0 is maximum angular momentum considered.
Number of orbitals: total =     16, bound =      4
**********************************************************************


***************************** Symmetries *****************************
Symmetry elements : (i) (Cinf) (sigma)
Symmetry group    : Kh
**********************************************************************

Input: [SpinComponents = unpolarized]
Input: [SmearingFunction = semiconducting]
Input: [SymmetrizeDensity = no]

******************************* States *******************************
Total electronic charge  =        1.000
Number of states         =        1
States block-size        =        1
**********************************************************************

Info: Using default spacing(1) [b] =  0.567
Info: Using default spacing(2) [b] =  0.567
Info: Using default spacing(3) [b] =  0.567
Input: [CurvMethod = curv_uniform]
Input: [DerivativesStencil = stencil_star]

************************** Parallelization ***************************
Octopus will run in *serial*
**********************************************************************

Info: Generating weights for finite-difference discretization of x-gradient
Info: Generating weights for finite-difference discretization of y-gradient
Info: Generating weights for finite-difference discretization of z-gradient
Info: Generating weights for finite-difference discretization of Laplacian

******************************** Grid ********************************
Simulation Box:
  Type = minimum
  Species =    Na     Radius =  13.228 b
  Octopus will run in 3 dimension(s).
  Octopus will treat the system as periodic in 0 dimension(s).
Main mesh:
  Spacing [b] = ( 0.567, 0.567, 0.567)    volume/point [b^3] =      0.18221
  - inner mesh =      52971
  - total mesh =      79699
  Grid Cutoff [H] =    15.354165    Grid Cutoff [Ry] =    30.708329
**********************************************************************

Info: states-block size = 0.6 MiB
Input: [StatesOrthogonalization = gram_schmidt]

****************************** Hartree *******************************
The chosen Poisson solver is 'interpolating scaling functions'
**********************************************************************


**************************** Theory Level ****************************
Input: [TheoryLevel = dft]

Exchange-correlation:
  Exchange
    Slater exchange (LDA)
    [1] PAM Dirac, Proceedings of the Cambridge Philosophical Society 26, 376 (1930)
    [2] F Bloch, Zeitschrift fuer Physik 57, 545 (1929)
  Correlation
    Perdew & Zunger (Modified) (LDA)
    [1] Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)
    [2] Modified to improve the matching between the low- and high-rs parts

Input: [SICCorrection = sic_none]
**********************************************************************

Input: [FilterPotentials = filter_none]
Info: Pseudopotential for Na
  Radii for localized parts:
    local part     =  3.5 b
    non-local part =  0.0 b
    orbitals       = 19.9 b

Input: [RelativisticCorrection = non_relativistic]
Input: [AbsorbingBoundaries = not_absorbing]

****************** Approximate memory requirements *******************
Mesh
  global  :       1.5 MiB
  local   :       1.8 MiB
  total   :       3.4 MiB

States
  real    :       0.6 MiB (par_kpoints + par_states + par_domains)
  complex :       1.2 MiB (par_kpoints + par_states + par_domains)

**********************************************************************

Info: Generating external potential
      done.
Info: Octopus initialization completed.
Info: Starting calculation mode.
Info: Allocating ground state wave-functions
Info: Blocks of states
      Block       1 contains       1 states:       1 -       1
Info: Ground-state allocation done.

** Warning:
**   Could not find 'restart/gs' directory for restart.
**   No restart information will be read.


** Warning:
**   Unable to read wavefunctions.
**   Starting from scratch!

Input: [MixField = density] (what to mix during SCF cycles)
Input: [TypeOfMixing = broyden]

**************************** Eigensolver *****************************
Input: [Eigensolver = cg]
Input: [Preconditioner = pre_filter]
Input: [SubspaceDiagonalization = standard]
**********************************************************************

Input: [LCAOStart = lcao_full]
Input: [LCAOScaleFactor = 1.000]
Input: [LCAOMaximumOrbitalRadius = 20.00 b]
Info: Single-precision storage for     1 extra orbitals will be allocated.
Info: Unnormalized total charge =      0.999665
Info: Renormalized total charge =      1.000000
Info: Setting up Hamiltonian.
Info: Performing initial LCAO calculation with      2 orbitals.
Info: Getting Hamiltonian matrix elements.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation
   1   --    -0.103066       1.000000
Info: Ground-state restart information will be written to 'restart/gs'.
Info: SCF using real wavefunctions.
Info: Starting SCF iteration.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    1 ************************
 etot = -1.84238245E-01 abs_ev   =  3.88E-04 rel_ev   =  3.75E-03
                        abs_dens =  6.54E-03 rel_dens =  6.54E-03
Matrix vector products:     27
Converged eigenvectors:      0
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103455       1.000000      (5.5E-05)

Elapsed time for SCF step     1:          0.08
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    2 ************************
 etot = -1.84238838E-01 abs_ev   =  1.02E-04 rel_ev   =  9.88E-04
                        abs_dens =  4.55E-03 rel_dens =  4.55E-03
Matrix vector products:     27
Converged eigenvectors:      0
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103353       1.000000      (1.3E-06)

Elapsed time for SCF step     2:          0.10
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    3 ************************
 etot = -1.84239600E-01 abs_ev   =  2.23E-04 rel_ev   =  2.16E-03
                        abs_dens =  4.35E-04 rel_dens =  4.35E-04
Matrix vector products:     19
Converged eigenvectors:      1
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103130       1.000000      (9.3E-07)

Elapsed time for SCF step     3:          0.09
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    4 ************************
 etot = -1.84239592E-01 abs_ev   =  2.81E-07 rel_ev   =  2.73E-06
                        abs_dens =  6.62E-04 rel_dens =  6.62E-04
Matrix vector products:     12
Converged eigenvectors:      1
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103130       1.000000      (7.2E-07)

Elapsed time for SCF step     4:          0.07
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    5 ************************
 etot = -1.84239609E-01 abs_ev   =  6.42E-06 rel_ev   =  6.23E-05
                        abs_dens =  9.12E-06 rel_dens =  9.12E-06
Matrix vector products:     15
Converged eigenvectors:      1
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103123       1.000000      (8.7E-07)

Elapsed time for SCF step     5:          0.07
**********************************************************************


             Info: Writing states. 2016/01/13 at 20:24:59


        Info: Finished writing states. 2016/01/13 at 20:24:59

Info: SCF converged in    5 iterations

Info: Finished writing information to 'restart/gs'.

             Calculation ended on 2016/01/13 at 20:24:59

                          Walltime:  01.  6s

Octopus emitted 2 warnings.
```
</pre>

Take now a look at the working directory. Besides the initial file ({{< file "inp" >}}) and the {{< file "out" >}} file, three new directories appear. In {{< file "static/" >}}, you will find the file {{< file "info" >}}, with information about the static calculation (it should be hopefully self-explanatory, otherwise please complain to the authors...). In {{< file "restart/" >}}, you will find the {{< file "gs" >}} directory that contains restart information about the ground-state, which is used if, for example, you want to start a time-dependent calculation afterwards. Finally, the {{< file "exec" >}} directory has information about the run of octopus; inside the {{< file "parser.log" >}} contains all the input variables parsed by octopus.

###  Exercises  

* Study how the total energy and eigenvalue of the sodium atom improve with the mesh spacing.
* Calculate the static polarizability of the sodium atom ({{< Variable2 "CalculationMode" >}} = em_resp). A {{< file "em_resp/freq_0.0000/alpha" >}} will be created containing the static polarizability tensor. 
* Calculate a few unoccupied states ({{< Variable2 "CalculationMode" >}} = unocc). The eigenspectrum will be in the file {{< file "static/eigenvalues" >}}. Why don't we find a Rydberg series in the eigenspectrum?
* Repeat the previous calculation with different exchange and correlation functionals like PBE, LB94, and exact exchange (see {{< Variable2 "XCFunctional" >}}).
* Perform a time-dependent evolution ({{< Variable2 "CalculationMode" >}} = td), to calculate the optical spectrum of the Na atom. Use a {{< Variable2 "TDDeltaStrength" >}} = 0.05, polarised in the ''x''-direction. The multipole moments of the density are output to the file {{< file "td.general/multipoles" >}}. You can process this file with the utility [[Manual:External_utilities:oct-propagation_spectrum | {{< command "oct-propagation_spectrum" >}}]] to obtain the optical spectrum. If you have computer time to waste, re-run the time-dependent simulation for some other xc choices. 

{{< manual_foot prev="Manual:Deprecated Utilities" next="Manual:Examples:Benzene" >}}
---------------------------------------------
