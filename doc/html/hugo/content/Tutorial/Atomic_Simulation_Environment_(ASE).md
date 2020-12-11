---
title: "Atomic Simulation Environment (ASE)"
tags: ["Tutorial", "Advanced"]
series: "Tutorial"
---


##  Overview  

About ASE: https://wiki.fysik.dtu.dk/ase/about.html

Installation: https://wiki.fysik.dtu.dk/ase/install.html

The Octopus interface is distributed as part of the ASE and does not require any special configuration or recompilation of Octopus.  It can be used to set up and run Octopus calculations from a Python script, and these can be combined with different algorithms for structure optimization, molecular dynamics, nudged-elastic-band method for saddle-point searches, vibration/phonon analysis, genetic algorithm, and other features of ASE.

Also, the ASE contains functions to generate structures like crystals, surfaces, or nanoparticles, as well as many standard molecule geometries.

The ASE interface works by writing input files to the disk and running Octopus as an external process.  That means some steps like structure optimizations can be less efficient because they require disk I/O to restore quantities such as densities or wavefunctions.  This disadvantage can be circumvented by using a ram-disk, if the amount of memory is not too large of course.  One can also simply use ASE to script the generation of input files (documentation may be slightly inadequate).

##  Exercises  

See https://wiki.fysik.dtu.dk/ase/ase/calculators/octopus.html

<span class=noprint><hr>
Back to [[Tutorials]]



---------------------------------------------
