---
title: "Visualization"
series: "Manual"
Weight: 102
---


Every given number of time iterations, or after ground-state calculations, some of the functions that characterise the system may be written to disk so that they may be analized. Files are written within {{< file "static/" >}} output directory after the self-consistent field, or within {{< file "td.x/" >}} directories, during evolution, where “x” stands for the iteration number at which each write is done. 

The function that you want to plot is selected by the {{< Variable2 "Output" >}} variable and the output format is chosen by the 
{{< Variable2 "OutputFormat" >}}.

###  dx  

This is an OpenDX network, aimed at the visualization of wavefunctions. To be able to use it, you need to have properly installed the [http://opendx.org OpenDX] program, as well as the Chemistry extensions developed at the Cornell Theory Center. Please take a look [[Releases-OpenDX|here]] to see how to obtain and install this software. Unfortunately, since this software is old and unmaintained, you are likely to have trouble installing.

Once these are working, you may follow the [tutorial for the benzene molecule](../Benzene_molecule).

###  XCrySDen  

Atomic coordinates (finite or periodic), forces, and functions on a grid can be plotted with the free program [http://www.xcrysden.org/ XCrySDen]. Its XSF format also can be read by [http://inac.cea.fr/sp2m/L_Sim/V_Sim/index.en.html V_Sim] and [http://www.geocities.jp/kmo_mma/crystal/en/vesta.html Vesta]. Beware, these all probably assume that your output is in Angstrom units (according to the [http://www.xcrysden.org/doc/XSF.html specification]), so use UnitsOutput = eV_Angstrom, or your data will be misinterpreted by the visualization software.

###  CUBE  

The Gaussian cube format (see http://gaussian.com/cubegen/, http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html and https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html for a more detailed description of the format) can be output, and can be read by VMD, XCrysDen, Avogadro, and other software. Note that CUBE files are always in atomic units, so the UnitsOutput input option will be ignored.

###  PDB  

Everything is supposed to be in Angstroms: http://deposit.rcsb.org/adit/docs/pdb_atom_format.html

###  XYZ  

Generally considered to be in Angstroms: http://openbabel.org/wiki/XYZ_%28format%29, https://en.wikipedia.org/wiki/XYZ_file_format, https://www.molpro.net/info/2012.1/doc/manual/node100.html, http://departments.icmab.es/leem/siesta/Documentation/Manuals/siesta-3.1-manual/node32.html

{{< manual_foot prev="Manual:Geometry Optimization" next="Manual:Advanced ways of running Octopus" >}}
---------------------------------------------
