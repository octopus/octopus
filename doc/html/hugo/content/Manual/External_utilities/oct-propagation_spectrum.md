---
title: "oct-propagation spectrum"
series: "Manual"
---


### NAME 
oct-propagate_spectrum - Calculates the absorption cross section tensor from the results of a time-propagation run.

### SYNOPSIS 
{{< command "oct-propagate_spectrum" >}}

[oct-propagate_spectrum does not read the standard input: all standard input
will be simply ignored. An input file named {{< file "inp" >}} must be present in the
running directory. Also, oct-propagate_spectrum accepts no command line
arguments, since there is not a standard way to do this with Fortran
90.]

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities.

This utility generates the dipole strength function of the given system. Its main input is the td.general/multipoles file. Output is written to a file called spectrum. This file is made of two columns: energy (in eV or a.u., depending on the units specified in the input file), and dipole strength function (in 1/eV, or 1/a.u., idem).

In the input file, the user may set the {{< Variable2 "PropagationSpectrumTransform" >}} (this should be set to “sine” for proper use), the {{< Variable2 "PropagationSpectrumDampMode" >}} (recommended value is “polynomial”, which ensures fulfilling of the N-sum rule), the {{< Variable2 "PropagationSpectrumStartTime" >}}, the {{< Variable2 "PropagationSpectrumEndTime" >}}, the {{< Variable2 "PropagationSpectrumEnergyStep" >}}, and the {{< Variable2 "PropagationSpectrumMaxEnergy" >}}.

{{< manual_foot prev="Manual:External utilities:oct-photoelectron_spectrum" next="Manual:External utilities:oct-run_periodic_table" >}}
---------------------------------------------
