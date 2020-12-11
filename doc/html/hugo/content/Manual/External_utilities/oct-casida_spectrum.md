---
title: "oct-casida spectrum"
series: "Manual"
---


### NAME 
oct-casida_spectrum - Broadens the linear-response spectra from the Casida run mode.

### SYNOPSIS 
{{< command "oct-casida_spectrum" >}}

[oct-casida_spectrum does not read standard input: all standard input will be
simply ignored. An input file named {{< file "inp" >}} must be present in the running
directory. Also, oct-casida_spectrum accepts no command line arguments, since
there is not a standard way to do this with Fortran 90.]

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities.

It post-processes the files {{< file "casida/casida" >}}, {{< file "casida/petersilka" >}} and
{{< file "casida/eps_diff" >}}, that contain the excitation energies and oscillator strengths.

The parameters of the spectrum can be set using the variables {{< Variable2 "CasidaSpectrumBroadening" >}}, {{< Variable2 "CasidaSpectrumMinEnergy" >}}, {{< Variable2 "CasidaSpectrumMaxEnergy" >}}, and {{< Variable2 "CasidaSpectrumEnergyStep" >}}.

{{< manual_foot prev="Manual:External utilities:oct-atomic_occupations" next="Manual:External utilities:oct-center-geom" >}}
---------------------------------------------
