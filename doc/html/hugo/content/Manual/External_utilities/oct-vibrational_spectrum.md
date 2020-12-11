---
title: "oct-vibrational spectrum"
series: "Manual"
---


This utility calculates the vibrational spectrum from a molecular-dynamics run. 

What this utility does is to read the velocity from the {{< file "td.general/coordinates" >}} and calculate the Velocity Autocorrelation Function:

$$
C\_{v}(t)=\\sum\_{i=1}^{N\_{atoms}}\\vec{v}\_i(t)\\cdot\\vec{v}\_i(t\_0)\\ ,
$$

afterward a cosinusoidal envelope is added, to make the periodic extension of the function continuous and then the spectrum is calculated by taking the Fourier transform of the function. On exit, two files are generated {{< file "td.general/velocity_autocorrelation" >}} and {{< file "td.general/vibrational" >}}.

This utility honours the variables {{< Variable2 "PropagationSpectrumStartTime" >}} and {{< Variable2 "PropagationSpectrumEndTime" >}} to control the time of sampling. Note that the velocity in the initial time must be different from zero, or $C_{v}$ will be identically zero.

As a discrete Fourier tranform is used, this utility can take several minutes to process a large run. If {{< octopus >}} was compiled with OpenMP support, this utility can run in several threads.

{{< manual_foot prev="Manual:External utilities:oct-vdW_c6" next="Manual:External utilities:oct-xyz-anim" >}}
---------------------------------------------
