---
title: "oct-xyz-anim"
series: "Manual"
---


### NAME 
oct-xyz-anim - Constructs an "animated" xyz file from a time-propagation file where the atoms were allowed to move

### SYNOPSIS 
{{< command "oct-xyz-anim" >}}

[oct-xyz-anim does not read the standard input: all standard input will
be simply ignored. An input file named {{< file "inp" >}} must be present in the
running directory. Also, oct-xyz-anim accepts no command line arguments,
since there is not a standard way to do this with Fortran 90.]

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities. It reads out the {{< file "td.general/coordinates" >}} file, and makes a movie in XYZ format, called 'movie.xyz'.

{{< manual_foot prev="Manual:External utilities:oct-vibrational_spectrum" next="Manual:Deprecated_Utilities" >}}
---------------------------------------------
