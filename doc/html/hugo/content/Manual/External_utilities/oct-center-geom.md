---
title: "oct-center-geom"
series: "Manual"
---


### NAME 
oct-center-geom - Centers a molecule's geometry

### SYNOPSIS 
{{< command "oct-center-geom" >}}

[oct-center-geom does not read the standard input: all standard input
will be simply ignored. An input file named {{< file "inp" >}} must be present in the
running directory. Also, oct-center-geom accepts no command-line
arguments, since there is not a standard way to do this with Fortran
90.]

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities.

It reads the coordinates defined in the {{< file "inp" >}} file, and constructs an output xyz file, that will be called
{{< file "adjusted.xyz" >}} file, that describes the same system but in which the
atomic coordinates are centered, and (optionally) has the axes aligned.

To control the orientation of the centered molecule there are two parameters: {{< Variable2 "MainAxis" >}} and {{< Variable2 "AxisType" >}}.

Do not forget then to change your input file to use this file instead of your old geometry (by changing {{< Variable2 "XYZCoordinates" >}}={{< value "'adjusted.xyz'" >}}). 

Be careful with units, this utility honours the {{< Variable2 "Units" >}}, {{< Variable2 "UnitsInput" >}} and {{< Variable2 "UnitsOutput" >}} variables, so the {{< file "adjusted.xyz" >}} file will be in the specified output units.

{{< manual_foot prev="Manual:External utilities:oct-casida_spectrum" next="Manual:External utilities:oct-check_deallocs" >}}
---------------------------------------------
