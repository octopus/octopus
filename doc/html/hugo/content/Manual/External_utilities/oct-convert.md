---
title: "oct-convert"
series: "Manual"
---


### Name 
oct-convert - Octopus utility to read obj files and write in many different file formats

### Description 

This executable gives the ability to read files written during the ground-state or time-dependent execution

### Example 

You can run ground-state and time-dependent execution of the [benzene example](../Benzene_molecule).

Then, we have to add this to the inp file, if we want to have the ground state density in DX format:

```text
 {{< Variable2 "Output" >}} = density
 {{< Variable2 "OutputFormat" >}} = dx
 {{< Variable2 "ConvertFolder" >}} = 'restart/gs'
 {{< Variable2 "ConvertFilename" >}} = 'density'
 {{< Variable2 "ConvertIterateFolder" >}} = no
```

To convert the restart wave-functions (from 1 to 10) of a td run:

```text
 {{< Variable2 "Output" >}} = density
 {{< Variable2 "OutputFormat" >}} = dx
 {{< Variable2 "ConvertIterateFolder" >}} = no
 {{< Variable2 "ConvertFilename" >}} = ' '
 {{< Variable2 "ConvertStart" >}} = 1
 {{< Variable2 "ConvertEnd" >}}   = 10
 {{< Variable2 "ConvertFolder" >}} = 'restart/td/'
```

If we want to convert the densities of the time-dependent executions, from files td.0000001 to td.0000010:

```text
 {{< Variable2 "Output" >}} = density
 {{< Variable2 "OutputFormat" >}} = dx
 {{< Variable2 "ConvertIterateFolder" >}} = yes
 {{< Variable2 "ConvertStart" >}} = 1
 {{< Variable2 "ConvertEnd" >}}   = 10
```


{{< manual_foot prev="Manual:External utilities:oct-conductivity" next="Manual:External utilities:oct-dielectric-function" >}}
---------------------------------------------
