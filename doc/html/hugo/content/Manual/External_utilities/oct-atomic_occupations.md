---
title: "oct-atomic occupations"
series: "Manual"
---


This script prints out to standard output the electronic configuration of a given atom, and the way in which this electronic configuration should be reflected in the {{< Variable2 "Occupations" >}} block of a corresponding {{< octopus >}} input file.

### Options 

{{< flag "-s species" >}}
species should be the atomic symbol (e.g. Na, Au, etc).

{{< flag "-h" >}}
Show a brief summary of command line options.

### Examples 
{{< command_line "oct-atomic_occupations -s Na" >}}

{{< command_line "oct-atomic_occupations -s Ti_sc" >}}

```text

$ for x in \$(cat /usr/share/octopus/PP/defaults | awk '{print \$1}')
> do oct-atomic_occupations -s $x
> done
```
</pre>

{{< manual_foot prev="Manual:External utilities:oct-analyze_projections" next="Manual:External utilities:oct-casida_spectrum" >}}
---------------------------------------------
