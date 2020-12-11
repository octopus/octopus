---
title: "oct-display partitions"
series: "Manual"
---


### NAME 
oct-display_partitions - Graphical representation of the mesh partitions.

### SYNOPSIS 
{{< command "oct-display_partitions" >}}

[oct-display_partitions does not read the standard input: all standard
input will be simply ignored. Also, oct-display_partitions accepts no
command line arguments.]

### DESCRIPTION 
This program is one of the octopus utilities.

It is a script to use gnuplot to plot the partitioning of the mesh. The
files {{< TODO >}} {{< file "mesh_partition.???" >}} as produced in {{< file "debug/mesh_partition" >}} by the
{{< octopus >}} debug mode have to be present in the current working directory
when this script is invoked. The plot can be found in the file
{{< file "mesh_partitions.png" >}}. This script generates a gnuplot-script called
{{< file "mesh_partitions_index.gp" >}} which is stored in the current working
directory and can be loaded into gnuplot manually. With:

{{< command_line "gnuplot> load mesh_partitions_index.gp" >}}
{{< command_line "gnuplot> set term x11" >}}
{{< command_line "gnuplot> replot" >}}

the plot can be reproduced and shown on the screen so that rotating and
zooming is possible.

{{< manual_foot prev="Manual:External utilities:oct-dielectric-function" next="Manual:External utilities:oct-harmonic-spectrum" >}}
---------------------------------------------
