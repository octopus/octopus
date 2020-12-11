---
title: "oct-run periodic table"
series: "Manual"
---


### NAME 
oct-run_periodic_table - run octopus for atoms

### SYNOPSIS 
{{< command  " oct-run_periodic_table [ option ] ..." >}}

### DESCRIPTION 
This script is one of the octopus utilities.

Simple script to run octopus for all atoms in the periodic table for
which pseudopotentials are available.

{{< flag "-h" >}}
Show this message

{{< flag "-n" >}}
Run in dry run mode (show what would be executed)

{{< flag "-s species" >}}
Run octopus for given species.

{{< flag "-a" >}}
Run octopus for all atoms

{{< flag "-r" >}}
Do a restart, i.e. do not use {{< Variable2 "fromScratch" >}}=yes.

{{< flag "-t temperature" >}}
Run with electronic temperature.

{{< flag "-e no_states" >}}
Use extra states

{{< flag "-x octopus_executable" >}}

{{< flag "-d defaults_file_for_pseudopotentials" >}}

### EXAMPLES 
The following command would run the octopus executable to
calculate the He atom.

{{< command_line " oct-run_periodic_table -s He -x octopus" >}}

{{< manual_foot prev="Manual:External utilities:oct-propagation_spectrum" next="Manual:External utilities:oct-run_regression_test.pl" >}}
---------------------------------------------
