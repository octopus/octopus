---
title: "Troubleshooting"
series: "Manual"
---


If {{< octopus >}} works properly on your system (i.e. you can recreate the results in the [[tutorials]]) but you have troubles using it for your own work, here are some things to try.  Please feel free to add your own ideas here.

#### Read parser.log 
If you look at the file {{< file "exec/parser.log" >}}, it will tell you the value of the variables that you set with the {{< file "inp" >}} file, as well as all the variables which are taking their default value.  This can sometimes be helpful in understanding the behavior of the program.

#### Use OutputDuringSCF 
If you add {{< Variable2 "OutputDuringSCF" >}} = {{< code_line "yes" >}} to your input file, you can examine the results of each iteration in the Self Consistent Field calculation.  So if you also have the variable {{< Variable2 "Output" >}} set to {{< code_line "Output" >}} = {{< code_line "density + potential" >}}, both the electron density and the Kohn-Sham, bare, exchange-correlation and Hartree potentials will be written to a folder (called, e.g., {{< code_line "scf.0001" >}}) after each SCF iteration.

#### Set Debug 
Set the variable {{< Variable2 "Debug" >}} to 1 for some extra diagnostic info and a stack trace with any fatal error, 2 to add a full stack strace, and 99 to get a stack trace from each MPI task when running in parallel.

{{< manual_foot prev="Manual:Symmetry" next="Manual:Ground State" >}}
---------------------------------------------
