---
Title: "Time propagation"
Weight: 1
---

In the new multisystem framework, the time propagation is handled by the routine {{< code "time_dependent_run_multisystem()" >}}

{{% expand "Expand for the source" %}}
```Fortran
#include_subroutine time_dependent_run_multisystem
```
{{% /expand %}}

The code is quite self-explanatory. 
The central part of the time propagation is the {{< code "do ... while" >}} loop ({{< code "! The full TD loop" >}}).
Here the function {{< code "system%dt_operation()" >}} is called successively. 
{{< notice note >}}
Note, that propagators in general will need several calls of this to finish one time step of the propagator!
{{< /notice >}} 

The {{< code "system%dt_operation()" >}} subroutine is handling all general tasks, such as starting and ending the SCF loop, or updating interaction, while
passing remaining operations down to the specific implementation of the systems.

{{% expand "Implementation of system_dt_operation" %}}
```Fortran
#include_subroutine system_dt_operation
```
{{% /expand %}}

