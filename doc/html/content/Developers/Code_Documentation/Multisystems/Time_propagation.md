---
Title: "Time propagation"
Weight: 1
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


In the new multisystem framework, the time propagation is handled by the routine {{< code "time_dependent_run_multisystem()" >}}

{{% expand "Expand for the source" %}}
```Fortran
#include_subroutine time_dependent_run_multisystem
```
{{% /expand %}}

The code is quite self-explanatory. 
The central part of the time propagation is the {{< code "do ... while" >}} loop (see {{< code "\"! The full TD loop\"" >}}).
Here the function {{< code "system%dt_operation()" >}} is called successively. 
{{< notice note >}}
Note, that propagators in general will need several calls of this to finish one time step of the propagator!
{{< /notice >}} 

In the multisystem framework, the highest level system is of type {{< code "multisystem_basic_t" >}}, which inherits the {{< code "dt_operation()" >}} routine from 
the abstract {{< code "multisystem_t" >}} type.
  
{{% expand "Implementation of multisystem_dt_operation" %}}
```Fortran
#include_subroutine multisystem_dt_operation
```
{{% /expand %}}

THis routine first calls the general {{< code "system%dt_operation()" >}} subroutine from the parent class {{< code "system_t" >}} and then loops over all
systems, which are part of the multisystem, and calls theis specific {{< code "dt_operations()" >}} routine.

The general {{< code "system%dt_operation()" >}} subroutine is handling all general tasks, such as starting and ending the SCF loop, or updating interaction. 


The steps implemented here are:

* {{< code "UPDATE_INTRERACTIONS" >}}
  - increment the **propagator**-clock by one tick (if all interactions can be updated)
  - try to update all interactions (this might not always be possible, if some interaction partners are still lagging behind)
  - if all interactions are updated succesfully, progress the propagator to the next step.
* {{< code "FINISHED" >}}
  - increment the **system**-clock by one tick.
  - rewind the propagator
* {{< code "START_SCF_LOOP" >}} and {{< code "END_SCF_LOOP" >}}
  - handles self-consistency and convergence for predictor-corrector algorithms.
  - no clocks are updated.

All remaining operations are passed down to the system specific routine {{< code "system%do_dt_operation()" >}}.


{{% expand "Implementation of system%dt_operation()" %}}
```Fortran
#include_subroutine system_dt_operation
```
{{% /expand %}}

This routine is triggering the update of the interactions, via the call to {{< code "system%update_interations()" >}},
which loops over the interactions associated with the system, and all required exposed quantities.

{{% expand "Implementation of system%update_interations()" %}}
```Fortran
#include_function system_update_interactions
```
{{% /expand %}}


Other extensions of the {{< code "system_t" >}} class might overload that routine in order to implement their system-specific steps.
The routine {{< code "system%do_dt_operation()" >}} is where the actual algorithm operations for that system have to be implemented.
