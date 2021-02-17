---
Title: "Time propagation"
Weight: 10
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


In the new multisystem framework, the time propagation is handled by the routine {{< code "time_dependent_run_multisystem()" >}}

### The main loop

{{% expand "Implementation of time_dependent_run_multisystem()" %}}
```Fortran
#include_subroutine time_dependent_run_multisystem
```
{{% /expand %}}

The code is quite self-explanatory. 
The central part of the time propagation is the {{< code "do ... while" >}} loop (see {{< code "\"! The full TD loop\"" >}}).
Here the function {{< code "system%dt_operation()" >}} is called successively, triggering one algorithmic step of the multisystem propagator.
{{< notice note >}}
The propagators in general will need several calls of this to finish one time step of the propagator!
{{< /notice >}} 
{{< notice info >}}
At the top level, there is only one system (in general called {{< code "." >}}), which is of type {{< code "multisystem_basic_t" >}}, which then contains other systems as subsystems. These subsystems can be multisystems themselves.
{{< /notice >}}
### Updating the system


In the multisystem framework, the highest level system is of type {{< code "multisystem_basic_t" >}}, which inherits the {{< code "dt_operation()" >}} routine from 
the abstract {{< code "multisystem_t" >}} type.  

{{% expand "Implementation of multisystem_dt_operation" %}}
```Fortran
#include_subroutine multisystem_dt_operation
```
{{% /expand %}}

This routine first calls the general {{< code "system_dt_operation()" >}} subroutine from the parent class {{< code "system_t" >}} and then loops over all
subsystems, which are part of the multisystem, and calls their specific {{< code "system%dt_operations()" >}} routine.

The general {{< code "system_dt_operation()" >}} subroutine is handling all general algorithmic operations, such as starting and ending the SCF loop, or updating interaction. 

{{% expand "Implementation of system_dt_operation()" %}}
```Fortran
#include_subroutine system_dt_operation
```
{{% /expand %}}

This routine fetches the next operation from the propagator and executes it.
The following operations are implemented here:

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

All remaining operations are passed down to the system specific routine {{< code "system%do_dt_operation()" >}}, which is the place where specific systems need to implement all algorithmic steps, required by the propagators, allowed for that system.

{{% notice note %}}
In a multisystem, the multisystem itself as well as all its subsystems are propagating according to their _own_ propagators.
In principle, they could be different, which could lead to unexpected behaviour. If propagators of the subsystems are not specified explicitely in the 
input file, they inherit automatically the propagator of the parent multisystem. Therefore, they should *not* be specified individually.
The multisystem container itself always uses an 'empty' propagator, which only implements the operations {{< code OP_UPDATE_INTERACTIONS >}} and {{< code OP_FINISHED >}}.
For {{< code multisystem_t >}} objects, the {{< code "do_td_operations()" >}} routine only implements the {{< code SKIP >}} operation. 
{{% /notice %}}

This routine is triggering the update of the interactions, via the call to {{< code "system%update_interations()" >}},
which loops over the interactions associated with the system, and all required exposed quantities.
### Updating the interactions

Each of the systems (i.e. the multisystem and its subsystems) are attempting to update their interactions with {{< code "system%update_interactions()" >}}.

{{% expand "Implementation of system%update_interations()" %}}
```Fortran
#include_function system_update_interactions
```
{{% /expand %}}

The first part makes sure that the {{< code "update_interactions_start()" >}} routine is only called when no interaction has been updated yet in this time step, iu.e. if their clocks are behind the system clock.

{{< notice note >}}
It is assumed that the interaction can never be ahead of time, compared to the propagator. It is therefore sufficient to ask whether the interaction time equals the propagator time to determine whether an interaction is up-to-date.
{{< /notice >}}

In the second part, we actually attempt to update the interactions, if needed. If an interaction is *not* up-to-date, it first needs to be ensured that the quantities, on which the interaction depends, are up-to-date. Once all quantities are updated, the {{< code "update()" >}} routine of the interaction can be called with the current propagator time as argument. 

Finally, if all interactions are updated successfully, there is a step {{< code "update_interactions_finish()" >}}, which can be necessary for some systems.

{{% expand "Implementation of interaction_with_partner_update()" %}}
```Fortran
#include_function interaction_with_partner_update
```
{{% /expand %}}

This function is to be called with for a specific {{< code requested_time >}}, which is the time at which the interaction is needed.

The function first tries to update the exposed quantities of the partner. If this is successfull, the function calls the {{< code "calculate()" >}} routine of the interaction, sets the interaction clock to the requested time and returns with the result {{< code ".true." >}}, indicating a successful update. In case the function was blocked by one of the exposed quantities, the interaction is not calculated and the function returns {{< code ".false." >}}.


As can be seen, it is, possible that an interaction cannot be updated at a given time. This can be the case when the interaction also depends on quantities of the other interaction partner, and that partner is not yet at the requested time.
 

{{% expand "Definition of system_update_exposed_quantities()" %}}
```Fortran
#include_function system_update_exposed_quantities
```
{{% /expand %}}

This routine demands some more comments:

Firstly, the routine is only allowed if the interaction is of type {{< code "interaction_with_partner_t" >}}, as otherwise there is no partner to update exposed quantities.
It is important to remember who called the function, and which quantities are supposed to be updated:

{{< code "system_update_exposed_quantities()" >}} is called from {{< code "Implementation of interaction_with_partner_update()" >}} for each interaction of a system for the interaction partner of that system with respect to a given interaction.

In the following situations, the routine is **not** allowd to update the exposed quantities:
* we are in an internal SCF loop of the propagation step
* one of the quantities cannot be updated becuase it is either
  - too far behind in time (where two far is more than one clock tick)
  - it is one clock tick behind, but the user requested {{< code "OPTION__INTERACTIONTIMING__TIMING_EXACT" >}}
In these cases, the routine will return a {{< code ".false." >}}. The propagator will continue, but the system in question will be held at this step.
Other systems, however, can continue to progress, and might be able to get to the time, where then this quantity can be updated.

If the timing conditions are fulfilled (and we are not in an internal SCF loop) the system will call the specific routines to update the exposed quantities, and copies their values to the interaction.

### Deadlocks and how to avoid them


