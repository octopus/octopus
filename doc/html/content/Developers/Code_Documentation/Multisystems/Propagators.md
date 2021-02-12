---
Title: Propagators
Weight: 2
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

Propagators are implemented as a [state machine](https://en.wikipedia.org/wiki/Finite-state_machine). In every iteration step of the main td-loop ({{< developers "Code_Documentation:Multisystems:Time_propagation" "see Time Propagation" >}}), each the propagator is progressed by one algorithmic step.


### The propagator class

{{% expand "Definition of propagator_t" %}}
```Fortran
#include_type_def propagator_t
```
{{% /expand %}}

The type {{< code propagator_t >}} is an extension of {{< code algorithm_t >}}. Therefore, it contains the list of operations, which define the propagator.
The elements of the propagator algorithm are defined as {{< developers "Code_Documentation:Multisystems:Algorithms#algorithmic-operations" "algorithmic operations" >}}.
The complete propagation algorithm is then defined by adding each step of the algorithm to the propagator (see the Verlet example, below).

{{< notice note >}}
Note, that at this level, progators are independent of the actual implementation of each step. These have to be implemented within the system, for which the propagator will be applied.
{{< /notice >}}



### Example: the Verlet propagator

The Verlet operator is represented by the type {{< code "propagator_verlet_t" >}} which extends {{< code "propagator_t" >}}:


{{% expand "Definition of propagator_verlet_t" %}}
```Fortran
#include_type_def propagator_verlet_t
```
{{% /expand %}}

Fort the Verlet propagator, we define the following operations:
```Fortran
#include_code_doc verlet_propagation_operations
```
These are used to define the algorithm, which is done in the constructor of the propagator:
```Fortran
#include_function propagator_verlet_constructor
```
Note, that the algorithm also uses steps, which are not specific to the Verlet algorithm, and are defined in the {{< code "propagator_oct_m" >}} module.

So far, the propagator is only defined in an abstract way. It is down to the individual systems (or better, their programmers) to implement the individual steps of the algorithm, in terms of the dynamic variables for that system. An easy examample to demonstrate this is the classical particle, as implemented in {{< code "classical_particles_t" >}}

{{% expand "definition of classical_particles_t" %}}
```Fortran
#include_type_def classical_particles_t
```
{{% /expand %}}
This describes an extremely simple system, consisting of a set of classical, neutral particles. The dynamic variables (i.e. the state) of the system are the positions, velocities and accelerations.

One ''tick'' of the propagator is defined in the function {{< code "classical_particle_do_td" >}}. As can be seen in the code below, this function implements all possible algorithmic steps for all propagators, allowed for that system.


{{% expand "Example implementation for classical particles" %}}
```Fortran
#include_subroutine classical_particles_do_td
```
{{% /expand %}}

