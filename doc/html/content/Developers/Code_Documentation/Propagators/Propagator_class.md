---
Title: Propagator class
Weight: 1
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


### The propagator class

{{% expand "Definition of propagator_t" %}}
```Fortran
#include_type_def propagator_t
```
{{% /expand %}}

The type {{< code propagator_t >}} is an extension of {{< code algorithm_t >}}. Therefore, it contains the list of operations, which define the propagator.
The elements of the propagator algorithm are defined as {{< developers "Code_Documentation:Multisystems:Algorithms#algorithmic-operations" "algorithmic operations" >}}.
The complete propagation algorithm is then defined by adding each step of the algorithm to the propagator (see the examples).

{{< notice note >}}
Note, that at this level, progators are independent of the actual implementation of each step. These have to be implemented within the system, for which the propagator will be applied.
{{< /notice >}}

Here, we define the following operations:
```Fortran
#include_code_doc general_propagation_operations
```
These operations are general and not bound to a specific propagator, or a specific system. Therefore, they are implemented in the {{< code "system_t" >}} class. For a discussion, see the section on {{< developers "Code_Documentation:Multisystems:Time_propagation" "time propagation" >}}.


The class procedures of {{< code "propagator_t" >}} are those. handling the internal state of the propagator.

Specific propagators are defined as classes extending {{< code "propagator_t" >}}. The necessary specific algorithmic steps are to be defined in the scope of the module file, containing the extending class.

Examples are:
* {{< developers "Code_Documentation:Propagators:Verlet" "Verlet algorithm" >}}
* {{< developers "Code_Documentation:Propagators:Beeman" "Beeman algorithm" >}}
* {{< developers "Code_Documentation:Propagators:Exp_mid" "exponential midpoint algorithm" >}}