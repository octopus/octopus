---
Title: Propagators
Weight: 5
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
The elements of the propagator algorithm are defined as {{< developers "Code_Documentation:Algorithms#algorithmic-operations" "algorithmic operations" >}}.
The complete propagation algorithm is then defined by adding each step of the algorithm to the propagator (see the Verlet example, below).

{{< notice note >}}
Note, that at this level, progators are independent of the actual implementation of each step. These have to be implemented within the system, for which the propagator will be applied.
{{< /notice >}}



### Example: Verlet

{{% expand "Definition of propagator_verlet_t" %}}
```Fortran
#include_type_def propagator_verlet_t
```
{{% /expand %}}


```Fortran
#include_function propagator_verlet_constructor
```
