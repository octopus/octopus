---
Title: "Exponential Midpoint"
Weight: 12
---

For the exponential midpoint propagator, we need to define the following operations:
```Fortran
#include_code_doc exp_mid_propagation_operations
```
These are used to define the algorithm, which is done in the constructor of the propagator:
```Fortran
#include_function propagator_exp_mid_constructor
```
### The timeline explained



{{< d3-sequence file="/develop/graph_data/propagation-exp-midpoint-3.json" viewContainers="yes" viewGhosts="yes" >}}

This graph illustrates how the state machine is stepping through the algorithm. Each system is picking the next algorithmic step from the propagator. For the containers (i.e. ''root'' and ''earth''), the only steps are ''Updating interactions'' and ''Finished''. The {{< emph real >}} systems, on the other hand, are progressing orderly through the operations, defined in the propagator.
