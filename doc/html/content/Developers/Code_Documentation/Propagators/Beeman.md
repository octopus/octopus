---
Title: "Beeman"
Weight: 11
---

[Beeman's algorithm](https://en.wikipedia.org/wiki/Beeman%27s_algorithm) is one of the simplest which allows a predictor-corrector implementation, and demonstrates an algorithm, including a self-consistent loop within one propagator time-step. The self-contistent predictor-corrector feature is implemented as optional.


For the Beeman propagator, we need to define the following operations:
```Fortran
#include_code_doc beeman_propagation_operations
```
These are used to define the algorithm, which is done in the constructor of the propagator:
```Fortran
#include_function propagator_beeman_constructor
```

### The timeline explained



{{< d3-sequence file="/develop/graph_data/propagation-beeman-3.json" viewContainers="yes" viewGhosts="yes"  >}}

This graph illustrates how the state machine is stepping through the algorithm. Each system is picking the next algorithmic step from the propagator. For the containers (i.e. ''root'' and ''earth''), the only steps are ''Updating interactions'' and ''Finished''. The {{< emph real >}} systems, on the other hand, are progressing orderly through the operations, defined in the propagator.
