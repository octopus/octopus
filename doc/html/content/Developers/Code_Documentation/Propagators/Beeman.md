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
