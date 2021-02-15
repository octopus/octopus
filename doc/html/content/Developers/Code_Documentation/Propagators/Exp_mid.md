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
