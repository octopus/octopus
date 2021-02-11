---
Title: "System clocks"
Weigth: 5
---

Clocks are essential to keep things synchronized. This is also the case in {{<octopus>}}, where instances of the {{< code clock_t >}} type are used to ensure that different systems stay in sync during time propagation.

The smallest unit of time in {{< octopus >}} is one {{< code tick >}}.


{{% expand "Definition of clock_t" %}}
```Fortran
#include_type_def clock_t
```
{{% /expand %}}
