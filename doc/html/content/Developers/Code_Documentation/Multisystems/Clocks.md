---
Title: "Clocks"
Weigth: 4
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


Clocks are essential to keep things synchronized. This is also the case in {{<octopus>}}, where instances of the {{< code clock_t >}} type are used to ensure that different systems stay in sync during time propagation. Clocks are used to determine whether quantities or interactions need to be updated in a propagation step.

The smallest unit of time in {{< octopus >}} is one {{< code tick >}}.


{{% expand "Definition of clock_t" %}}
```Fortran
#include_type_def clock_t
```
{{% /expand %}}
