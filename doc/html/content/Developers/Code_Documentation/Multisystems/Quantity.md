---
Title: Exposed quantities
Weight: 3
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


Systems can expose quantities that can be used to calculate interactions
with other systems.

Some quantities are dynamical variables of the system. Such quantities are
usually updated by the propagation algorithm and cannot be calculated
on-demand. Such quantities must be marked as "protected".

The module {{< source "multisystem/quantity.F90" "quantity.F90" >}} defines the parameters, which act as index to an exposed quantity within a system.
```Fortran
#include_code_doc quantity
```

Any system, through its base class {{< code "interaction_partner_t" >}} owns an array of type {{< code "quantity_t" >}}

```Fortran
#include_type_def quantity_t
```

This determines whether a quanity is required for a given system, and also associates a specific clock with each quantity.

{{< notice info >}}
''Protected'' quantities deserve some extra mention: Protected quantities are updated in the propagation step of the system, and do not need to be (cannot be) updated explicitely by {{< code "update_exposed_quantity()" >}}. 
{{< /notice >}}