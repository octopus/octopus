---
Title: Initialization
Weight: 5
---

### Constructing the systems

In the multisystem framework, it is important to detail how the various components are initialized.
The top level system is created in {{< source "main/run.F90" >}} by calling the {{< code multisystem_basic_t >}} constructor
{{% expand "Definition of multisystem_basic_constructor()" %}}
```Fortran
#include_function multisystem_basic_constructor
```
{{% /expand %}}

which itself calls {{< code "multisystem_basic_init()" >}} and {{< code "multisystem_init()" >}}.

{{% expand "Definition of multisystem_basic_init()" %}}
```Fortran
#include_subroutine multisystem_basic_init
```
```Fortran
#include_subroutine multisystem_init
```
{{% /expand %}}

Finally, {{< code "multisystem_create_system()" >}} loops over the subsystems and calls the system factory to create each of them.
{{% expand "Definition of multisystem_create_system()" %}}
```Fortran
#include_subroutine multisystem_create_system
```
{{% /expand %}}

### Initializing the propagators

