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

The propagators for the system(s) are initialized in {{< code "time_dependent_run_multisystem()" >}} by a call to {{< code "systems%init_propagators()" >}}.
As {{< code "time_dependent_run_multisystem()" >}} is called with {{< code "systems" >}} being of type {{< code "multisystem_basic_t" >}}, this translates to a call to
{{< code "multisystem_init_propagator()" >}}.
{{% expand "Definition of multisystem_init_propagator()" %}}
```Fortran
#include_subroutine multisystem_init_propagator
```
{{% /expand %}}

This routine creates a minimal propagator for the multisystem itself, which is defined by {{< code "propagator_constructor()" >}}:
{{% expand "Definition of propagator_constructor()" %}}
```Fortran
#include_function propagator_constructor
```
{{% /expand %}}

Then, {{< code "multisystem_init_propagator()" >}} loops over the subsystems, and calls {{< code "system%propagator_init()" >}} for each.

In general (i.e. if the function is not overloaded by a derived class), systems use the function {{< code "system_propagator_init()" >}}: 
{{% expand "Definition of system_init_propagator()" %}}
```Fortran
#include_subroutine system_init_propagator
```
{{% /expand %}}

{{< code "system_propagator_init()" >}} parses the input file, looking for the variable {{< variable "TDSystemPropagator" >}}. Here it is important to remember
the order in which the parser treats namespaces.

If the system is specified in the input file as:
```text
%Systems
"System_A" | <type_A>
"System_B" | <type B>
%
```
the top-level multisystem will internally be called {{< code "." >}} and we have in total three namespaces:
```text
. 
./System_A
./System_B
```

The {{< code "system%propagator_init()" >}} routines for {{< code "System_A" >}} and {{< code "System_B" >}} will first look whether {{< variable "TDSystemPropagator" >}} is defined in their respective namespace. If this is _not_ the case, they will look for the variable in the parent namespace, which here is the global namespace {{< code "." >}}.

{{% notice warning %}}
The code does allow to specify different propagators for the multisystem, and the subsystems. While this _might_ work if the subsystems do _not_ interact, it will most likely fail for interacting systems. Therefore, it is highly recommended **not** to specify the propagators for the subsystems separately, unless one knows exactly what one is doing.
{{% /notice %}}


