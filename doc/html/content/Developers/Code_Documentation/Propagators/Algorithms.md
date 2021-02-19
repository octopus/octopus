---
Title: "Algorithms"
Weight: 3
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


In {{< octopus >}}, many operations, such as time propagations, geometry optimization, etc. are implemented in terms of algorithms.
An algorithm, in general, contains a set of instructions, which are performed in a well-defined order. 


### Algorithm container

The {{< code algorithm_t >}} class itself, is only an extention of the {{< developers "Miscellanea:Linked_list" "linked list" >}}.

{{% expand "Definition of algorithm_t" %}}
```Fortran
#include_type_def algorithm_t
```
{{% /expand %}}

[Algorithmic operations](#algorithmic-operations) can be added with the function {{< code "add_operation()" >}}. Examples are discussed in the section 
{{< developers "Code_Documentation:Propagators" "Propagators" >}}.

### Algorithm iterator

{{% expand "Definition of algorithm_iterator_t" %}}
```Fortran
#include_type_def algorithm_iterator_t
```
{{% /expand %}}

### Algorithmic operations 

{{% expand "Definition of algorithmic_operation_t" %}}
```Fortran
#include_type_def algorithmic_operation_t
```
{{% /expand %}}

Some global algorithmic steps are defined in {{< source "multisystem/propagator.F90" >}}:
```Fortran
#include_code_doc general_propagation_operations
```

Derived propagators can then add their own steps, such as e.g. in {{< source "multisystem/propagator_verlet.F90" >}}:
```Fortran
#include_code_doc verlet_propagation_operations
```

{{% notice note %}}
It is important to stress here, that one algorithmic step, in general, does not advance any clock. Clocks are only advanced in steps which update the corresponding entity, which could be a system, an exposed quantity or an interaction. It is therefore quite common, that at a specific state of the algorithm, clocks or different entities have different values.
{{% /expand %}}