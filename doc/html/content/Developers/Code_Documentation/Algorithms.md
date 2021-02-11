---
Title: "Algorithms"
Weight: 10
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


In {{< octopus >}}, many operations, such as time propagations, geometry optimization, etc. are implemented in terms of algorithms.
An algorithm, in general, contains a set of instructions, which are performed in a well-defined order. 

The {{< code algorithm_t >}} class itself, is only an extention of the {{< developers "Miscellanea:Linked_list" "linked list" >}}.

{{% expand "Definition of algorithm_t" %}}
```Fortran
#include_type_def algorithm_t
```
{{% /expand %}}

{{% expand "Definition of algorithmic_operation_t" %}}
```Fortran
#include_type_def algorithmic_operation_t
```
{{% /expand %}}

{{% expand "Definition of algorithm_iterator_t" %}}
```Fortran
#include_type_def algorithm_iterator_t
```
{{% /expand %}}
