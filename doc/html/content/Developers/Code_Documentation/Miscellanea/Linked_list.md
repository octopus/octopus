---
Title: "Linked list"
Weight: 1
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


The {{< code "linked_list_t" >}} type implements a linked list of unlimited polymorphic
values. This allows the storage of any type of data. Iterating over the
list is done using the associated iterator. These two classes are not meant
to used as is, but rather to be extended and by providing an add method to the
list and a get_next method to the iterator.

{{% expand "Definition of linked_list_t" %}}
```Fortran
#include_type_def linked_list_t
```
{{% /expand %}}

{{% expand "Definition of list_node_t" %}}
```Fortran
#include_type_def list_node_t
```
{{% /expand %}}

{{% expand "Definition of linked_list_iterator_t" %}}
```Fortran
#include_type_def linked_list_iterator_t
```
{{% /expand %}}
