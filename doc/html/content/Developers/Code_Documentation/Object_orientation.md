---
Title: "Object orientation in Octopus"
Weight: 1
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

### Motivation 

Since version 10, {{< octopus >}} is gradually being converted into a fully object oriented code, using the OOP features of Fortran 2003.

The main features of OOP used in {{< octopus >}} are:

* Encapsulation
* Inheritance

The benefits of OOP design are:

* less code duplication
* better readable code
* code can be closer to the physics
* low level aspects can be hidden


### What is an object?

An object refers to a data-structure, which is bundled with the functions (methods) acting on that data. Usually, it is good practice to declare the data as private, 
and allow access only through so-called access functions (getters and setters). This is called encapsulation, and allows later to change details of the implementation
without affecting the rest of the code.

Objects are instances of classes, which define the data structures and the methods for its objects.





