---
Title: "Multisystem introduction"
weight: 1
---

## Octopus with multisystem support

{{< octopus >}} now supports the calculation of compound systems. While the aim of this development is to allow coupling matter to Maxwell fields, or coupling a full DFT system to a system, treated by a tight-binding Hamiltonian, in this tutorial, we will start with a very different model, namely celestial bodies. Due to their simplicity, they are well suited to explain how to set up a system, containing various subsystems.



New variables and concepts in this tutorial: 

* {{< variable "Systems" >}}:
* {{< variable "Interactions" >}}
* Namespaces

### Systems

{{< variable "Systems" >}} describes the nature of a subsystem. Examples are

* classical particles
* charged particles
* a Maxwell system
* a tight binding system
* electrons (not yet implemnented in the multisystem framework)

{{< code-block >}}
%{{< variable "Systems">}}
"Name-1" | {{< emph "system-type-1" >}}
"Name-2" | {{< emph "system-type-2" >}}
...
%
{{< /code-block >}}
The possible system types can be found in the reference of {{< variable "Systems" >}}.
The names, given to a system, define namespaces, which will be used in the remainder of the input file, in order to define system specific values for variables.
This can be seen in the example below.

If a (sub-)system is a multisystem, the components need to be specified in another {{< variable "Systems" >}} block.


### Interactions

Different subsystems are pretty boring, if they cannot interact with each other. These mutual interactions are defined in the {{< variable "Interactions" >}} block.

{{< code-block >}}
%{{< variable "Interactions" >}}
{{< emph "interaction_type" >}} | {{< emph "interaction_mode" >}} [ | {{< emph "partners" >}} ]
...
%
{{< /code-block >}}

### Namespaces

Many possible systems share the same variables. In the following examples of celestial bodies, planets are all classical particles (albeit quite big ones), which are specified by their mass. Obviously, they need to be able to have different masses. For that reason, the system-specific variables are defined by their variable name, with the system name added as prefix (or namespace).

In our example, that would be:

{{< code-block >}}
%{{< variable "Systems">}}
"Earth" | classical_particle
...
%

Earth.{{< variable "ParticleMass" >}} = 5.97237e24
%Earth.{{< variable "ParticleInitialPosition" >}}
 -147364661998.16476 | -24608859261.610123 | 1665165.2801353487
%
%Earth.{{< variable "ParticleInitialVelocity" >}}
 4431.136612956525 | -29497.611635546345 | 0.343475566161544
%
{{< /code-block >}}
