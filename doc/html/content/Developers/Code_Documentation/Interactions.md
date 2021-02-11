---
Title: Interactions
section: Developers
Weight: 3
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

### Introduction

In physics, interactions are in general symmetric. However, the magnitude with which it affects the interaction partners can be very different.
Also, the way the interactions have to be implemented can differ depending on the ''direction'' of the interaction.

Therefore, in {{< octopus >}} interactions are broken down into their forward and their backward parts, which are treated seperately.

{{% mermaid %}}
flowchart LR
    A[partner A] --forward--> B
    B[partner B] --backward--> A
{{% /mermaid %}}

As we will see, in some approximations, only one of them (usually the forward interaction) are used.

The ''partners'' are instances of {{< code interaction_partner_t >}} derived types. The interactions (from now on uni-directional) are implemented as derived types of the abstract {{< code interaction_t >}} class.

{{% mermaid %}}
flowchart LR
    A[partner A] --interaction--> B[partner B]
{{% /mermaid %}}

An interaction ''belongs'' to partner A and acts on partner B.

### Abstract classes

#### {{< code interaction_t >}}

{{% expand "Definition of interaction_t" %}}
```Fortran
#include_type_def interaction_t
```
{{% /expand %}}

This is the minimal data structure for interactions, which only contains information relating to partner A, who owns the interacion.
In particular, the interaction has a list of quantities, which partner A needs to expose, so that the interaction can be calculated.
Furthermore, this abstract type already contains the clock for the interaction, and defines the interfaces for the {{< code "update()" >}} 
and {{< code "calculate()" >}} routines.

The purpose of the {{< code "update()" >}} routine is to perform the checks *whether* any quantities need to be updated, and is implemented on a high level (in {{< code interaction_with_partner_t >}}) while {{< code "calculate()" >}} is specific to the interactions and will perform the update of the quantities.

Currently, there is no interaction type directly derived from {{< code interaction_t >}}. This type was mainly introduced for clarity, and also in order to break the problem of a circular dependency, which would have occurred otherwise.

#### {{< code interaction_with_partner_t >}}

Curretly, all interactions are derived from {{< code interaction_with_partner_t >}}, which extends {{< code interaction_t >}} and adds the information about partner B, fron here on only called ''partner''.

{{% expand "Definition of interaction_with_partner_t" %}}
```Fortran
#include_type_def interaction_with_partner_t
```
{{% /expand %}}

This adds a list of IDs of the quantities, the partner exposes to the interaction.

Furthermore, at this level, we can implement the {{< code "update()" >}} procedure.

#### {{< code force_interaction_t >}}

The next level of specialization of interaction, are all interactions which create a force on the partner. Here we only add the actual force vector, acting on the partner system.
{{% expand "Definition of force_interaction_t" %}}
```Fortran
#include_type_def force_interaction_t
```
{{% /expand %}}

</br>

### Specific classes:

#### Ghost interaction

{{% expand "Definition of ghost_interaction_t" %}}
```Fortran
#include_type_def ghost_interaction_t
```
{{% /expand %}}

#### Gravity

{{% expand "Definition of gravity_t" %}}
```Fortran
#include_type_def gravity_t
```
{{% /expand %}}

#### Coulomb force

{{% expand "Definition of coulomb_force_t" %}}
```Fortran
#include_type_def coulomb_force_t
```
{{% /expand %}}

#### Lorentz force

{{% expand "Definition of lorentz_force_t" %}}
```Fortran
#include_type_def lorentz_force_t
```
{{% /expand %}}
