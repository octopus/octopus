---
Title: "Multisystems"
Weight: 10
Description: "How to contribute to the code."
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


Support for multisystem in {{< octopus >}} is implemented through an object-oriented framework.

Currently, two major modes are implemented: 

1. The legacy mode, which only allows for one ''matter'' system, consisting of electrons, ions and external fields.
2. The new multisystem framework, which allows for several coupled systems, e.g. maxwell, charged particles, etc. 

At the time of writing (Feb. 2021), electrons and ions are not yet available as separate systems.

### Legacy mode

The ''legacy'' mode of the code is used whenever the input file does _not_ have a {{< variable "Systems">}} block.
In this case, the top level system is initialized to be of (the old) {{< code "electrons_t" >}} type. This type describes the combined electron-ion system.
It is planned for te future, that this will be split into the new {{< code electrons_t >}} and {{< code ions_t >}}, which will descibe the electrons and ions as separate systems.
The current {{< code electrons_t >}} is to be replaced by the {{< code matter_t >}}, which is a multisystem, containing electrons and ions.
### Multisystem mode

If the input file contains the {{< variable "Systems" >}} block, the code uses the new multisystems mode.
In this multisystem mode, from the user perspective, the highest level system is a ''multisystem''. Multisystems are containers which can host other system types, including other multisystems. From the code perspective, the {{< code multisystem_t >}} type is a special case of the {{< code system_t >}} type (i.e. it {{< emph extends >}}  {{< code system_t >}}).

The following chapters will discuss in more detail:

* [the multisystem classes](system_classes/)
* [interactions](interactions/)
* [exposed quantities](quantities/)
* [clocks](clocks/)
* [initialization](initialization/)
* [time propagation](time_propagation/)


