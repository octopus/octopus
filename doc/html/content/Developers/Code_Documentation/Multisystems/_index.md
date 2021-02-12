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

### Multisystem mode

In the new multisystem mode, from the user perspective, the highest level system is a ''multisystem''. Multisystems are containers which can host other system types, including other multisystems. From the code perspective, the {{< code multisystem_t >}} type is a special case of the {{< code system_t >}} type (i.e. it {{< emph extends >}}  {{< code system_t >}}).

