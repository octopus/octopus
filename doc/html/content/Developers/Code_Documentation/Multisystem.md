---
Title: "Multisystem support"
Weight: 1
---

### Introduction

Support for multisystem in {{< octopus >}} is implemented through an object-oriented framework.

Currently, two major modes are implemented: 

1. The legacy mode, which only allows for one ''matter'' system, consisting of electrons, ions and external fields.
2. The new multisystem framework, which allows for several coupled systems, e.g. maxwell, charged particles, etc. 

At the time of writing (Feb. 2021), electrons and ions are not yet available as separate systems.

### Legacy mode

### Multisystem mode

In the new multisystem mode, from the user perspective, the highest level system is a ''multisystem''. Multisystems are containers which can host other system types, including other multisystems. From the code perspective, the {{< code multisystem_t >}} type is a special case of the {{< code system_t >}} type (i.e. it {{< emph extends >}}  {{< code system_t >}}).

#### {{< code system_t >}}

The {{< code system_t >}} type is the abstract type for all systems. 

{{% expand "Definition of system_t" %}}
```Fortran
#include_type_def system_t
```
{{% /expand %}}

#### {{< code multisystem_t >}}

{{% expand "Definition of multisystem_t" %}}
```Fortran
#include_type_def multisystem_t
```
{{% /expand %}}


<!--
{{< mermaid >}}
 classDiagram
    class interaction_partner_t{
    + namespace_t: namespace
    + clock_t    : clock
    + integer_list_t: supported_interactions_as_partner
    + quantity_t   : quantities[MAX_QUANTITIES]

    + update_exposed_quantities()
    + update_exposed_quantity()
    + copy_quantities_to_interaction() 


    }
    
    class system_t{
    - integer : accumulated_loop_ticks
    + space_t: space
    + propagator_t: prop 
    + integer : interaction_timing  
    + integer_list_t: supported_interactions
    + interaction_list_t: interactions 
    + mpi_grp_t : grp  

    + dt_operation()
    + reset_clocks()
    + update_exposed_quantities()
    + init_propagator()
    + init_all_interactions()
    + init_parallelization()
    + update_interactions()
    + update_interactions_start()
    + update_interactions_finish()
    + propagation_start()
    + propagation_finish()
    + has_reached_final_propagation_time
    + output_start()
    + output_write()
    + output_finish()
    + process_is_slave()
    + init_interaction()
    + initial_conditions()
    + do_td_operation()
    + iteration_info()
    + is_tolerance_reached()
    + update_quantity()

    }
    
    class multisystem_t {
    }

    class interaction_t{
    - integer : n_system_quantities
    - integer : system_quantities(:) 
    + clock_t: clock 
    + character(len=:) : label
    }

    class interaction_with_partner_t{
    + interaction_partner_t: partner
    + update(requested_time)
    }

    interaction_partner_t <|-- system_t
    system_t <|-- multisystem_t

    interaction_t <|-- interaction_with_partner_t
{{< /mermaid >}}
-->

