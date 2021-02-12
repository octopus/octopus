---
Title: "Multisystem classes"
Weight: 5
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}



### Abstract classes


#### {{< code interaction_partner_t >}}

An interaction_partner in {{< octopus >}} is anything, which can have an interaction with any other system.
For instance electrons and ions are interaction partners, but also the photonic system described by the Maxwell system, or external potentials.
Therefore, it is the base class of all possible systems and multisystems.

{{% expand "Definition of interaction_partner_t" %}}
```Fortran
#include_type_def interaction_partner_t
```
{{% /expand %}}

Each {{< name interaction_partner >}} is associated with a {{< name namespace >}}, owns a {{< name clock >}}, as well as a list of {{< developers "code_documentation/multisystems/interactions" "interactions" >}} in which it can be a partner, and a list of physical quantities, which are exposed to other systems (through the interactions) and a list
or quantities, which it exposes to the interactions.

It also provides the basic functions to update exposed quantities, and copy them to the interaction. More details about this mechanism are described in the section on 
{{< developers "code_documentation/multisystems/interactions" "interactions" >}}.

#### {{< code system_t >}}

The {{< code system_t >}} type is the abstract type for all systems. 
As all possible systems are potential partners of some interaction, the {{< code system_t >}} type itself extends the abstract {{< code interaction_partner_t >}} type.

{{% expand "Definition of system_t" %}}
```Fortran
#include_type_def system_t
```
{{% /expand %}}

{{% expand "Definition of system_dt_operation" %}}
```Fortran
#include_subroutine system_dt_operation
```
{{% /expand %}}

#### {{< code multisystem_t >}}


{{% expand "Definition of multisystem_t" %}}
```Fortran
#include_type_def multisystem_t
```
{{% /expand %}}

{{< notice note >}}
{{< code multisystem_t >}} is an abstract class and cannot be used as such to describe a set of systems in the code.
The {{< code type >}} to be used for combined systems is {{< code multisystem_basic_t >}}, which extends {{< code multisystem_t >}}.
{{< /notice >}}




### Specific classes

#### {{< code multisystem_basic_t >}}

{{% expand "Definition of multisystem_basic_t" %}}
```Fortran
#include_type_def multisystem_basic_t
```
{{% /expand %}}

{{< code "multisystem_basic_t" >}} is a specific (i.e. non-abstract) container type, which can host other systems. Its propagator

{{% expand "Definition of classical_particles_t" %}}
```Fortran
#include_type_def classical_particles_t
```
{{% /expand %}}

{{% expand "Definition of electrons_t" %}}
```Fortran
#include_type_def electrons_t
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

