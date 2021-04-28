---
Title: "Verlet Propagation"
tutorials: "Multisystem"
weight: 4
---

We are now ready to perform the time propagation of the system. 
The following input file 

{{% expand "input file" %}}
{{< code-block >}}
#include_input testsuite/multisystem/04-propagation_verlet.02-three_body.inp
{{< /code-block >}}
{{% /expand %}}

The input file specifies the propagator as:
{{< code-block >}}
{{< variable "TDSystemPropagator" >}} = verlet
{{< /code-block >}}

According to [Wikipedia](https://en.wikipedia.org/wiki/Verlet_integration), the Verlet algorithm is defined as:

- Calculate $\vec{x}(t + \Delta t) = \vec{x}(t) + \vec{v}(t) \Delta t + \tfrac12 \vec{a}(t) \Delta t^2$.
- Calculate interaction potential using $\vec{x}(t + \Delta t)$.
- Derive $\vec{a}(t + \Delta t)$ from the interaction potential using $\vec{x}(t + \Delta t)$.
- Calculate $\vec{v}(t + \Delta t) = \vec{v}(t) + \tfrac12 \big(\vec{a}(t) + \vec{a}(t + \Delta t)\big)\Delta t$.

These 4 steps have to be performed for each time step.

{{% expand "Detailled description of the algorithmic steps" %}}
The following graph depicts all algorithmic steps of the state machine (see the {{% developers "Code_Documentation/Propagators/" "description of the propagator implementation" %}}).
{{< d3-sequence file="/develop/graph_data/propagation-3body-verlet-equal-step.json" >}}

Each box in this diagram represents a function call of Octopus, which is related to the algorithmic steps of the propagation.
It can be seen how all three systems step through the above algorithmic steps of the Verlet algorithm. You can use the (+) and (-) buttons to manually step through the algorithm to see in which order the steps are performed, and in which steps the various clocks of the systems are advanced.
The interactive graph also allows to display the (unphysical) container systems and the so-called ''ghost'' interactions, which are always present and are an internal tool to ensure that systems stay synchronized, even if they do not have a physical interaction between them.

You can set the variable {{% code-inline %}}{{< variable "Debug" >}} = info{{% /code-inline %}} in the input file and (re-)run the example. Look for the ''Debug:'' lines in the output can compare to the graph above.
{{% /expand %}}

