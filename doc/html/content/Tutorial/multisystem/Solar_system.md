---
Title: "Example: Celestial dynamics"
tutorials: "Multisystem"
weight: 2
---


In the simplest example we define a system, which consists of similar subsystems. The input file below shows a model describing the motion of the earth, the moon and the sun:

Here, all three subsystems are treated equally. We have two levels of systems: 

* the solar system, which can be tought of as a container
* the celestial bodies (sub, earth and moon)

{{< expand "Expand for input file with two levels of nesting" >}}
{{< code-block >}}
#include_input testsuite/multisystem/01-nested_systems.01-two_levels.inp
{{< /code-block >}}
{{< /expand >}}

While this is a perfectly fine definition of the system, sometimes it is beneficial to group systems together.
For instance, we might define a "Earth" as the combined system of "Terra" and "Luna", which is orbiting the sun:

{{< expand "Expand for input file with three levels of nesting" >}}
{{< code-block >}}
#include_input testsuite/multisystem/01-nested_systems.02-three_levels.inp
{{< /code-block >}}
{{< /expand >}}

