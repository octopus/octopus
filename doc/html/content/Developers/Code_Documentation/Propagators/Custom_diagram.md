---
Title: "Custom Diagram"
Weight: 20
---

This page allows custom propagation diagrams, for arbitrary time propagation runs.

To generate the required input files, currently a special branch of Octopus needs to be compiled and run. Later, the functionality will be included in main branches.
You need to checkout the ''multisystem_debug'' branch from the {{< octopus-git >}} repository.

In order to enable the output, you need to set
{{< code-block >}}
{{< variable "Debug">}} = info
{{< /code-block >}}

{{< notice note >}}
Please, keep in mind that a lot of output is produced, and lower the number of time steps to the minimum to complete a few steps for all subsystems
{{< /notice >}}

Use the button below to navigate to the {{< file "propagation.txt" >}} of your calculation.

{{< d3-sequence file="load"  viewContainers="yes" viewGhosts="yes" >}}
