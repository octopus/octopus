---
Title: "Propagators"
Weight: 12
Description: "How to contribute to the code."
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


Propagators are implemented as a [state machine](https://en.wikipedia.org/wiki/Finite-state_machine). In every iteration step of the main td-loop ({{< developers "Code_Documentation:Multisystems:Time_propagation" "see Time Propagation" >}}), each the propagator is progressed by one algorithmic step.

* {{< developers "Code_Documentation:Propagators:Propagator_class" "Description of the propagator class " >}} 
* {{< developers "Code_Documentation:Propagators:Algorithms" "Description of the algorithm class " >}} 
* Examples are:
  * {{< developers "Code_Documentation:Propagators:Verlet" "Verlet algorithm" >}} (details on how to implement propagators)
  * {{< developers "Code_Documentation:Propagators:Beeman" "Beeman algorithm" >}} (example of predictor-corrector propagators)
  * {{< developers "Code_Documentation:Propagators:Exp_mid" "exponential midpoint algorithm" >}}