---
title: "Casida"
series: "Manual"
weight: 4
description: "Alternative to linear response"
---


Mark Casida's formulation of Linear-Response TDDFT allows calculations of the excitation energies of a finite system.  For small molecules, this is normally the fastest way to calculate them.

To perform a Casida calculation you first need a {{< Manual "Ground State" "ground-state" >}} calculation and a calculation of unoccupied states; then you run the code with {{< Variable2 "CalculationMode" >}}{{< code "=casida" >}}.

This is an input file for a linear-response calculation of the nitrogen dimer, found at {{< file "SHARE/testsuite/linear_response/01-casida.*" >}}.

```text
                                                                                
%CalculationMode
gs | unocc | casida
"ground_state_" |  "unocc_states_" | "lrtddft_"
1 | 2 | 3
%
                                                                                
FromScratch = yes
                                                                                
bond_length = 2.0744
                                                                                
%Coordinates
"N" |  -bond_length/2 |  0.0 |  0.0 | no
"N" |   bond_length/2 |  0.0 |  0.0 | no
%
                                                                                
%Species
"N" | 14.0067000 | tm2 | 7 | 2 | 0
%
                                                                                
BoxShape = sphere
                                                                                
Radius = 12.0
Spacing = 0.36
                                                                                
SpinComponents = unpolarized
                                                                                
XFunctional = lda_x
CFunctional = lda_c_vwn
                                                                                
MaximumIter = 200
ConvRelDens = 1e-5
                                                                                
LCAOStart = lcao_full
LCAODimension = 18
                                                                                  
EigenSolver = cg_new
EigenSolverInitTolerance = 1e-2
EigenSolverFinalTolerance = 1e-5
EigenSolverFinalToleranceIteration = 6
EigenSolverMaxIter = 20
unocc_states_EigenSolverMaxIter = 1000
                                                                                
TypeOfMixing = broyden
                                                                                
NumberUnoccStates = 9
                                                                                
PoissonSolver = fft_corrected
```
</pre>


{{< manual_foot prev="Manual:Time-Dependent" next="Manual:Linear_Response" >}}
---------------------------------------------
