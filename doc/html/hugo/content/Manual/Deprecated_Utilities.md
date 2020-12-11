---
title: "Deprecated Utilities"
series: "Manual"
Weight: 100
---


Some utilities present in older version of octopus have superseded by other utilities or integrated into the main code, in this section you can find how to replace them.

###  octopus_cmplx  

Not exactly a utility, now there is no complex version of the code, the main executable {{< command "octopus" >}} can calculate real or complex wavefunctions depending of the requirements of input file.

###  oct-sf  

This utility was replaced by [[Manual:External utilities:oct-propagation_spectrum | {{< command "oct-propagation_spectrum" >}}]].

###  oct-cross-section  

This utility was replaced by [[Manual:External utilities:oct-propagation_spectrum | {{< command "oct-propagation_spectrum" >}}]].

###  oct-hs-mult and oct-hs-acc  

These utilities were replaced by [[Manual:External utilities:oct-harmonic-spectrum | {{< command "oct-harmonic-spectrum" >}}]].

###  oct-excite  

This utility was used to calculate the excitation spectrum within linear response. It no longer exists because the linear-response formalism calculations are now done by the main code, by making use of the "casida" run mode. See Manual:Calculation Modes:Casida

###  oct-make-st  

This utility was used to replace some of the Kohn-Sham states in the restart directory by Gaussians wave packets. This sort of operation can now be done using the {{< Variable2 "UserDefinedStates" >}} variable.

###  oct-rotatory_strength and oct-rsf  

This utility has been removed. It is replaced by a new {{< Variable2 "PropagationSpectrumType" >}} in {{< file "oct-propagation_spectrum" >}}.

{{< manual_foot prev="Manual:External utilities:oct-xyz-anim" next="Manual:Examples:Hello world" >}}
---------------------------------------------
