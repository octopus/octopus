---
title: "Output"
series: "Manual"
weight: 3
description: "Understanding the output"
---


At first you may be quite happy that you have mastered the input file, and {{< octopus >}} runs without errors. However, eventually you (or your thesis advisor) will want to learn something about the system you have managed to describe to {{< octopus >}}.

#### Ground-State DFT 
{{< Octopus >}} sends some relevant information to the standard output (which you may have redirected to a file).  Here you will see energies and occupations of the eigenstates of your system.  These values and other information can also be found in the file {{< file "static/info" >}}.

However {{< octopus >}} also calculates the wavefunctions of these states and the positions of the nuclei in your system.  Thus it can tell you the density of the dipole moment, the charge density, or the matrix elements of the dipole moment operator between different states.  Look at the values that the {{< Variable2 "Output" >}} variable can take to see the possibilities.

For example, if you include  
```text
 <tt>
 Output = wfs_sqmod + potential 
 </tt>
```
in your {{< file "inp" >}} file, {{< Octopus >}} will create separate text files in the directory {{< file "static" >}} with the values of the square modulus of the wave function and the local, classical, Hartree, and exchange/correlation parts of the Kohn-Sham potential at the points in your mesh.  

You can specify the formatting details for these input files with the {{< Variable2 "OutputFormat" >}} variable and the other variables in the {{< Variable2 "Output" >}} section of the Reference Manual.  For example, you can specify that the file will only contain values along the x, y, or z axis, or in the plane x=0, y=0, or z=0.  You can also set the format to be readable by the graphics programs [http://www.opendx.org/ OpenDX], [http://www.gnuplot.info/ gnuplot] or [http://www.mathworks.com/products/matlab/ MatLab].  OpenDX can make plots of iso-surfaces if you have data in three-dimensions.  However gnuplot can only make a 3-d plot of a function of two variables, i.e. if you have the values of a wavefunction in a plane, and 2-d plots of a function of one variable, i.e. the value of the wavefunction along an axis.

#### Time-Dependent DFT 

##### Optical Properties 

A primary reason for using a time-dependent DFT program is to obtain the optical properties of your system.  You have two choices for this, linear-response theory <i>a la</i> [http://scitation.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=JCPSA6000104000013005134000001&idtype=cvips&gifs=yes Jamorski, Casida & Salahub]
[^footnote-1]
, 
or explicit time-propagation of the system after a perturbation, a la [http://prola.aps.org/abstract/PRB/v54/i7/p4484_1  Yabana & Bertsch]
[^footnote-2]
.  
You may wish to read more about these methods in the paper by [http://www3.interscience.wiley.com/cgi-bin/abstract/112650991/ABSTRACT?CRETRY=1&SRETRY=0 Castro et al.]
[^footnote-3]


##### Linear-Response Theory 
Linear-response theory is based on the idea that a small (time-dependent) perturbation in an externally applied electric potential 
$\delta v (r, \omega )$ will result in a (time-dependent) perturbation of the electronic density 
$ \delta \rho (r, \omega )$ which is linearly related to the size of the perturbation:
$\delta \rho (r, \omega ) =  \int d^{3} r' \chi (r, r'; \omega) \delta v (r, \omega )$.  Here, obviously, the time-dependence is Fourier-transformed into a frequency-dependence, $ \omega $. The susceptibility,
$\chi(r, r'; \omega) $, is a density-density response function, because it is the response of the charge density to a potential that couples to the charge density of the system.  Because of this, it has poles at the excitation energies of the many-body system, meaning that the induced density also has these poles.  One can use this analytical property to find a related operator whose eigenvalues are these many-body excitation energies.  The matrix elements of the operator contain among other things:  1) occupied and unoccupied Kohn-Sham states and energies (from a ground state DFT calculation) and 2) an exchange-correlation kernel, 
$ f_{xc}(r, r', \omega) = {{\delta v_{xc}[n(r,\omega)]}\over{\delta n(r',\omega)}}\mid_{\delta v_{ext} = 0}$.

Casida's equations are a full solution to this problem (for real wavefunctions). The Tamm-Dancoff approximation uses only occupied-unoccupied transitions. The Petersilka approximation uses only the diagonal elements of the Tamm-Dancoff matrix, <i>i.e.</i> there is no mixing of transitions.
[^footnote-4]
```text
 It takes only a little more time to calculate the whole matrix, so Petersilka is provided mostly for comparison.
```

These methods are clearly much faster (an order of magnitude) than propagating
in time, but it turns out that they are very sensitive to the quality
of the unoccupied states. This means that it is very hard to converge the
excitation energy, because one requires a very large simulation box (much
larger than when propagating in real time).

##### Electronic Excitations by Means of Time-Propagation 
See {{< Manual "Time_Dependent-Delta_kick:_Calculating_an_absorption_spectrum" "Time_Dependent-Delta_kick:_Calculating_an_absorption_spectrum" >}}.

#### References 
<references/>


{{< manual_foot prev="Manual:Discretization" next="Manual:Troubleshooting" >}}
---------------------------------------------
[^footnote-1]: {{< Article title="Dynamic polarizabilities and excitation spectra from a molecular implementation of time-dependent density-functional response theory: N2 as a case study" authors="Christine Jamorski, Mark E. Casida, and Dennis R. Salahub " journal="J. Chem. Phys." year="1996" vol="104" issue="13" pages="5134-5147" url="http://scitation.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=JCPSA6000104000013005134000001&idtype=cvips&gifs=yes" doi="10.1063/1.471140 " >}}

[^footnote-2]: {{< Article title="Time-dependent local-density approximation in real time " authors="K. Yabana, G. F. Bertsch" journal="Phys. Rev. B" year="1996" vol="54" issue="7" pages="4484 - 4487" url="http://link.aps.org/abstract/PRB/v54/p4484" doi="10.1103/PhysRevB.54.4484" >}}

[^footnote-3]: {{< Article title="octopus: a tool for the application of time-dependent density functional theory" authors="Alberto Castro, Heiko Appel, Micael Oliveira, Carlo A. Rozzi, Xavier Andrade, Florian Lorenzen, M. A. L. Marques, E. K. U. Gross, Angel Rubio " journal="physica status solidi (b)" year="2006" volume="243" issue="11" pages="2465-2488" url="http://www3.interscience.wiley.com/cgi-bin/abstract/112650991/ABSTRACT?CRETRY=1&SRETRY=0" doi="10.1002/pssb.200642067" >}}

[^footnote-4]: {{< Article title="Excitation Energies from Time-Dependent Density-Functional Theory" authors="Petersilka, M.  and Gossmann, U. J. and Gross, E. K. U." journal="Phys. Rev. Lett." volume="76" number="8" pages="1212--1215" numpages="3" year="1996" month="Feb" doi="10.1103/PhysRevLett.76.1212" publisher="American Physical Society" url="http://prola.aps.org/abstract/PRL/v76/i8/p1212_1" >}}

