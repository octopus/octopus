---
title: "About Octopus"
series: "Manual"
weight: 1
description: " "
---


{{< Octopus >}} is a software package for {{< name "density-functional theory" >}} (DFT), and {{< name "time-dependent density functional theory" >}} (TDDFT).

##  Developers  

The main development team of this program  is composed of:

* Joseba Alberdi
* Xavien Andrade
* Florian Buchholz
* Alberto Castro
* Tilman Dannert
* Umberto De Giovannini
* Alain Delgado
* Nicole Helbig
* Hannes Huebener
* Joaquim Jornet-Somoza
* Ask Larsen
* Irina Lebedeva
* Miguel A. L. Marques
* Fernando Nogueira
* Micael Oliveira
* Carlo Andrea Rozzi
* Angel Rubio
* Ravindra Shinde
* Jose R. F. Sousa
* David Strubbe
* Iris Theophilou
* Alejandro Varas
* Matthieu Verstraete
* Philipp Wopperer

Former developers:

* Heiko Appel
* Fulvio Berardi
* Johanna Fuks
* David Kammerlander
* Kevin Krieger
* Florian Lorenzen
* Danilo Nitsche
* Roberto Olivares-Amaya
* Arto Sakko
* Axel Thimm
* Jessica Walkenhorst
* Jan Werschnik

Other contributors are:

* Sebastien Hamel: parallel version of oct-excite
* Eugene S. Kadantsev: linear response code

##  Introduction  

{{< octopus >}} is a pseudopotential real-space package aimed at the simulation of the electron-ion dynamics of one-, two-, and three-dimensional ﬁnite systems subject to time-dependent
electromagnetic ﬁelds. The program is based on time-dependent density-functional theory (TDDFT) in the Kohn-Sham scheme. All quantities are expanded in a regular mesh
in real space, and the simulations are performed in real time. The program has been
successfully used to calculate linear and non-linear absorption spectra, harmonic spectra,
laser induced fragmentation, etc. of a variety of systems. The fundamentals of DFT and
TDDFT can be found, e.g., in the books 
[^footnote-1]
and[^footnote-2]
. 
All information about the octopus
package can be found in its homepage, http://www.tddft.org/programs/octopus/, and in the articles 
[^footnote-3]

and 
[^footnote-4]
.

The main advantage of real-space methods is the simplicity and
intuitiveness of the whole procedure. First of all, quantities like the density or
the wave-functions are very simple to visualize in real space. Furthermore,
the method is fairly simple to implement numerically for 1-, 2-, or 3-dimensional
systems, and for a variety of different boundary conditions. For example, one
can study a finite system, a molecule, or a cluster without the need of a super-cell,
simply by imposing that the wave-functions are zero at a surface far enough from the system.
In the same way, an infinite system, a polymer, a surface, or bulk material can be
studied by imposing the appropriate cyclic boundary conditions. Note also that
in the real-space method there is only one convergence parameter, namely the grid-spacing, and that decreasing the grid spacing always improves the result.

Unfortunately, real-space methods suffer from a few drawbacks. For example,
most of the real-space implementations are not variational, i.e., we may
find a total energy lower than the true energy, and if we reduce the grid-spacing
the energy can actually increase. Moreover, the grid breaks translational
symmetry, and can also break other symmetries that the system may possess. This can
lead to the artificial lifting of some degeneracies, to the appearance of
spurious peaks in spectra, etc. Of course all these problems can be minimized
by reducing the grid-spacing.

##  History  

{{< octopus >}} is based on a fixed-nucleus code written by George F. Bertsch and K. Yabana to perform real-time dynamics in clusters 
[^footnote-5]

and on a condensed matter real-space plane-wave based code written by A. Rubio, X. Blase and S.G. Louie 
[^footnote-6]
. 
The code was afterwards extended to handle periodic systems by G.F. Bertsch, J.I. Iwata, A. Rubio, and K. Yabana 
[^footnote-7]
. 
Contemporaneously there was a major rewrite of the original cluster code to handle a vast majority of finite systems. At this point the cluster code was named {{< name "tddft" >}}.

This version was consequently enhanced and beautified by A. Castro (at the time Ph.D. student of A. Rubio), originating a fairly verbose 15,000 lines of {{< name "Fortran 90/77" >}}. In the year 2000, M. Marques (aka Hyllios, aka António de Faria, corsário português), joined the A. Rubio group in Valladolid as a postdoc. Having to use {{< name "tddft" >}} for his work, and being petulant enough to think he could structure the code better than his predecessors, he started a major rewrite of the code together with A. Castro, finishing version 0.2 of {{< name "tddft" >}}. But things were still not perfect: due to their limited experience in {{< name "Fortran 90" >}}, and due to the inadequacy of this language for anything beyond a {{< name "HELLO WORLD" >}} program, several parts of the code were still clumsy. Also the idea of GPLing the almost 20,000 lines arose during an alcoholic evening. So after several weeks of frantic coding and after getting rid of the {{< name "Numerical Recipes" >}} code that still lingered around, {{< octopus >}} was born.

The present released version has been completely rewritten and keeps very little relation to the old version (even input and output files) and has been enhanced with major new flags to perform various excited-state dynamics in finite and extended systems. The code will be updated frequently and new versions can be found here.

If you find the code useful for you research we would appreciate if you give reference to this work and previous ones.

##  Contributing to Octopus  

If you have some free time, and if you feel like taking a joy ride with {{< name "Fortran 90" >}}, just drop us an email. You can also send us patches, comments, ideas, wishes, etc. They will be included in new releases of octopus.

If you found a have a bug, please report it to our Bug Tracking System: http://www.tddft.org/trac/octopus/newticket

##  The {{< octopus >}}  Copying Conditions  

This program is “free”; this means that everyone is free to use it and free to redistribute it on a free basis. What is not allowed is to try to prevent others from further sharing any version of this program that they might get from you.

Specifically, we want to make sure that you have the right to give away copies of the program, that you receive source code or else can get it if you want it, that you can change this program or use pieces of them in new free programs, and that you know you can do these things.

To make sure that everyone has such rights, we have to forbid you to deprive anyone else of these rights. For example, if you distribute copies of the program, you must give the recipients all the rights that you have. You must make sure that they, too, receive or can get the source code. And you must tell them their rights.

Also, for our own protection, we must make certain that everyone finds out that there is no warranty for this program. If these programs are modified by someone else and passed on, we want their recipients to know that what they have is not what we distributed, so that any problems introduced by others will not reflect on our reputation.

The precise conditions of the license are found in the General Public Licenses that accompany it.

Please note that {{< octopus >}} distribution normally comes with some external libraries that are not covered by the GPL license, please see the [[Manual:Appendix:Copying | Copying ]] Appendix for the copying conditions or these packages.

<references />

{{< manual_foot prev="Manual:Octopus" next="Manual:Installation" >}}
---------------------------------------------
[^footnote-1]: {{< Book title="A Primer in Density Functional Theory" author="C. Fiolhais, F. Nogueira, and M.A.L. Marques (editors)" publisher="Springer Berlin Heidelberg New York" isbn="3-540-03082-2" series="Lecture Notes in Physics" year="2006" issn="0075-8450" >}}

[^footnote-2]: {{< Book title="Time-dependent Density Functional Theory" author="M. A. L. Marques and C. A. Ullrich and F. Nogueira and A. Rubio and K. Burke and E. K. U. Gross (editors)" publisher="Springer Berlin Heidelberg New York" isbn="3-540-35422-0" series="Lecture Notes in Physics" year="2006" issn="0075-8450" >}}

[^footnote-3]: {{< Article title="octopus: a first principles tool for excited states electron-ion dynamics" authors="M.A.L. Marques, A. Castro, G. F. Bertsch, and A. Rubio" journal="Comp. Phys. Comm." volume="151" pages="60 " year="2003" >}}

[^footnote-4]: {{< Article title="octopus: a tool for the application of time-dependent density functional theory" authors="A. Castro, H. Appel, M. Oliveira, C. A. Rozzi, X. Andrade, F. Lorenzen, M.A.L. Marques, E. K. U. Gross, and A. Rubio" journal="Phys. Stat. Sol. (b)" volume="243" pages="2465" year="2006" >}}

[^footnote-5]: {{< Article authors="G.F. Bertsch and K. Yabana" title="Time-dependent local-density approximation in real time " journal="Phys. Rev. B" volume="54" pages="4484" year="1996" >}}

[^footnote-6]: {{< Article authors="A. Rubio, X. Blase, and S.G. Louie" title="Ab Initio Photoabsorption Spectra and Structures of Small Semiconductor and Metal Clusters" journal="Phys. Rev. Lett." volume="77" pages="247" year="1996" >}}

[^footnote-7]: {{< Article authors="G.F. Bertsch, J.I. Iwata, A. Rubio, and K. Yabana" title="Real-space, real-time method for the dielectric function" journal="Phys. Rev. B" volume="62" pages="7998" year="2000" >}}

