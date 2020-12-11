---
title: "Units"
series: "Manual"
---


Before entering into the physics in {{< octopus >}} we have to address a very important issue: {{< emph "units" >}}. There are different unit systems that can be used at the atomic scale: the most used are atomic units and what we call "convenient" units. Here we present both unit systems and explain how to use them in {{< octopus >}}.

###  Atomic Units  

Atomic units are a Gaussian system of units (by "Gaussian" it means that the
vacuum dielectric constant has no dimensions and is set to be
$\epsilon_0 = {1 \over {4\pi}}$),
in which the numerical values of the Bohr radius, the electronic
charge, the electronic mass, and the reduced Planck's constant are set to one:

$$
(1)\\qquad a\_0 = 1; e^2 = 1; m\_e = 1; \\hbar = 1.
$$

This simplifies formulae (although some may feel it presents a serious hazard for dimensional analysis,
interpretation and understanding of formulae, and physics in general.
But this is just a personal taste).
This sets directly two fundamental units, the atomic units of length and of mass:

$$
(2)\\qquad {\\rm au}\_{\\rm length} = a\_0 = 5.2917721\\times 10^{-11}~{\\rm m};\\quad 
{\\rm au}\_{\\rm mass} = m\_e = 9.1093819\\times 10^{-31}~{\\rm kg}.
$$

Since the squared charge must have units of energy times length, we can thus
set the atomic unit of energy

$$
(3)\\qquad {\\rm au}\_{\\rm energy} = {e^2 \\over a\_0} = 4.3597438\\times 10^{-18}~{\\rm J},
$$

which is called Hartree, {{< hartree >}}. And, since the energy has units of mass times
length squared per time squared, this helps us get the atomic unit of time:

$$
(4)\\qquad {\\rm Ha} = m\_e { a\_0^2 \\over {\\rm} {\\rm au}\_{\\rm time}^2} \\to 
{\\rm au}\_{\\rm time} = a\_0 \\sqrt{m\_e \\over {\\rm Ha}} = {a\_0 \\over e} \\sqrt{m\_e a\_0}
= 2.4188843\\times 10^{-17}~{\\rm s}.
$$

Now the catch is: what about Planck's constant? Its dimensions are of energy
times time, and thus we should be able to derive its value by now. But at the
beginning we set it to one! The point is that the four physics constants
used ($a_0, m_e, e^2, \hbar$) are not independent, since:

$$
(5)\\qquad a\_0 =  { \\hbar^2 \\over {m\_e \\; {e^2 \\over {4 \\pi \\epsilon\_0} } } }.
$$

In this way, we could actually have derived the atomic unit of time in an
easier way, using Planck's constant:

$$
(6)\\qquad \\hbar = 1\\; {\\rm Ha}\\,{\\rm au}\_{\\rm time} \\Rightarrow {\\rm au}\_{\\rm time} = { \\hbar \\over {\\rm Ha}} = 
{ {\\hbar a\_0} \\over e^2}\\,.
$$

And combining (6) and (5) we retrieve (4).

###  Convenient Units  

Much of the literature in this field is written using Ångströms and electronvolts
as the units of length and of energy, respectively. So it may be "convenient"
to define a system of units, derived from the atomic system of units, in which
we make that substitution. And so we will call it "convenient".

The unit mass remains the same, and thus the unit of time must change, being
now $\hbar /{\rm eV}\,$,
with $\hbar = 6.582\,1220(20)\times 10^{-16}~\rm eV\,s$.

###  Units in {{< octopus >}}  

Except where otherwise noted, {{< octopus >}} expects all values in the input file to be in atomic units. If you prefer to use other units in the input file, the code provides some handy conversion factors. For example, to write some length value in Ångströms, you can simply multiply the value by {{< code "angstrom" >}}:
```text
 Spacing = 0.5*angstrom
```
A complete list of units {{< octopus >}} knows about can be found in the {{< Variable2 "Units" >}} variable description. 

By default {{< octopus >}} writes all values atomic units. You can switch to convenient units by setting the variable {{< Variable2 "UnitsOutput" >}} to {{< value "ev_angstrom" >}}.

###  Mass Units  

An exception for units in {{< octopus >}} is mass units. When dealing with the mass of ions, [http://en.wikipedia.org/wiki/Amu atomic mass units (amu)] are always used. This unit is defined as $1/12$ of the mass of the <sup>12</sup>{{< name "C" >}} atom. In keeping with standard conventions in solid-state physics, effective masses of electrons are always reported in units of the electron mass (''i.e.'' the atomic unit of mass), even in the eV-Å system.

###  Charge Units  

In both unit systems, the charge unit is the electron charge ''e'' (''i.e.'' the atomic unit of charge).

###  Unit Conversions  

Converting units can be a very time-consuming and error-prone task when done by hand, especially when there are implicit constants set to one, as in the case of atomic units. That is why it's better to use as specialized software like [http://www.gnu.org/software/units/units.html Gnu Units].

In some fields, a very common unit to express the absorption spectrum is Mb. To convert a strength function from 1/eV to Mb, multiply by $\pi h c r_e\,$, with $r_e=e^2/(m_e c^2)\,$. The numerical factor is 109.7609735.


{{< manual_foot prev="Manual:Running Octopus" next="Manual:Physical System" >}}

---------------------------------------------
