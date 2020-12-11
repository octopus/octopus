---
title: "Input file"
series: "Manual"
weight: 2
description: "About the input file"
---


{{< Octopus >}} uses a single input file from which to read user instructions to know what to calculate and how. This page explains how to generate that file and what is the general format. The Octopus parser is a library found in the {{< file "liboct_parser" >}} directory of the source, based on bison and C. You can find two (old) separate release versions of it at the bottom of the [[Releases]] page.

###  Input file  

Input options should be in a file called {{< file "inp" >}}, in the directory {{< octopus >}} is run from. This is a plain {{< name "ASCII" >}} text file, to create or edit it you can use any text editor like {{< name "emacs" >}}, {{< name "vi" >}}, {{< name "jed" >}}, {{< name "pico" >}}, {{< name "gedit" >}}, etc. For a fairly comprehensive example, just look at the tutorial page [[Tutorial:Nitrogen_atom]].

At the beginning of the program, the parser reads the input file, parses it, and generates a list of variables that will be read by {{< octopus >}} (note that the input is case-independent). There are two kind of variables: scalar values (strings or numbers), and blocks (that you may view as matrices).

###  Scalar Variables  

A scalar variable {{< code "var" >}} can be defined by:

{{< code_line  " var <nowiki>=</nowiki> exp" >}}

{{< code "var" >}} can contain any alphanumeric character plus _, and {{< code "exp" >}} can be a quote-delimited string, a number (integer, real, or complex), a variable name, or a mathematical expression. Complex numbers are defined as <tt>{real, imag}</tt>. Real numbers can use scientific notation with <tt>e</tt> or <tt>E</tt> (no <tt>d</tt> or <tt>D</tt>, Fortran people), such as <tt>6.02e23</tt>. Variable names are not case-sensitive, and you must not redefine a previously defined symbol -- especially not the reserved variables <tt>x, y, z, r, w, t</tt> which are used in space- or time-dependent expressions, where <tt>w</tt> is the 4th space coordinate when operating in 4D.

###  Mathematical expressions  

The parser can interpret expressions in the input file, either to assign the result to a variable, or for defining functions such as a potential in the {{< Variable2 "Species" >}} block or a time-dependent function in the {{< Variable2 "TDFunctions" >}} block. The arguments can be numbers or other variables.

####  Arithmetic operators  
; {{< code "a+b" >}}: addition
; {{< code "a-b" >}}: subtraction
; {{< code "-a" >}}: unary minus
; {{< code "a*b" >}}: multiplication
; {{< code "a/b" >}}: division
; {{< code "a^b" >}}: exponentiation

####  Logical operators  

Logical operation will return 0 for false or 1 for true. You can exploit this to define a piecewise expression, e.g. <tt>"2 * (x <= 0) - 3 * (x > 0)"</tt> (although the <tt>step</tt> function may also be used). The comparison operators (except <tt>==</tt>) use only the real part of complex numbers.

; {{< code "a < b" >}}: less than
; <tt>a <= b</tt>: less than or equal to ($\le$)
; {{< code "a > b" >}}: greater than
; <tt>a >= b</tt>: greater than or equal to ($\ge$)
; <tt>a == b</tt>: equal to
; {{< code "a && b" >}}: logical and
; <tt>a || b</tt>: logical or
; {{< code "!a" >}}: logical not

####  Functions  

;{{< code "sqrt(x)" >}}: The square root of {{< code "x" >}}.
;{{< code "exp(x)" >}}: The exponential of {{< code "x" >}}.
;{{< code "log(x)" >}} or {{< code "ln(x)" >}}: The natural logarithm of {{< code "x" >}}.
;{{< code "log10(x)" >}}: Base 10 logarithm of {{< code "x" >}}.
;{{< code "logb(x, b)" >}}: Base {{< code "b" >}} logarithm of {{< code "x" >}}.

;{{< code "{x, y&-125;" >}}: The complex number $x + iy$.
;{{< code "arg(z)" >}}: Argument of the complex number {{< code "z" >}}, $\arg(z)$, where $-\pi < \arg(z) <= \pi$.
;{{< code "abs(z)" >}}: Magnitude of the complex number {{< code "z" >}}, $|z|$.
;{{< code "abs2(z)" >}}: Magnitude squared of the complex number {{< code "z" >}}, $|z|^2$.
;{{< code "logabs(z)" >}}: Natural logarithm of the magnitude of the complex number {{< code "z" >}}, $\log|z|$. It allows an accurate evaluation of $\log|z|$ when $|z|$ is close to one. The direct evaluation of {{< code "log(abs(z))" >}} would lead to a loss of precision in this case.
;{{< code "conjg(z)" >}}: Complex conjugate of the complex number {{< code "z" >}}, $z^* = x - i y$.
;{{< code "inv(z)" >}}: Inverse, or reciprocal, of the complex number {{< code "z" >}}, $\frac{1}{z} = \frac{x - i y}{x^2 + y^2}$.

;{{< code "sin(x)" >}}, {{< code "cos(x)" >}}, {{< code "tan(x)" >}}, {{< code "cot(x)" >}}, {{< code "sec(x)" >}}, {{< code "csc(x)" >}}: The sine, cosine, tangent, cotangent, secant and cosecant of {{< code "x" >}}.
;{{< code "asin(x)" >}}, {{< code "acos(x)" >}}, {{< code "atan(x)" >}}, {{< code "acot(x)" >}}, {{< code "asec(x)" >}}, {{< code "acsc(x)" >}}: The inverse (arc-) sine, cosine, tangent, cotangent, secant and cosecant of {{< code "x" >}}.
;{{< code "atan2(x,y)" >}}: = $\mathrm{atan}(y/x)$.
;{{< code "sinh(x)" >}}, {{< code "cosh(x)" >}}, {{< code "tanh(x)" >}}, {{< code "coth(x)" >}}, {{< code "sech(x)" >}}, {{< code "csch(x)" >}}: The hyperbolic sine, cosine, tangent, cotangent, secant and cosecant of {{< code "x" >}}.
;{{< code "asinh(x)" >}}, {{< code "acosh(x)" >}}, {{< code "atanh(x)" >}}, {{< code "acoth(x)" >}}, {{< code "asech(x)" >}}, {{< code "acsch(x)" >}}: The inverse hyperbolic sine, cosine, tangent, cotangent, secant and cosecant of {{< code "x" >}}.

;{{< code "min(x, y)" >}}: The minimum of {{< code "x" >}} and {{< code "y" >}}.
;{{< code "max(x, y)" >}}: The maximum of {{< code "x" >}} and {{< code "y" >}}.
;{{< code "step(x)" >}}: The Heaviside step function in {{< code "x" >}}. This can be used for piecewise-defined functions.
;{{< code "erf(x)" >}}: The error function $\mathrm{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x dt e^{-t^2}$.

;{{< code "realpart(z)" >}}: The real part of the complex number {{< code "z" >}}.
;{{< code "imagpart(z)" >}}: The imaginary part of the complex number {{< code "z" >}}.
;{{< code "floor(x)" >}}: The largest integer less than the real number {{< code "x" >}}.
;{{< code "ceiling(x)" >}}: The smallest integer greater than the real number {{< code "x" >}}.

These mathematical operations are all based on the GSL library and are defined in {{< code "symbols.c" >}} and {{< code "grammar.y" >}}.

####  References  
* https://www.gnu.org/software/gsl/manual/html_node/Properties-of-complex-numbers.html
* https://www.gnu.org/software/gsl/manual/html_node/Complex-arithmetic-operators.html
* https://www.gnu.org/software/gsl/manual/html_node/Error-Function.html
* https://www.gnu.org/software/libc/manual/html_node/Rounding-Functions.html
* https://www.gnu.org/software/gsl/manual/html_node/Representation-of-complex-numbers.html

###  Predefined variables  

There are some predefined constants for your convenience:

;{{< code "pi" >}}: {{< code "3.141592653589793" >}}.
;{{< code "e" >}}: The base of the natural logarithms.
;{{< code "false" >}} or {{< code "no" >}}: False.
;{{< code "true" >}} or {{< code "yes" >}}: True.
;{{< code "i" >}}: The imaginary unit $i$, ''i.e.'' {{< code "{0, 1" >}}}

Since version 7.0 there are also some predefined units that can be found by searching for the decimal point {{< code "." >}} in {{< code "share/variables" >}}:

;{{< code "angstrom" >}}: {{< code "1.8897261328856432" >}}.
;{{< code "pm" >}} or {{< code "picometer" >}}: {{< code "0.018897261328856432" >}}.
;{{< code "nm" >}} or {{< code "nanometer" >}}: {{< code "18.897261328856432" >}}.
;{{< code "ry" >}} or {{< code "rydberg" >}}: {{< code "0.5" >}}.
;{{< code "eV" >}} or {{< code "electronvolt" >}}: {{< code "0.03674932539796232" >}}.
;{{< code "invcm" >}}: {{< code "4.5563353e-06" >}}.
;{{< code "kelvin" >}}: {{< code " 3.1668105e-06" >}}.
;{{< code "kjoule_mol" >}}: {{< code " 0.00038087988" >}}.
;{{< code "kcal_mol" >}}: {{< code " 0.0015936014" >}}.
;{{< code "as" >}} or {{< code "attosecond" >}}: {{< code "0.0413413737896" >}}.
;{{< code "fs" >}} or {{< code "femtosecond" >}}: {{< code "41.3413737896" >}}.
;{{< code "ps" >}} or {{< code "picosecond" >}}: {{< code "41341.3737896" >}}.
;{{< code "c" >}}: {{< code "137.035999139" >}}.

###  Blocks  

Blocks are defined as a collection of values, organised in row and column format.
The syntax is the following:

{{code_line|%var}}
{{code_line|<nowiki>  exp | exp | exp | ...</nowiki>}}
{{code_line|<nowiki>  exp | exp | exp | ...</nowiki>}}
{{< code_line "  ..." >}}
{{code_line|%}}

Rows in a block are separated by a newline, while columns are separated by the character | or by a tab. There may be any number of lines and any number of columns in a block. Note also that each line can have a different number of columns. Values in a block don't have to be of the same type.

###  Comments  

Everything following the character {{< code "-" >}} until the end of the line is considered a comment and is simply cast into oblivion.

###  Includes  

With <code>include FILENAME</code> it is possible to include external files into the input file. To illustrate the usage of this command we can split the input file from the [methane tutorial](../Methane molecule ) in two.

{| class="wikitable" style="background:white;"
!inp
!geometry.oct
|-
|
```text
 {{< Variable2 "CalculationMode" >}} = gs
 {{< Variable2 "UnitsOutput" >}} = eV_Angstrom
  
 {{< Variable2 "Radius" >}} = 3.5*angstrom
 {{< Variable2 "Spacing" >}} = 0.22*angstrom
 
 
 include geometry.oct
```
|
```text
 CH = 1.2*angstrom
 %{{< Variable2 "Coordinates" >}}
   "C" |           0 |          0 |           0
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
```
|} 

###  Environment variables  

You can also [[Manual:Advanced_ways_of_running_Octopus-Passing_arguments_from_environment_variables | set variables using the environment]], which can be helpful in scripting.

###  Default values  

If {{< octopus >}} tries to read a variable that is not defined in the input file, it automatically assigns to it a default value (there are some cases where {{< octopus >}} cannot find a sensible default value and it will stop with an error). All variables read (present or not in the input file) are output to the file {{< file "exec/parser.log" >}}. The variable that are not defined in the input file will have a <tt>-default</tt> comment to it. If you are not sure of what the program is reading, just take a look at it.

We recommend you to keep the variables in the input file to a minimum: ''do not write a variable that will be assigned its default value''. The default can change in newer versions of {{< octopus >}} and old values might cause problems. Besides that, your input files become difficult to read and understand.

###  Documentation  

Each input variable has (or should have) its own documentation explaining what it does and the valid values it may take. This documentation can be obtained [http://octopus-code.org/doc/{{< octopus_version >}}/html/vars.php online] or it can also be accessed by the {{< Manual "External_utilities:oct-help" "oct-help" >}} command.

###  Experimental features  

Even in the stable releases of {{< Octopus >}} there are many features that are being developed and are not suitable for production runs. To protect users from inadvertly using these parts they are declared as ''Experimental''.

When you try to use one of these experimental functionalities {{< octopus >}} will stop with an error. If you want to use it you need to set the variable {{< Variable2 "ExperimentalFeatures" >}} to {{< code "yes" >}}. Now {{< octopus >}} will only emit a warning.

By setting {{< Variable2 "ExperimentalFeatures" >}} to {{< code "yes" >}} you will be allowed to use parts of the code that are not complete or not well tested and most likely produce wrong results. If you want to use them for production runs you should contact the Octopus developers first.

###  Good practices  

In order to ensure compatibility with newer versions of {{< octopus >}} and avoid problems, keep in mind the following rules of good practice when writing input files:

* Although input variables that take an option as an input can also take a number, the number representation makes the input file less readable and it is likely to change in the future. So '''avoid using numbers instead of values'''. For example
```text
 UnitsOutput = ev_angstrom
```
''must always'' be used instead of
```text
 UnitsOutput = 3
```
* '''Do not include variables that are not required in the input file''', especially declarations of values that are just a copy of the default value. This makes the input file longer, less readable and, since defaults are likely to change, it makes more probable that your input file will have problems with newer versions of {{< octopus >}}. Instead rely on default values.

* '''Avoid duplicating information in the input file'''. Use your own variables and the mathematical-interpretation capabilities for that. For example, you should use:
```text
 m = 0.1
 c = 137.036
 E = m*c^2
```
instead of
```text
 m = 0.1
 c = 137.036
 E = 1877.8865
```
In the second case, you might change the value of {{< code "m" >}} (or {{< code "c" >}} if you are a cosmologist) while forgetting to update {{< code "E" >}}, ending up with an inconsistent file.
```text
 
```
{{< manual_foot prev="Manual:Installation" next="Manual:Running Octopus" >}}
---------------------------------------------
