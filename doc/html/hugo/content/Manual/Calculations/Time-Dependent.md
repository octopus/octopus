---
title: "Time-Dependent"
series: "Manual"
weight: 2
description: "Time evolution"
---


###  Time evolution  

When {{< Variable2 "CalculationMode" >}} = {{< code "td" >}}, the code performs the time-propagation of the electronic orbitals and &ndash; if required &ndash; the ionic positions. The latter task does not pose major algorithmical problems (the usual Verlet algorithms deal with that task); however the best way to propagate a Schr&ouml;dinger-like equation is still unclear. Due to this fact, we provide with a rather excessive selection of possibilities for that purpose. Before describing the set of variables necessary to specify the way in which the time evolution is to be performed, it is worth making a brief introduction to the problem.

We are concerned with a set of Schrödinger-like equations for the electronic orbitals $\psi_j(t)\,\!$:

$$
i{\\partial \\over \\partial t}\\psi\_j(t) = H(t)\\psi\_j(t)\\,\\!
$$

$$
\\psi\_j(t=0) = \\psi\_j^{(0)}\\,\\!
$$

Because this equation is linear (the time derivative and the Hamiltonian are both linear operators), one may formally define a linear “evolution” operator, $U(T, t)\,\!,$ which transforms the initial vector into the solution at time $T\,\!$:

$$
\\psi\_j(T) = U(T, 0)\\psi\_j^{(0)}
$$

Moreover, there is the formally exact expression for the evolution operator

$$
\\psi\_j(T) = \\mathcal{T}\\!\\!\\exp\\left\\{ -i\\!\\!\\int\_0^{T}d\\tau H(\\tau)\\right\\} \\psi\_j^{(0)}\\,\\!,
$$

where $\mathcal{T}\!\!\exp\,\!$ is the time-ordered exponential, which is a short-hand for:

$$
\\psi\_j(T) = \\left\\{\\sum\_{n=0}^{\\infty} \\frac{\\left(-i\\right)^n}{n!} \\int\_0^t d\\tau\_1 \\cdots \\int\_0^t d\\tau\_n H(\\tau\_1) \\cdots H(\\tau\_n)\\right\\}\\psi\_j^{(0)}\\,\\!
$$

If the Hamiltonian commutes with itself at different times, we can drop the time-ordering product, and leave a simple exponential. If the Hamiltonian is time-independent, which makes it trivially self commuting, the solution is simply written as:

$$
\\psi\_j(T) = \\exp\\left\\{ -iTH\\right\\} \\psi\_j^{(0)}\\,\\!.
$$

Unfortunately, this is not the case for TDDFT when the system is exposed to external time-dependent perturbations like electric and magnetic fields or pulsed lasers. But even without an external time-dependency, there remains the intrinsic time-dependency of the Kohn-Sham Hamiltonian, which is built &ldquo;self-consistently&rdquo; from the varying electronic density.

The first step to tackle this problem is to split the propagation of the long interval $[0, T]\,\!$ into $N\,\!$ smaller steps by utilizing the group property

$$
U(T, t) = U(T, t')U(t', t)\\,\\!
$$

of the time-evolution operator. This yields the following time discretization:

$$
U(T, 0) = \\prod\_{i=0}^{N-1}U(t\_i+\\Delta t, t\_i)\\,\\!,
$$

where $t_0=0\,\!$, $t_N=T\,\!$, $\Delta t = T/N\,\!$. So at each time step we are dealing with the problem of performing the short-time propagation:

$$
\\psi\_j(t+\\Delta t) = U(t+\\Delta t, t)\\psi\_j(t) = \\mathcal{T}\\!\\!\\exp\\left\\{ -i\\int\_{t}^{t+\\Delta t}d\\tau H(\\tau)\\right\\} \\psi\_j(t)\\,\\!.
$$

In this way, one can monitor the evolution in the interior of $[0, T]\,\!$. In fact, the possibility of monitoring the evolution is generally a requirement.

This requirement imposes a natural restriction on the maximum size of $\Delta t\,\!$: if $\omega_{\mathrm{max}}\,\!$ is the maximum frequency that we want to discern, $\Delta t\,\!$ should be no larger than $\approx 1/\omega_{\mathrm{max}}\,\!$. Below this $\Delta t_{\mathrm{max}}\,\!$, we are free to choose $\Delta t\,\!$ considering performance reasons: technically, the reason for the discretization is two-fold: the time-dependence of $H\,\!$ is alleviated, and the norm of the exponential argument is reduced (the norm increases linearly with $\Delta t\,\!$).

Since we cannot drop the time-ordering product, the desired algorithm cannot be reduced, in principle, to the calculation of the action of the exponential of an operator over the initial vector. Some algorithms tailored to approximate the evolution operator, in fact, do not even need to peform such operator exponentials. Most of them, however, do rely on the calculation of one or more exponentials, such as the ones used by {{< Octopus >}}. This is why in principle we need to specify two different issues: the &ldquo;evolution method&rdquo;, and the &ldquo;exponential method&rdquo;. In other words: we need an algorithm to approximate the evolution operator $U(t+\Delta{t},t)\,\!$ &ndash; which will be specified by variable {{< Variable2 "TDPropagator" >}} &ndash; and, if this algorithm requires it, we will also need an algorithm to approximate the exponential of a matrix operator $\exp\left\{ A\right\}\,\!$ &ndash; which will be specified by variable {{< Variable2 "TDExponentialMethod" >}}.

####  Propagators for the time-dependent Kohn-Sham equations  

In the following, we describe the propagators available in {{< Octopus >}} for the time-dependent Kohn-Sham equations. These propagators solve the problem of approximating the orbitals $\psi_j(t+\Delta t)\,\!$ from the knowledge of $\psi_j(\tau)\,\!$ and $H(\tau)\,\!$for $0\le\tau\le t\,\!$. Some methods require the knowledge of the Hamiltonian at some points $\tau\,\!$ in time between $t\,\!$ and $t+\Delta t\,\!$.

This quantity can be approximated by the following procedure:

- Approximate $H(\tau)\,\!$ through extrapolation of a polynomial fit to a certain number of previous time steps.
- Propagate $\psi_j(t)\,\!$ to get $\psi_j(t+\Delta t)\,\!$.
- Calculate $H(t+\Delta t)\,\!$ from the orbitals $\psi_j(t+\Delta t)\,\!$.
- Interpolate the required $H(\tau)\,\!$ from $H(t)\,\!$ and $H(t+\Delta t)\,\!$.
- Repeat the steps 2 to 4 until self-consistency is reached.

In {{< Octopus >}}, however, the above scheme is dropped for performance reasons and only step 1 is implemented, via a second-order extrapolation, except for the first two steps where the extrapolation obviously cannot be trusted. Instead, we rely on a sufficiently small $\Delta t\,\!$.

#####  Midpoint rules  

The implicit midpoint rule or Crank-Nicolson method ({{< Variable2 "TDPropagator" >}}={{< code "crank_nicolson" >}}) calculates the exponential of the Hamiltonian by a first-order Pad&eacute; approximation. The Hamiltonian is evaluated at $t + \Delta t/2\,\!$ and the integral is dropped:

$$
U\_{\\mathrm{CN}}(t + \\Delta t, t) = \\frac{1- i \\frac{\\Delta t}{2}H(t+\\Delta t/2)}{1+ \\frac{\\Delta t}{2}H(t+\\Delta t/2)}\\,\\!
$$

The calculation of the matrix fraction is transformed into the solution of the linear system

$$
L \\psi\_j(t+\\Delta t) = b\\,\\!
$$

with the known quantities

$$
L = 1+i\\frac{\\Delta t}{2}H(t+\\Delta t/2)\\,\\!
$$

and

$$
b = \\left\\{1-i\\frac{\\Delta t}{2}H(t+\\Delta t/2)\\right\\}\\psi\_j(t).\\,\\!
$$

The Crank-Nicolson scheme is unitary and preserves time-reversal symmetry.

A simpler but similar scheme is the exponential midpoint rule
({{< Variable2 "TDPropagator" >}}={{< code "exp_mid" >}}), which is unitary and time-reversal-symmetry-preserving, but has the drawback that it requires fairly small time steps:

$$
U\_{\\mathrm{EM}}(t + \\Delta t, t) = \\exp\\left\\{-i\\Delta t H(t + \\Delta t/2)\\right\\}\\,\\!
$$

#####  Magnus expansions  

Magnus expansion ({{< Variable2 "TDPropagator" >}}={{< code "magnus" >}}) is the most sophisticated propagation implemented in {{< Octopus >}}. According to Magnus, there exists an operator $\Omega(t+\Delta t, t)\,\!$ for which holds

$$
U(t+\\Delta t, t) = \\exp\\left\\{\\Omega(t+\\Delta t, t)\\right\\}\\,\\!
$$

given by the series

$$
\\Omega(t + \\Delta t, t) = \\sum\_{k=1}^\\infty \\Omega\_k(t + \\Delta t, t)\\,\\!
$$

that converges at least for some local environment of $t\,\!$.

The operators $\Omega_k(t+\Delta t, t)\,\!$ are generated with the recursive procedure

$$
\\Omega\_k(t+\\Delta t, t) = \\sum\_{l=0}^{k-1} \\frac{B\_l}{l!}\\int\_t^{t+\\Delta t} S^l\_k(\\tau)d\\tau,\\,\\!
$$

$$
S\_1^0(\\tau) = -iH(\\tau), \\quad S\_k^0=0\\ \\mathrm{for}\\ k\>1,\\,\\!
$$

$$
S\_k^j(\\tau) = \\sum\_{m=1}^{k-l}\\left\[\\Omega\_m(t+\\Delta t, t), S\_{k-m}^{l-1}(\\tau)\\right\]\\ \\mathrm{for}\\ 1\\le l\\le k-1,\\,\\!
$$

where $B_l\,\!$ are Bernoulli numbers.

In {{< Octopus >}} we have implemented a fourth-order Magnus expansion, which means that the operator series is truncated to second order and the integrals are calculated by a second-order quadrature formula. The resulting operator is

$$
\\Omega\_{M(4)}(t + \\Delta t, t) = -i\\frac{\\Delta t}{2}\\left\[H(t\_1)+ H(t\_2)\\right\] - \\frac{\\sqrt{3}\\Delta t^2}{12}\\left\[H(t\_2), H(t\_1)\\right\]\\,\\!
$$

with the quadrature sampling points $t_{1,2} = t + \left[\frac{1}{2}\mp \frac{\sqrt{3}}{6}\right]\Delta t\,\!$.

With a modified Hamiltonian

$$
H\_{M(4)}(t, \\Delta t) = \\overline{H}(t, \\Delta t) + i \\left\[T + v\_{\\mathrm{ext}}^{\\mathrm{nonlocal}}(t), \\overline{\\Delta V}(t, \\Delta t)\\right\],\\,\\!
$$

$$
\\overline{H}(t, \\Delta t) = T + \\frac{1}{2}\\left\[V\_{\\mathrm{KS}}(t\_1)+V\_{\\mathrm{KS}}(t\_2)\\right\],\\,\\!
$$

$$
\\overline{\\Delta V}(t, \\Delta t) = \\frac{\\sqrt{3}}{12}\\Delta t\\left\[V\_{\\mathrm{KS}}(t\_1)-V\_{\\mathrm{KS}}(t\_1)\\right\]\\,\\!
$$

the propagator takes the form

$$
U\_{M(4)}(t+\\Delta t, t) = \\exp\\left\\{-i \\Delta t H\_{M(4)}(t, \\Delta t)\\right\\}.\\,\\!
$$

Note that only the nonlocal components of the Kohn-Sham Hamiltonian contribute to the commutator and, furthermore, that we make the assumption that the nonlocal potential $v_{\mathrm{ext}}^{\mathrm{nonlocal}}(t)\,\!$ does not vary significantly in the interval $[t, t+\Delta t]\,\!$. This nonlocal component stems from the ionic pseudopotentials. Consequently, its variation is caused by the ionic movement, which is negligible on the electronic time scale determining $\Delta t\,\!$.

#####  Time-reversal-symmetry based propagation  

Due to time-reversal symmetry, propagating backwards by $\Delta t/2\,\!$ starting from $\psi_j(t+\Delta t)\,\!$ should lead to the same result as propagating forwards by $\Delta t/2\,\!$ starting from $\psi_j(t)\,\!$. Using the simplest approximation to the evolution operator, we end up with the condition

$$
\\exp\\left\\{+i\\frac{\\Delta t}{2}H(t +\\Delta t)\\right\\} \\psi\_j(t+\\Delta t) = \\exp\\left\\{-i\\frac{\\Delta t}{2}H(t)\\right\\} \\psi\_j(t)\\,\\!
$$

Rearranging the terms gives an approximation for the propagator:

$$
U\_{\\mathrm{ETRS}}(t + \\Delta t, t) = \\exp\\left\\{-i \\frac{\\Delta t}{2} H(t + \\Delta t)\\right\\}
\\exp\\left\\{-i \\frac{\\Delta t}{2}H(t)\\right\\}\\,\\!
$$

This ''enforced time-reversal symmetry'' method is available in {{< Octopus >}} by setting ({{< Variable2 "TDPropagator" >}}={{< code "etrs" >}}).

It is worthwhile to give a sidenote on the approximation of $H(t+\Delta t)\,\!$: for the {{< code "etrs" >}} method the Hamiltonian at time $t+\Delta t\,\!$ is not extrapolated, but rather built from the density $n' = \sum_j|\psi'_j(t+\Delta t)|^2\,\!$, which, in turn, is calculated from the orbital estimate

$$
\\psi'\_j(t+\\Delta t) = \\exp\\left\\{-i\\Delta t H(t)\\right\\}\\psi\_j(t).\\,\\!
$$

There exists, however, also a variant {{< Variable2 "TDPropagator" >}}={{< code "aetrs" >}} ''(approximated enforced time-reversal symmetry)'' that performs a second-order polynomial extrapolation of $H(t-k\Delta t)\,\!$, $k=0,1,2\,\!$, to get $H(t+\Delta t)\,\!$, which is about 40&nbsp;% faster than {{< code "etrs" >}}.

####  Approximations to the exponential of an operator  

Most of the evolution methods described in the previous section require the calculation of the exponential of a matrix. Before going into the details of the {{< Variable2 "TDExponentialMethod" >}} that {{< Octopus >}} provides, one general difficulty deserves mentioning.

In principle, we would like to calculate $\exp\left\{A\right\}v\,\!$ by first calculating $\exp\left\{A\right\}\,\!$ and then applying this result to any vector $v\,\!$. Unfortunately, the size of our Hamiltonian matrix is of the order $\approx 10^5\,\!$ which forbids its full storage in matrix form. Apart from that, methods to calculate the exponential of a matrix explicitly are limited to a few thousand matrix elements. For this reason, we have to turn to iterative solutions that calculate $\exp\left\{A\right\}v\,\!,$ for a particular vector $v\,\!$.

#####  Polynomial expansions  

The simplest possibility is to use the definition of the exponential of a matrix

$$
\\exp\\left\\{A\\right\\} = \\sum\_{k=0}^\\infty \\frac{1}{k!}A^k\\,\\!
$$

to get the approximation of order $k\,\!$({{< Variable2 "TDExponentialMethod" >}}={{< code "taylor" >}}):

$$
\\mathrm{taylor}\_N\\left\\{A, v\\right\\} = \\sum\_{k=0}^N\\frac{1}{k!}A^kv\\,\\!
$$

$\mathrm{taylor}_N\left\{A, v\right\}\,\!$ amounts to the expansion of the exponential in the standard polynomial base $\left\{1, x, x^2, \ldots\right\}\,\!$.

Experience shows that the choice $k=4\,\!$ gives particularly good results for our TDDFT implementation.

The standard polynomial basis is not the only possibility; {{< Octopus >}} also implements the expansion of the exponential in the Chebyshev basis ({{< Variable2 "TDExponentialMethod" >}}={{< code "chebyshev" >}}):

$$
\\mathrm{cheb}\_N\\left\\{A, v\\right\\} = \\sum\_{k=0}^Nc\_kT\_k(A)kv\\,\\!
$$

with $T_k\,\!$ being the Chebyshev polynomial of order $k\,\!$.

The truncation $N\,\!$ for both these expansions can be set by the input variable {{< Variable2 "TDExpOrder" >}}.

#####  Krylov-subspace projection  

The Krylov subspace of order $N\,\!$ for a given operator $A\,\!$ and vector $v\,\!$ is defined as

$$
\\mathcal{K}\_N\\left\\{A, v\\right\\} = \\mathrm{span}\\left\\{v, Av, A^2v, \\ldots, A^{N-1}v\\right\\}.\\,\\!
$$

The idea is to find the element of $\mathcal{K}_N\left\{A, v\right\}\,\!$ that optimally approximates $\exp\left\{A\right\}\,\!$ in a least-squares sense. It can be proven that this is given by

$$
\\beta V\_N(V\_N^T\\exp\\left\\{A\\right\\}V\_N)e\_1\\,\\!
$$

with $V_N = [v_1, \ldots, v_N]\,\!$ being an orthonormal basis of $\mathcal{K}_N\left\{A, v\right\}\,\!$ calculated with the Lanczos procedure. To get rid of the exponential, the approximation

$$
V\_N^T\\exp\\left\\{A\\right\\}V\_N \\approx \\exp\\left\\{V\_N^TAV\_N\\right\\}\\,\\!
$$

is introduced. The object $H_N=V_N^TAV_N\,\!$ is also a result of the Lanczos procedure and of the small dimension $k\times k\,\!$ which allows for direct calculation of the remaining exponential. So, we have

$$
\\mathrm{lanczos}\_N\\left\\{A, v\\right\\} = V\_N\\exp\\left\\{H\_N\\right\\}e\_1,\\,\\!\</math\>

which can be selected by setting {{< Variable2 "TDExponentialMethod" >}}={{< code "lanczos" >}}. When actually calculating the exponential, {{< Octopus >}} recursively generates the quantities \<math\>H\_N\\,\\!\</math\> and \<math\>V\_N\\,\\!\</math\> until some convergence criterion &ndash; controlled by {{< Variable2 "TDLanczosTol" >}} &ndash; is met or the order exceeds the maximum allowed order controlled by {{< Variable2 "TDExpOrder" >}}. In the latter case, a warning is written. Care should be taken in choosing the convergence parameter small enough to avoid faulty evolutions.

###  Performing a time-propagation  

To actually perform a time-propagation, a necessary prerequisite is to have the ground-state density (see \[\[Manual:Ground\_State|Ground State\]\]) of the system under consideration.

Then, the {{< Variable2 "CalculationMode" >}} has to be set to {{< code "td" >}}.

Besides the already mentioned input variables concerning the propagation method, two more important parameters have to be set: {{< Variable2 "TDTimeStep" >}} and {{< Variable2 "TDMaxSteps" >}}. The first one sets \<math\>\\Delta t\\,\\!\</math\>, the second sets the number of iterations. It is convenient to define these two variables like this:


 T  = 0.1     - Length of propagation.
 dt = 0.002   - Length of time step.
 
 {{< Variable2 "TDMaxSteps" >}} = T/dt
 {{< Variable2 "TDTimeStep" >}} = dt


In the absence of a perturbation, the system's total energy should be a constant. This check to determine a reasonable time step has to be done before any other calculation.

The relevant parts of {{< Octopus >}}' output might look like

\<pre\>
  Iter           Time        Energy     Elapsed Time
      1       0.002000   -541.542073         4.800
      2       0.004000   -541.542073         4.560
      3       0.006000   -541.542073         3.880
      .
      .
      .
     48       0.096000   -541.542073         7.420
     49       0.098000   -541.542073         7.620
     50       0.100000   -541.542073         7.490
\</pre\>

The energy is printed in the third column and should remain constant.
In general, larger time steps are desirable to shorten computational propagation time but keep in mind that the size of the time step limits to the maximum frequency you will be able to observe, and, consequently, also any external perturbation to which you might wish to expose the system.

###  External fields  



#####  Delta kick: Calculating an absorption spectrum  

To obtain the linear optical absorption spectrum of the system, we follow the scheme proposed by Yabana and Bertsch, and
excite all frequencies of the system by giving some small
momentum (\<math\>{\\mathcal K}\</math\>) to the electrons. This is achieved by transforming the ground-state wavefunctions according to:

$$
  \\varphi\_i(r, \\delta t) = e^{i {\\mathcal K} z} \\varphi\_i(r, 0)
  \\,,
$$

and then propagating these wavefunctions for some (finite) time.
The spectrum can then be obtained from the expression for the dipole
strength function $S(\omega)$:

$$
  S(\\omega) = \\frac{2 \\omega}{\\pi} {\\mathcal Im} \\alpha(\\omega)
  \\,,
$$

where the dynamical polarizability, $\alpha(\omega)$, is
essentially the Fourier transform of the dipole moment of the system
$d(t)$:

$$
  \\alpha(\\omega) = \\frac{1}{\\mathcal{K}}\\int\\!\\! dt \\; 
  e^{i\\omega t} \\left\[d(t)-d(0)\\right\]\\,.
$$

With this definition, the Thomas-Reiche-Kuhn ''f''-sum rule 
for the number of electrons, $N$, is given by the
integral:

$$
 N = \\int\\!\\! d\\omega\\;S(\\omega)\\,.
$$

This sum rule can be used to check the quality of the calculations. Another check
is energy conservation, which TDDFT respects when no external field is applied.

To obtain a spectrum with such a recipe in {{< Octopus >}} one has to follow the steps:

- Choose a system of linearly independent (can be non-orthogonal) axes (defined by the variable {{< Variable2 "TDPolarization" >}}). The defaults are the 3 Cartesian axes.
- Add an electric kick using the variable {{< Variable2 "TDDeltaStrength" >}}. Its value should be small so that we remain in the linear regime, but large enough to avoid numerical errors.
- Run a time-dependent simulation for the 3 independent directions. This involves running {{< octopus >}} 3 times, changing the value of {{< Variable2 "TDPolarizationDirection" >}} each time. After each run, move the file <tt>td.general/multipoles</tt> to <tt>td.general/multipoles.X<tt>, where <tt>X</tt> is 1, 2, or 3.
- Run the utility {{< Manual "External_utilities:oct-propagation_spectrum" "oct-propagation_spectrum" >}}.

#####  Lasers  

To calculate non-linear optical properties, we follow
the evolution of the system under the influence of a laser field that is treated 
in the dipole approximation (although this constraint can be removed). 
The harmonic emission spectrum can then be calculated from
the acceleration of the dipole moment:

$$
  H(\\omega) \\propto \\left| \\int\\!\\! dt \\; e^{i\\omega t} \\frac{d^2}{dt^2} d(t) \\right|^2\\,.
$$

During the propagation, charge density is absorbed at the boundaries
of the simulation region, either by an imaginary absorbing potential
or a mask function. In the first case, we add to the Kohn-Sham potential a term

$$
  V\_{\\rm eff}(r, t) = V\_{\\rm KS}(r,t) - i V\_{\\rm abs}(r)
  \\,,
$$

where $V_{\rm abs}$ is zero in the inner region of the simulation box, and
rises smoothly till the edges. By adjusting both the height and 
the shape of the potential, we can select which momenta are absorbed
and prevent the unwanted reflections at the boundary.
When using a mask, the wavefunction is multiplied in each time-step
by a function which is 1 in the inner simulation region and gradually 
goes to 0 at the borders:

$$
  \\varphi(r, t) \\rightarrow M(r) \\varphi(r, t)
  \\,.
$$

The absorbed charge can be interpreted as an ionization probability
and can be used to estimate the photo-electron spectra.
The box size has to be big enough so that the physical system is not
perturbed by the absorbing boundaries. Note that the wavefunctions are no 
longer normalized as the system slowly gets charged.


###  Symmetries  

The dynamic polarizability (trivially related to optical absorption) is, in its most general form, a 3x3 tensor. The reason is that we can shine light on the system polarized in any of the three Cartesian axes, and for each of these three cases measure how the dipole of the molecule oscillates along the three Cartesian axes. This usually means that to obtain the full dynamic polarizability of the molecule we usually need to apply 3 different perturbations along $x, y, z\,$, by setting {{< Variable2 "TDPolarizationDirection" >}} to 1, 2, or 3.

However, if the molecule has some symmetries, it is in general possible to reduce the total number of calculations from 3 to 2, or even 1. This is explained in detail [http://nano-bio.ehu.es/node/576 here]. To use this formalism in {{< octopus >}}, you can use the variables {{< Variable2 "TDPolarization" >}}, {{< Variable2 "TDPolarizationDirection" >}}, {{< Variable2 "TDPolarizationEquivAxes" >}}, and {{< Variable2 "TDPolarizationWprime" >}}.  The block {{< Variable2 "TDPolarization" >}} defines a basis (not necessarily orthogonal) chosen to maximize the number of equivalent axes, and {{< Variable2 "TDPolarizationDirection" >}}, and {{< Variable2 "TDPolarizationWprime" >}} are vectors specified in that basis.
Let us give a couple of examples.

The methane molecule has $T_d$ symmetry, which means that the response is identical for all directions. This means that we only need one propagation to obtain the whole tensor. This propagation can be performed in any direction we wish. So we could use the input

```text
 %{{< Variable2 "TDPolarization" >}}
  1 | 0 | 0
  0 | 1 | 0
  0 | 0 | 1
 %
 {{< Variable2 "TDPolarizationDirection" >}} = 1
 {{< Variable2 "TDPolarizationEquivAxes" >}} = 3
 %{{< Variable2 "TDPolarizationWprime" >}}
  0 | 0 | 1
 %
```

Note that we could have omitted the blocks {{< Variable2 "TDPolarization" >}} and {{< Variable2 "TDPolarizationWprime" >}} in the previous input file, as these are their default values. 

Now let us look at a linear molecule. In this case, you might think that we need two calculations to obtain the whole tensor, one for the direction along the axis of the molecule, and another for the axis perpendicular to the molecule. The fact is that we need only one, in a specially chosen direction, so that our field has components both along the axis of the molecule and perpendicular to it. Let us assume that the axis of the molecule is oriented along the $x$-axis. Then we can use

```text
 %{{< Variable2 "TDPolarization" >}}
  1/sqrt(2) | -1/sqrt(2) | 0
  1/sqrt(2) |  1/sqrt(2) | 0
  1/sqrt(2) |  0         | 1/sqrt(2)
 %
 {{< Variable2 "TDPolarizationDirection" >}} = 1
 {{< Variable2 "TDPolarizationEquivAxes" >}} = 3
 %{{< Variable2 "TDPolarizationWprime" >}}
  0 | 0 | 1
 %
```

You should try to convince yourself that the three axes are indeed equivalent and linearly independent.

Finally, let us look at a general planar molecule (in the ''xy'' plane). In principle we need only two calculations (reduced to one if more symmetries are present, as for benzene). In this case we chose one of the polarization axes on the plane, and the other two rotated 45 degrees out of plane:

```text
 %{{< Variable2 "TDPolarizationDirection" >}}
  1/sqrt(2) | 0 | 1/sqrt(2)
  1/sqrt(2) | 0 |-1/sqrt(2)
  0         | 1 | 0
 %
 
 {{< Variable2 "TDPolarizationEquivAxes" >}} = 2
```

In this case, we need two runs, one for {{< Variable2 "TDPolarizationDirection" >}} equal to 1, and another for it equal to 3. Note that if there are fewer than 3 equivalent axes, {{< Variable2 "TDPolarizationWprime" >}} is irrelevant.


{{< Manual_foot prev="Manual:Ground_State" next="Manual:Casida" >}}
---------------------------------------------
