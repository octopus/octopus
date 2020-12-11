---
title: "Advanced ways of running Octopus"
series: "Manual"
---


##  Parallel {{< Octopus >}}  

{{< octopus >}} can be run in multi-tentacle (parallel) mode if you have configured it that way when you compiled and installed it.  There are four possible parallelization strategies: domains, states, spin/''k''-points, and electron-hole pairs (only for Casida).  You do not need to specify in the input file what kind of parallelization strategy you want to use as {{< octopus >}} can find reasonable defaults. However, the experienced user may want to tune the corresponding variables.

### Parallel in States 
The first strategy is parallel in states. This means that each processor takes some states and propagates them independently of the others. 

To use this include 

```text
 {{< Variable2 "ParStates" >}} = auto
```

in your input file.

This is clearly the most efficient (and
simplest) way of parallelizing the code if the number of states is much
larger than the number of processors available (which is usually the
case for us). You can expect
a more or less linear speed-up with number of processors until you reach
around 5 states/processor. Beyond this it starts to saturate.

### Parallel in Domains 

When running parallel in "domains", {{< octopus >}} divides the simulation region (the box) into separate regions (domains) and assigns each of of these to a different processor. This allows it not only to speed up the calculation, but also to divide the memory among the different processors. The first step of this process, the splitting of the box, is in general a very complicated process. Note that we are talking about a simulation box of an almost arbitrary shape, and of an arbitrary number of processors. Furthermore, as the communication between processors grows proportionally to the surface of the domains, one should use an algorithm that divides the box such that each domain has the same number of points, and at the same time minimizes the total area of the domains. {{< octopus >}} uses the [http://www-users.cs.umn.edu/~karypis/metis/ METIS] library to make this domain decomposition of the space, METIS is included in the source code of {{< octopus >}} and will be used by default.

So for example, if you are planning to run the ground state of Na there's only one choice:

```text
 {{< Variable2 "ParDomains" >}} = auto
```

{{< octopus >}} will then partition the real-space mesh and compute the resulting
pieces on different nodes. Later when you plan to do time propagations
you have more choices at hand. You could run only parallel in states

```text
 {{< Variable2 "ParStates" >}} = auto
 {{< Variable2 "ParDomains" >}} = no
```

or in domains (as above), or you can decide to employ both parallelization
strategies

```text
 {{< Variable2 "ParStates" >}} = auto
 {{< Variable2 "ParDomains" >}} = auto
```

### Problems with parallelization 

The parallelization is only effective if a sufficiently large amount of
load is given to each node. Otherwise communication dominates. In that
case one should do the calculation in serial. This happens, for example,
if you attempt to do parallelization in real-space domains, but your grid
does not have too many points. If the code detects that the number of
points per node is not large enough, it gives the warning:

```text
 ** Warning:
 ** From node =    0
 ** I have less elements in a parallel group than recommended.
 ** Maybe you should reduce the number of nodes 
```

### Compiling and running in parallel 
The steps to compile and use {{< octopus >}} in parallel are:

1) use the --enable-mpi in the configuration script 

2) give to {{< octopus >}} the libraries, etc. necessary to compile a parallel
code. When you install mpich in your system, it creates a series of
programs called mpicc, mpif77, mpif90, etc. These are wrappers that call
the specific compiler with the appropriate flags, libraries, etc. This
means that to compile {{< octopus >}} in parallel you'll need mpif90. Please check
that this is present in your system (perhaps in /usr/local/bin if you
compiled mpich from the sources). {{< octopus >}} will try to find this compiler,
but it may not be able to find it for several reasons:

* mpif90 is _not_ in the path.

* there is another f90 compiler specified by the user/system. This means that if you have the variable FC  defined when you run the configure script, {{< octopus >}} will use the fortran compiler specified by that variable, and will not try to find mpif90. Then, the configuration fails.
```text
 
 export FC=<whatever path>/mpif90
 export FCFLAGS=<whatever flags>
 ./configure --enable-mpi
```

It is wise to also set the <tt>FCFLAGS</tt>, as {{< octopus >}} cannot tell which fortran
compiler you are using if you compile with mpif90. 

3)  Now you just have to run {{< octopus >}} in parallel (this step depends on your actual system, you may have to use mpirun or mpiexec to accomplish it). 

In some run modes (e.g., <tt>td</tt>), you can use the multi-level parallelization, i.e., to run in parallel in more than one way at the same time. In the <tt>td</tt> case, you can run parallel in states and in domains at the same time. In order to fine-tune this behavior, please take a look at the variables {{< Variable2 "ParStates" >}}, {{< Variable2 "ParDomains" >}}, {{< Variable2 "ParKPoints" >}}, and {{< Variable2 "ParOther" >}}. In order to check if everything is OK, take a look at the output of {{< octopus >}} in the "Parallelization" section. This is an example:

```text
 ************************** Parallelization ***************************
 Octopus will run in *parallel*
 Info: Number of nodes in par_states  group:     8 (      62)
 Info: Octopus will waste at least  9.68% of computer time
 **********************************************************************
```

In this case, {{< octopus >}} runs in parallel only in states, using 8 processors (for 62 states). Furthermore, some of the processors will be idle for 9.68% of the time (this is not that great, so maybe a different number of processors would be better in this case).

##  Passing arguments from environment variables   

{{< octopus >}} can also read input variables from the environment variables of the shell. This is especially useful if you want to call {{< octopus >}} from scripts or if you want to pass one-time arguments when running {{< octopus >}}.

To activate this feature you must first set {{< code "OCT_PARSE_ENV" >}} to some value, for example in sh/bash

{{< command_line "<nowiki>export OCT_PARSE_ENV=1</nowiki>" >}}

After that, to pass an input variable to {{< octopus >}}, you must prepend the name of the variable by {{< code "'OCT_'" >}}, for example:

{{< command_line "<nowiki>export OCT_Radius=10</nowiki>" >}}

After that you can run {{< octopus >}} as usual.

You can also pass the variables in the same command line:

{{< command_line "<nowiki>OCT_PARSE_ENV=1 OCT_Radius=10 octopus</nowiki>" >}}

If you set a variable both from the input file and as a environment variable, the environment variable takes precedence. Be careful with this behaviour. Blocks are not supported (suggestions as how to elegantly put a block inside of environment variables are welcome).

There is an additional environment variable <tt>OCTOPUS_SHARE</tt> that can be set to specify the location of the <tt>variables</tt> file used by the parser. If this environment variable is not set, the default is as set by the configure script, namely <tt>prefix/share/octopus</tt>. You can set this to something else if necessary, for example to handle moving an executable from one machine to another, where the original path does not exist on the new machine.

###  Examples  

For example this is a simple script that helps to determine the optimal Radius of the simulation box (an {{< file "inp" >}} file with the rest of the parameters has to be provided):

```text

export OCT_PARSE_ENV=1
for OCT_Radius in `seq 5 16`
do
  export OCT_Radius
  octopus >& out-$OCT_Radius
  energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
  echo $OCT_Radius $energy
done
```
</pre>

{{< manual_foot prev="Manual:Visualization" next="Manual:External utilities:oct-analyze_projections" >}}
---------------------------------------------
