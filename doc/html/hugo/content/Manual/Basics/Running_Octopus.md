---
title: "Running Octopus"
series: "Manual"
weight: 1
description: "The command line"
---


###  Input  

In order to run, {{< octopus >}} requires an input file that ''must'' be called {{< file "inp" >}}. Depending on your input file there are other files that you might need, like pseudopotentials or coordinate files (we will discuss this later in this manual).

The rest of the files that are required are part of the {{< octopus >}} installation; if {{< octopus >}} is correctly installed they should be available and {{< octopus >}} should know where to find them. With {{< octopus >}} you can't just copy the binary between systems and expect it to work.

###  Executing  

To run {{< octopus >}} you just have to give the {{< command "octopus" >}} command in the directory where you have your input file. While running, {{< octopus >}} will display information on the screen. If you want to capture this information you can send the output to a file, let's say {{< file "out.log" >}}, by executing it like this:

<!--{{% command_line "octopus > out.log" %}}-->
```bash
octopus > out.log
```

This captures only the normal output. If there is a warning or an error, it will be still printed on the screen. To capture everything to a file, run

<!--{{< command_line "octopus >& out.log" >}}-->
```bash
octopus >& out.log
```

If you want to run {{< octopus >}} in the background, append {{< command "&" >}} to the last command.

###  Output  

While running, {{< Octopus >}} will create several output files, all of them inside subdirectories in the same directory where it was run. The files that contain the physical information depend on what {{< octopus >}} is doing and they will be discussed in the next chapter. 

One directory that is always created is {{< file "exec/" >}}, this file contains information about the {{< octopus >}} run. Here you will find the {{< file "parser.log" >}} file, a text file that contains {{< emph "all" >}} the input variables that were read by {{< octopus >}}, both the variables that are in the input file and the variables that took default values; in the second case they are marked by a comment as {{< code "-default" >}}. This file is very useful if you want to check that {{< octopus >}} is correctly parsing a variable or what are the default values that it is taking.

###  Clean Stop  

If you create a file called {{< file "stop" >}} in the running directory, {{< octopus >}} will exit gracefully after finishing the current iteration. May not work for gcm, invert_ks, casida run modes. You can use this to prevent possible corruption of restart information when your job is killed by a scheduler, by preemptively asking {{< octopus >}} to quit automatically in a job script like this:

```text
 -PBS -l walltime=4:10:00
 mpirun $HOME/octopus/bin/octopus &> output &
 sleep 4h
 touch stop
 wait
```

or more sophisticatedly like this:

```text
 MIN=`qstat -a $PBS_JOBID | awk '{wall=$9} END {print $9}' | awk -F: '{min=($1*60)+$2; print min}'`
 sh ~/sleep_stop.sh stop $((MIN-10))m > sleepy &
 mpirun $HOME/octopus/bin/octopus &> output
```

with auxiliary script {{< code "sleep_stop.sh" >}}:

```text
 -!/bin/bash
 if [ $- -ne 2 ]; then
    echo "Usage: sleep_stop.sh FILENAME TIME"
 else
    echo "Time until $1 created: $2"
    rm -f $1
    sleep $2
    touch $1
 fi
```

###  Restarting  

Another directory that will created is {{< file "restart/" >}}; in this directory {{< octopus >}} saves the information from the calculation that it is doing. This information can be used in the following cases:

* If {{< octopus >}} is stopped without finishing by some reason, it can restart from where it was without having to do all work again.
* If after the calculation is done (or even if it was stopped), the user wants to do the same simulation with some different parameters, {{< octopus >}} can save some work by starting from the restart information.
* There are some calculations that require the results of other type of calculation as an input; in this case it uses the files written in {{< file "restart/" >}} by the previous calculation (we will discuss this case later, when we talk about the different calculation modes).

Sometimes it's not desired to restart a calculation, but to start it from the very beginning. {{< octopus >}} can be instructed to do so by setting the input variable {{< Variable2 "fromScratch" >}} to {{< code "yes" >}}.

{{< manual_foot prev="Manual:Input file" next="Manual:Units" >}}
---------------------------------------------
