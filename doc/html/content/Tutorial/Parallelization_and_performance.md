---
title: "Parallelization and performance"
tags: ["Tutorial", "Advanced"]
#series: "Tutorial"
---

{{< notice note >}}
This tutorial page was set up for the Benasque TDDFT school 2014. The specific references to the supercomputer used at that time will have to be adapted for others to use this tutorial.
{{< /notice >}}

In this tutorial we will see how to run a relatively big system in the [Hopper supercomputer](https://www.nersc.gov/users/computational-systems/hopper)  (at NERSC in California), and how to measure its performance. There are a few key things you need to know about how to interact with the machine. To log in, run {{< code "ssh trainX@hopper.nersc.gov" >}} in your terminal, substituting the actual name of your training account. Be aware that since this machine is far away, you should not try running X-Windows programs! You submit jobs by the {{< code "qsub" >}} command, ''e.g.'' {{< code "qsub job.scr" >}}, which will put them in the queue for execution when there is free space. You can see what jobs you currently have in the queue by executing {{< code "qstat -u $USER" >}}, so you can see when your job finishes. A status code will be shown: Q = waiting in the queue, R = running, C = complete. You can cancel a job by {{< code "qdel" >}} + the job number, as written by {{< code "qstat" >}}. The job script (''e.g.'' {{< code "job.scr" >}}) specifies parameters to the PBS/Torque queuing system about how many cores to use, what commands to run, etc.

#  Running the ground state  

We will need the input file, job submission script, coordinates file, and a pseudopotential for Mg (for the other elements we will use the default ones that come with Octopus). The pseudopotential is available in the Octopus directory at {{< code "jube/input/Mg.fhi" >}}. You can copy the files on hopper directly from {{< code "/global/homes/d/dstrubbe/octopus_tutorial" >}} to your scratch directory as follows:

```bash
 cd $SCRATCH
 cp -r /global/homes/d/dstrubbe/octopus_tutorial .
```

The input file ({{< file inp >}}):
{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 
 #### System size and parameters
 
 {{< variable "Spacing" >}} = 0.20
 {{< variable "Radius" >}} = 4.0
 
 {{< variable "Units" >}} = ev_angstrom
 {{< variable "XYZCoordinates" >}} = "xyz"
 %{{< variable "Species" >}}
  "Mg" | 24.305 | spec_ps_fhi | 12 | 3 | 2
 %
 {{< variable "ExcessCharge" >}} = 0
 
 {{< variable "XCFunctional" >}} = gga_x_pbe + gga_c_pbe
 
 {{< variable "ExtraStates" >}} = 18
 {{< variable "Eigensolver" >}} = rmmdiis
 {{< variable "LCAOAlternative" >}} = yes
 {{< variable "SmearingFunction" >}} = fermi_dirac
 {{< variable "Smearing" >}} = 0.1
 {{< variable "Mixing" >}} = 0.15
               
 #### GS
 {{< variable "MaximumIter" >}} = 300
 {{< variable "EigensolverTolerance" >}} = 1e-8
 {{< variable "ConvRelDens" >}} = 5e-8
 
 #### Saving memory
                    
 {{< variable "SymmetriesCompute" >}} = no
 {{< variable "PartitionPrint" >}} = no          
 {{< variable "MeshPartitionPackage" >}} = metis
 
 # Additional options
 {{< variable "ExperimentalFeatures" >}} = yes
{{< /code-block >}}

Submission script {{< file job.scr >}} with 24 CPU processor cores: 
```bash
 #!/bin/bash
 #PBS -q regular
 #PBS -l mppwidth=24
 #PBS -l advres=benasque.348
 #PBS -l walltime=0:30:00
 #PBS -N testing_chl
 #PBS -V
 
 module load octopus/4.1.2
 cd $PBS_O_WORKDIR
 aprun -n 24 octopus_mpi &> output_gs_24
```

To run:
```bash
 qsub job.scr
```

Coordinates file {{< file xyz >}}. Take a look at it (on your local machine) with visualization software such as {{< code "xcrysden" >}} to see what kind of molecule we are dealing with.
{{% expand "Expand: xyz file" %}}
{{< code-block >}}
 148
 units: A
      C                  -10.053414   -2.206790   -1.658417
      C                   -9.302986   -0.895144   -1.860284
      O                   -8.207781   -0.752638   -1.308974
      H                  -10.995492   -2.174759   -2.206494
      H                  -10.255731   -2.349212   -0.596194
      H                   -9.446853   -3.033807   -2.027068
      N                   -9.793020    0.106019   -2.626896
      C                   -9.183715    1.417441   -2.564153
      H                  -10.613638    0.026372   -3.226328
      H                   -9.703953    2.093287   -3.243002
      H                   -8.135921    1.345699   -2.855849
      H                   -9.251150    1.801765   -1.546540
      Mg                  -7.589714   -0.406247    0.589093
      C                   -5.742517   -3.346473    0.598957
      C                   -4.800687    1.285199   -0.432991
      C                   -9.267821    2.551225    1.101757
      C                  -10.244120   -2.188907    2.138434
      N                   -5.505702   -0.986494    0.214199
      C                   -5.032092   -2.255827    0.120854
      C                   -3.658251   -2.282504   -0.557964
      C                   -3.334778   -0.768793   -0.749853
      C                   -4.616448   -0.057153   -0.306840
      C                   -2.157453   -0.329206    0.109120
      C                   -3.729999   -3.115235   -1.844687
      C                   -2.766872   -2.716411   -2.960238
      C                   -1.331215   -2.744537   -2.525131
      O                   -0.791040   -3.518163   -1.775651
      O                   -0.671785   -1.759611   -3.177730
      N                   -7.068068    1.624336    0.467327
      C                   -5.959997    2.072749   -0.116268
      C                   -6.119618    3.519474   -0.405906
      C                   -7.367643    3.886460    0.028358
      C                   -8.016810    2.676853    0.594467
      C                   -5.092331    4.342863   -1.061147
      C                   -7.969215    5.193689   -0.096865
      C                   -8.764346    5.768120    0.809053
      N                   -9.475374    0.106282    1.501210
      C                   -9.944387    1.380878    1.591879
      C                  -11.236629    1.354534    2.222881
      C                  -11.533164   -0.000127    2.475692
      C                  -10.409531   -0.760031    2.032458
      C                  -12.065705    2.482921    2.565847
      O                  -11.852095    3.614620    2.163993
      C                  -12.773981   -0.530453    3.073700
      C                  -13.971292   -0.223432    2.172826
      N                   -7.915114   -2.308227    1.405664
      C                   -9.126921   -2.901012    1.866373
      C                   -8.885865   -4.386110    1.948373
      C                   -7.612956   -4.594110    1.506344
      C                   -7.041245   -3.265509    1.152599
      C                   -9.873356   -5.349056    2.414656
      C                   -6.572822   -5.562522    1.244981
      O                   -6.595619   -6.753651    1.424254
      C                   -5.368183   -4.810959    0.605175
      C                   -4.115033   -5.056118    1.400777
      O                   -3.994084   -5.204225    2.587175
      O                   -3.071063   -5.083305    0.537626
      C                   -1.772257   -5.277910    1.082827
      C                    0.753800   -1.655738   -3.032461
      C                    1.062861   -0.433850   -2.240797
      C                    1.344919   -0.409138   -0.931698
      C                    1.342960   -1.602548   -0.036439
      C                    1.703393    0.883925   -0.254470
      C                    3.229775    1.030190   -0.125447
      C                    3.578204    2.129802    0.882753
      C                    5.078428    2.494313    0.882648
      C                    5.257557    3.825237    1.627432
      C                    5.895268    1.372849    1.548820
      C                    7.408164    1.625845    1.487162
      C                    8.167514    0.444435    2.106308
      C                    9.683652    0.673563    2.265625
      C                    9.978504    1.909940    3.121920
      C                   10.398360    0.775333    0.905797
      C                   10.942484   -0.577896    0.428481
      C                   12.283515   -0.895300    1.103652
      C                   12.802866   -2.293828    0.723383
      C                   14.303660   -2.390580    1.022546
      C                   12.056040   -3.380248    1.508281
      H                   -3.971743    1.871588   -0.839072
      H                   -9.894483    3.467231    1.143013
      H                  -11.139444   -2.723143    2.491943
      H                   -2.902664   -2.745806    0.121901
      H                   -3.127452   -0.534048   -1.818442
      H                   -1.989258    0.753032    0.043950
      H                   -1.225073   -0.813801   -0.218322
      H                   -2.303939   -0.577629    1.169379
      H                   -3.543942   -4.185902   -1.581612
      H                   -4.763611   -3.094413   -2.254413
      H                   -2.880990   -3.436020   -3.813951
      H                   -3.035086   -1.724063   -3.383246
      H                   -5.426480    5.379329   -1.240079
      H                   -4.814253    3.940870   -2.049085
      H                   -4.174428    4.410148   -0.457160
      H                   -7.707449    5.719082   -1.025637
      H                   -9.043716    5.319894    1.750247
      H                   -9.186822    6.753651    0.680885
      H                  -12.935051    2.315258    3.226144
      H                  -12.713878   -1.625588    3.247884
      H                  -12.955458   -0.101758    4.088044
      H                  -14.185704    0.848513    2.128895
      H                  -13.784329   -0.571622    1.147629
      H                  -14.870493   -0.715865    2.568269
      H                   -9.798172   -6.334093    1.921673
      H                  -10.917255   -5.017070    2.290652
      H                   -9.746714   -5.558010    3.501605
      H                   -5.239174   -5.189527   -0.447800
      H                   -1.114846   -4.997929    0.244201
      H                   -1.652845   -6.330646    1.359351
      H                   -1.603686   -4.627444    1.947904
      H                    1.089674   -1.567319   -4.088044
      H                    1.172188   -2.585399   -2.602590
      H                    1.049836    0.486560   -2.828784
      H                    0.921746   -2.501009   -0.512396
      H                    0.749437   -1.428138    0.873723
      H                    2.361521   -1.853834    0.298450
      H                    1.296005    1.752865   -0.810301
      H                    1.230975    0.930527    0.748394
      H                    3.680330    0.071964    0.195869
      H                    3.670214    1.253311   -1.114605
      H                    3.279266    1.810989    1.899143
      H                    2.984923    3.037266    0.660345
      H                    5.423686    2.621092   -0.169554
      H                    6.303818    4.153543    1.610281
      H                    4.955503    3.747690    2.675806
      H                    4.666802    4.624833    1.168771
      H                    5.581171    1.260899    2.604502
      H                    5.663558    0.409177    1.057908
      H                    7.657453    2.557380    2.028025
      H                    7.731488    1.788373    0.443677
      H                    8.007397   -0.459379    1.487656
      H                    7.741534    0.221822    3.106361
      H                   10.097389   -0.215807    2.810860
      H                    9.433721    1.873203    4.072707
      H                    9.711620    2.837351    2.607095
      H                   11.049875    1.968786    3.352603
      H                    9.720700    1.204097    0.146879
      H                   11.237355    1.495880    0.987191
      H                   10.208178   -1.378522    0.628804
      H                   11.075839   -0.562083   -0.668671
      H                   13.030624   -0.126520    0.820076
      H                   12.191260   -0.815694    2.203489
      H                   12.646233   -2.454757   -0.370577
      H                   14.870493   -1.631152    0.474361
      H                   14.499359   -2.243318    2.092463
      H                   14.704426   -3.370353    0.747273
      H                   10.974607   -3.312656    1.362529
      H                   12.253579   -3.281398    2.585541
      H                   12.380065   -4.380739    1.209889
{{< /code-block >}}
{{% /expand %}}

When your job finishes, take a look at the output to see what happened and make sure it completed successfully. Then we can do time-propagation.

#  Running the time-dependent profiling  

We change the input file accordingly. Change the {{< variable "CalculationMode" >}} from {{< code gs >}} to {{< code td >}}, and add the following lines:

{{< code-block >}}
 ##### TD
                    
 T = 18
 dt = 0.003
 {{< variable "TDPropagator" >}} = aetrs
 {{< variable "TDTimeStep" >}} = dt
 
 # Profiling
 
 {{< variable "ProfilingMode" >}} = prof_memory
 {{< variable "TDMaxSteps" >}} = 30
 {{< variable "FromScratch" >}} = yes
{{< /code-block >}}

Now it is the time to do exactly the same TD run, changing the number of CPU processor cores. You have to change the XXX to powers of 2, 2^x. Start at 64 (which will be fastest) and divide by 2, in steps down to 4. (Running on 2 or 1 cores may not work.)

```bash
 #!/bin/bash
 #PBS -q regular
 #PBS -l mppwidth=XXX
 #PBS -l advres=benasque.348
 #PBS -l walltime=0:30:00
 #PBS -N testing_chl
 #PBS -V 
 
 module load octopus/4.1.2
 cd $PBS_O_WORKDIR
 aprun -n XXX octopus_mpi &> output_td_XXX
```

Different {{< code "profiling.000xxx" >}} folders will be created with each execution. We need to process them, mainly to be able to plot the information they contain. For that we can run the next script. It runs fine without any argument, but we can have more control in the files that it is going to process by using the following arguments: "analyze.sh 64 000004 2". The first argument is the biggest number of CPU processor cores that is going to be considered. The second optional argument is the number of the reference file that is going to be used. The third one is the starting number, i.e. the smallest number of CPU cores to consider.

```bash
 #!/bin/bash
 
 ## Copyright (C) 2012,2014 J. Alberdi-Rodriguez
 ##
 ## This program is free software; you can redistribute it and/or modify
 ## it under the terms of the GNU General Public License as published by
 ## the Free Software Foundation; either version 2, or (at your option)
 ## any later version.
 ##
 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.
 ##
 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, write to the Free Software
 ## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 ## 02111-1307, USA.
 ##
 ## analyze.sh
 
 # Define the biggest number of processors.
 if [ -z "$1" ]; then
     last_proc=64
 else
     last_proc=$1
 fi
 
 # Define the reference file/folder
 if [ -z "$2" ]; then
     ref=000004
 else
     ref=$2
 fi
 
 # Define the starting value
 if [ -z "$3" ]; then
     start=1
 else
     start=$3
 fi
 
 #Initialise the output file
 echo "-  " > profile_$start
 for ((num=$start;num<=$last_proc;num*=2)); do
     echo $num >> profile_$start
 done
 
 rm -f tmp1
 
 # Analyze all profiling.XXXXXX/time.000000 to get the time per subroutine
 count=0
 for function_name in $(less profiling.$ref/time.000000  |  awk '{print $1}')
 do
     if [ $count -lt 4 ]; then
 	count=$((count + 1))
     else
 	echo $function_name >> tmp1
        -iterate over power of two profilings
 	for ((num=$start;num<=$last_proc;num*=2)); do
 	    folder=`printf 'profiling.%06d\n' $num `
 	    x=$(less $folder/time.000000 | grep "^$function_name " | awk '{print $3}' )
 	    zero=_"$x"_
 	    if [ "$zero" != "__" ]; then
 		echo $x >> tmp1
 	    else
 		echo "0" >> tmp1
 	    fi
 	done
 	paste profile_$start tmp1 > tmp2
 	rm tmp1
 	cp tmp2 profile_$start
     fi
 done
 
 echo "The result is in the \"profile_$start\" file"
```

At this point we should run "analyze.sh 64 000004 2". Thus, we will create files named "profile_2". You can take a look at the following columns in the profiling data:

* TIME_STEP; the iteration time. It has a good scaling.
* COMPLETE_DATASET; the whole time of the execution. In general it decreases, it is more obvious in a real execution, where the initialization time is the same and execution one is bigger.
* SYSTEM_INIT; initialization time. We were able to stop the increasing time, and now is almost constant independently of the number of processes.
* POISSON_SOLVER; execution time for the Poisson. It is somehow constant in this case, but now with the other solvers and domain parallelization.
* RESTART_WRITE; time for writing the restart files. It depends much in the system status, more than in the number of running processes. Could be heavily decreased if it is written to the local drive.

Now we can plot it using the following script:
```bash
 #!/bin/bash
 
 ## Copyright (C) 2014 J. Alberdi-Rodriguez
 ##
 ## This program is free software; you can redistribute it and/or modify
 ## it under the terms of the GNU General Public License as published by
 ## the Free Software Foundation; either version 2, or (at your option)
 ## any later version.
 ##
 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.
 ##
 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, write to the Free Software
 ## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 ## 02111-1307, USA.
 ##
 ## plot_function.sh
 
 if [ $- -eq 0 ]; then
    function="TIME_STEP"
 else
    function=$1
 fi
 echo $function
 
 column_number=$( awk -v fun=$function '                                                                   
 { for(i=1;i<=NF;i++){                                                                                     
    if ($i == fun)                                                                                        
       {print i+1 }                                                                                       
    }                                                                                                     
 }' profile_2 )
 
 script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
 
 sed  "s/REF/$column_number/g" $script_dir/plot_ref > plot_base
 sed -i "s/FUNCTION/$function/g" plot_base
 gnuplot plot_base
```

We also need this auxiliary file:
```bash
 ## Copyright (C) 2014 J. Alberdi-Rodriguez
 ##
 ## This program is free software; you can redistribute it and/or modify
 ## it under the terms of the GNU General Public License as published by
 ## the Free Software Foundation; either version 2, or (at your option)
 ## any later version.
 ##
 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.
 ##
 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, write to the Free Software
 ## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 ## 02111-1307, USA.
 ##
 ## plot_ref
 
 set t postscript eps enhanced color solid
 set output "gnuplot.eps"
 
 set xlabel "MPI processes"
 set ylabel "t (s)"
 set logscale yx 2
 
 plot "profile_2" u 1:REF w linespoint t "FUNCTION 2^x"
```

Something else you can try is 12, 24, 48 and 96 cores, because each of the nodes has 24 cores. In this case, you would need "analyze.sh 96 000003 3" to make "profile_3", and then in the plotting script,

```bash
 plot "profile_2" u 1:REF w linespoint t "FUNCTION 2^x", "profile_3" u 1:REF w lp t "FUNCTION 3·(2^x)"
```

##  Parallelization in domains vs states  

We can divide up the work among the processors in different ways, by dividing up the points into domains for each processor, or dividing the states into groups for each processor, or a combination of both. Try out different combinations by adding to your input file

{{< code-block >}}
 {{< variable "ParStates" >}} = 2
 {{< variable "ParDomains" >}} = 12
{{< /code-block >}}

and run on 24 cores, with different numbers in the first two fields whose product is the total number of processors (''e.g.'' 6 x 4, 3 x 8, ...).

##  PFFT Poisson solver  

Another thing you can try is to compare the PFFT (parallel FFT) Poisson solver against the one we were using before (look in the output file to see which one it was). You will need to use this aprun line in your job script instead of the previous one:

```bash
 aprun -n XXX /global/homes/j/joseba/octopus/bin/octopus_mpi &> output_td_XXX
```

in order to use a different octopus compilation that uses that library, and add these lines to your input file:

{{< code-block >}}
 {{< variable "PoissonSolver" >}} = fft
 {{< variable "FFTLibrary" >}} = pfft
{{< /code-block >}}

Compare some runs against one on a similar number of processors that you did previously. How does the time for this solver compare? You can also try ordinary FFTW (not parallel) with

{{< code-block >}}
 {{< variable "FFTLibrary" >}} = fftw
{{< /code-block >}}

##  Parallelization of the ground state  

We can also try different parameters and algorithms to see their effect on the speed of the ground-state calculation, for 24, 48, or 96 processors. Look each up in the variable reference to see what they mean, and see which of the options you were using in the previous runs.
* parallelization in domains vs states (as above)
* Eigensolver = rmmdiis, plan, cg, cg_new, lobpcg.
* StatesOrthogonalization = cholesky_serial, cholesky_parallel, mgs, qr
* SubspaceDiagonalization = standard, scalapack
* linear-combination of atomic orbitals (LCAO) for initial guess: {{< variable "LCAOAlternative" >}} = yes, no. In this case, add {{< variable "MaximumIter" >}} = 0 to do just LCAO rather than the whole calculations.

