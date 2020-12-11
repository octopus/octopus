---
title: "BerkeleyGW"
tags: ["Tutorial", "Bulk", "GW"]
series: "Tutorial"
---


'''NOTE''': This tutorial page is set up for the [http://www.benasque.org/2018tddft Benasque TDDFT school 2018].

* Your instructors: David Strubbe and Adriel González.

* Test your knowledge with a [http://faculty.ucmerced.edu/dstrubbe/BerkeleyGW_GW_BSE_quiz.pdf quiz]! 

##  Interacting with Cori  

The [[BerkeleyGW]] tutorial is done on the [https://www.nersc.gov/users/computational-systems/cori Cori supercomputer] (at NERSC in California). There are a few key things you need to know about how to interact with the machine:
* To log in, run {{< code "ssh trainXXX@cori.nersc.gov" >}} in your terminal, substituting the actual name of your training account for {{< code "XXX" >}}.
* Be aware that since this machine is far away, running X-Windows programs will be very slow.
* For visualization with XCrySDen or other tools, installing and running on your laptop is recommended.
* We submit jobs using the SLURM queue manager, using a "reservation" for Aug 22 (benasque2018_1), 23 (benasque2018_2), and 24 (benasque2018_3). The reservations are usable from 15:00 to 19:00 each day. If you want to run outside of these times, use the debug queue. The tutorials are designed to be run using an "interactive job", in which you have a node continually available for your use and then executables launched from the command line are run immediately. Start an interactive job like this (do not submit more than one):
```text
 salloc -N 1 -q regular --reservation=benasque2018_2 -t 03:00:00 -C haswell
```
If you are using the debug queue rather than a reservation, use this:
```text
 salloc -N 1 -q debug -t 00:30:00 -C haswell
```
<!-- To use a reservation, you need to add the --reservation=<reservation name> flag into your submission. So for example for the second day you could either put the following line in your batch script:
```text
 -SBATCH --reservation=benasque2018_2
```
or add the flag on the command line when you submit your script:
```text
 sbatch --reservation=benasque2018_2 ./myscript.sl
```
It also works with interactive jobs:
```text
 salloc --reservation=benasque2018_2 <other arguments> -->
```

* Here is an example script using 32 cores. If you are using an interactive job, don't use this.

```text
 -!/bin/bash
 
 -SBATCH -J test_pulpo
 -SBATCH -N 1
 -SBATCH -C haswell
 -SBATCH -p debug
 -SBATCH -t 00:30:00
 -SBATCH --export=ALL
 
 module load octopus/8.2
 srun -n 32 octopus &> output
```

* To copy files from Cori to your local machine, in a terminal on your local machine, write {{< code "scp trainXXX@cori.nersc.gov:FULL_PATH_TO_YOUR_FILE ." >}} (filling in the username and filename) and enter your password when prompted. For very small ASCII files, you may find cut and paste more convenient.
* To see if you have jobs running, do {{< code "squeue -u trainXXX" >}}. You should not have more than one in the queue; if you do, cancel them with {{< code "scancel JOBID" >}}, filling in the number for JOBID from the output of {{< code "squeue" >}}.
* Accounts will expire on September 11. Feel free to copy the files off the machine before that to somewhere else for your future reference.

##  Documentation and resources  

* [http://www.tddft.org/programs/octopus/doc/generated/html/vars.php Octopus variable reference] 
* [http://oldsite.berkeleygw.org/releases/manual_v1.2.0.html BerkeleyGW manual for version 1.2.0]
* [http://www.benasque.org/2018tddft/talks_contr/232_BerkeleyGW_octopus_2018.pdf Lecture: Practical calculations with the GW approximation and Bethe-Salpeter equation in BerkeleyGW]
* [http://faculty.ucmerced.edu/dstrubbe/BerkeleyGW_BSE_2018.pdf Lecture: Practical BSE Calculations with BerkeleyGW + Octopus]
* [http://arxiv.org/abs/1111.4429 BerkeleyGW implementation paper] on arxiv
* More extensive [https://sites.google.com/site/berkeleygw2018/about lecture slides and other examples] from a longer tutorial devoted solely to BerkeleyGW in January 2018
* [https://drive.google.com/file/d/1mLCZrvpZrwsfsnvSBfZs0aKpKFWvSZIU/view?usp=drive_web slides on developing an interface for a new DFT code]

##  Instructions  

The first time you log in, execute these lines which will help you see color-coding for what is a link, executable, or directory:

```text
 echo 'alias ls="ls --color"' >> ~/.bashrc.ext
 . ~/.bashrc
```

Each time you log in, you should do this:

```text
 - Load modules
 module load berkeleygw/1.2
```

```text
 - Go to the scratch directory, where all runs should happen.
 cd $SCRATCH
```

To begin with the examples,

```text
 - List all examples available
 ls /project/projectdirs/mp149/Benasque2018
```

```text
 - Copy 1-boron_nitride example to your directory
 cp -R /project/projectdirs/mp149/Benasque2018/1-boron_nitride .
```

```text
 - Go to your local folder and follow instructions
 cd 1-boron_nitride
 less README
```

##  Schedule  

* Day 1
** 1-boron_nitride, GW
** 2-benzene, GW
* Day 2
** 2-benzene, Bethe-Salpeter
** 3-xct_LiCl, exciton visualization (download [http://faculty.ucmerced.edu/dstrubbe/xctLiCl.zip here] or from Cori, unpack the archive with "unzip", and follow README)
** 4-silicon, Bethe-Salpeter
* Day 3 (optional)
** Stretch goals from 2-benzene
** Any of the plane-wave BerkeleyGW examples, on Cori at {{< code "/project/projectdirs/mp149/BGW-2018/" >}}: silicon with PARATEC; boron nitride sheet, benzene, or sodium with Quantum ESPRESSO.

##  Errata  

* In {{< code "1-boron_nitride/1-mf/2.1-epsilon" >}}, the reference to {{< code "./01-calculate_wfn.qsub" >}} should be {{< code "./02-calculate_wfn.qsub" >}}, and you should run {{< code "01-get_kgrid.sh" >}} first and look at the input and output files.
* In the output of {{< code "kgrid.x" >}}, in the {{< code "KPointsGrid" >}} block, the value for a small q-shift is not quite right, and should be divided by the number of k-points in that direction. The value given in the Octopus input files is the correct one.
* In {{< code "1-boron_nitride/2-bgw/1-epsilon" >}}, the script {{< code "./01-run_epsilon.qsub" >}} should have "32" rather than "48" in the execution line, or you will get some OMP errors.
* In {{< code "1-boron_nitride/2-bgw/2-sigma" >}}, the script {{< code "./01-run_sigma.qsub" >}} should have "32" rather than "48" in the execution line, or you will get some OMP errors.

##  For historical interest  
* [2014 version of this tutorial](../BerkeleyGW (2014) )
* [2016 version of this tutorial](../BerkeleyGW (2016) )

<span class=noprint><hr>
Back to [[Tutorials]]




---------------------------------------------
