---
title: "BerkeleyGW (2014)"
tags: ["Tutorial", "Bulk", "GW"]
#series: "Tutorial"
hidden: True
---


'''NOTE''': This tutorial page was set up for the Benasque TDDFT school 2014. The specific references to the supercomputer used at that time will have to be adapted for others to use this tutorial. More recent version: [[Tutorial:BerkeleyGW]]

##  Interacting with Hopper  

The [[BerkeleyGW]] tutorial is done on the [https://www.nersc.gov/users/computational-systems/hopper Hopper supercomputer] (at NERSC in California). There are a few key things you need to know about how to interact with the machine:
* To log in, run {{< code "ssh trainX@hopper.nersc.gov" >}} in your terminal, substituting the actual name of your training account.
* Be aware that since this machine is far away, you should not try running X-Windows programs!
* You submit jobs by the {{< code "qsub" >}} command, ''e.g.'' {{< code "qsub job.scr" >}}, which will put them in the queue for execution when there is free space.
* You can see what jobs you currently have in the queue by executing {{< code "qstat -u $USER" >}}, so you can see when your job finishes. A status code will be shown: Q = waiting in the queue, R = running, C = complete.
* You can cancel a job by {{< code "qdel job" >}}, where {{< code "job" >}} is the job id as written by {{< code "qstat" >}}.
* The job script (''e.g.'' {{< code "01-calculate_scf.qsub" >}}) specifies parameters to the PBS/Torque queuing system about how many cores to use, what commands to run, etc.
* To copy files from Hopper to your local machine, in a terminal on your local machine, write {{< code "scp trainX@hopper.nersc.gov:FULL_PATH_TO_YOUR_FILE ." >}} (filling in the username and filename) and enter your password when prompted. For very small ASCII files, you may find cut and paste more convenient.

##  Getting started in the tutorial  

###  Day 1  
To obtain the files for the boron nitride and benzene examples for the first day of the tutorial:

```text
 cd $SCRATCH
 /project/projectdirs/m1694/BGW-tddft/copy_day_1.sh
```

* In each case, enter your copy of the directory, and look at {{< code "README" >}} and follow instructions given there.
* Start by running {{< code "2-benzene/1-mf/1-scf" >}} and then {{< code "2-benzene/1-mf/2-wfn" >}}. This will take a little while, so while this runs, do the BN example.

###  Day 2  

To obtain the files for the silicon and benzene examples for the second day of the tutorial:

```text
 cd $SCRATCH
 /project/projectdirs/m1694/BGW-tddft/copy_day_2.sh
```

The solution for the benzene example is available at

```text
  /project/projectdirs/m1694/BGW-tddft/2-benzene_run
```

You can copy the necessary files (WFNs, bse*mat, eps*mat, eqp*) from there.

Other instructions:

* We will work on the following directories: 2-benzene and 3-silicon. We will not work on the 1-boron_nitride example!
* Start with the example 3-silicon.
* In each case, enter your copy of the directory, and look at {{< code "README" >}} and follow instructions given there.
* There is additional example for XCrySDen, which is available in the shared folder on imac01 (see instructions on blackboard). If for some reason you are not able to copy it from there, it can also be downloaded [http://civet.berkeley.edu/~jornada/files/xct_LiCl.zip here]. Note: to use XCrySDen on the iMacs, run {{< code "/sw/bin/xcrysden" >}}.

The examples are available for download here: [http://web.mit.edu/~dstrubbe/www/2-benzene.tar.gz 2-benzene.tar.gz], [http://web.mit.edu/~dstrubbe/www/3-silicon.tar.gz 3-silicon.tar.gz].

##  General workflow  

* Finish all basic goals from both the boron nitride and benzene examples before starting any stretch goal.

##  Documentation and resources  

* [http://www.tddft.org/programs/octopus/doc/generated/html/vars.php?page=sections Octopus variable reference] for the pre-release development version. 
* [http://www.berkeleygw.org/releases/manual_v1.0.6.html BerkeleyGW manual].
* [http://benasque.org/2014tddft/talks_contr/115_BerkeleyGW_octopus.pptx.pdf Intro slides from first day]
* [http://arxiv.org/abs/1111.4429 BerkeleyGW implementation paper] on arxiv.
* More extensive [http://www.nersc.gov/users/training/nersc-training-events/berkeleygw2013 lecture slides] from a longer tutorial devoted solely to BerkeleyGW in November 2012.
* The [http://benasque.org/2014tddft/talks_contr/128_Felipe_BSE_Presentation.pdf slides for the second day of tutorial].
Note that we are using the pre-release development version of Octopus in this tutorial, rather than the current release 4.1.2 which lacks full support for BerkeleyGW output. There are some small differences in output from 4.1.2.

<span class=noprint><hr>
Back to [[Tutorials]]




---------------------------------------------
