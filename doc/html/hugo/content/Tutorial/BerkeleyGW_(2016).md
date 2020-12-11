---
title: "BerkeleyGW (2016)"
tags: ["Tutorial", "Bulk", "GW"]
series: "Tutorial"
---


'''NOTE''': This tutorial page is set up for the [http://www.benasque.org/2016tddft Benasque TDDFT school 2016]. More recent version: [[Tutorial:BerkeleyGW]]

##  Interacting with Cori  

The [[BerkeleyGW]] tutorial is done on the [https://www.nersc.gov/users/computational-systems/cori Cori supercomputer] (at NERSC in California). There are a few key things you need to know about how to interact with the machine:
* To log in, run {{< code "ssh trainXXX@cori.nersc.gov" >}} in your terminal, substituting the actual name of your training account for {{< code "XXX" >}}.
* Be aware that since this machine is far away, you should not try running X-Windows programs.
* We will work with "debug" jobs which can submitted via the provided scripts in each directory. We are using the SLURM queue manager on Cori.
* To copy files from Cori to your local machine, in a terminal on your local machine, write {{< code "scp trainXXX@cori.nersc.gov:FULL_PATH_TO_YOUR_FILE ." >}} (filling in the username and filename) and enter your password when prompted. For very small ASCII files, you may find cut and paste more convenient.
* To see if you have jobs running, do {{< code "squeue -u trainXXX" >}}. You should not have more than one in the queue; if you do, cancel them with {{< code "scancel JOBID" >}}, filling in the number for JOBID from the output of {{< code "squeue" >}}.
* Accounts will expire on October 7.
* Cori is actually down for an upgrade since Sept 19, but the files we were using are available if you login to {{< code "edison.nersc.gov" >}} and go to {{< code "$CSCRATCH" >}}.

##  Documentation and resources  

* [http://www.tddft.org/programs/octopus/doc/generated/html/vars.php Octopus variable reference] 
* [http://www.berkeleygw.org/releases/manual_v1.2.0.html BerkeleyGW manual]
* [http://www.benasque.org/2016tddft/talks_contr/154_BerkeleyGW_octopus.pptx.pdf Intro slides from first day]
* [http://www.benasque.org/2016tddft/talks_contr/164_128_Felipe_BSE_Presentation.pdf Intro slides from second day]
* [http://arxiv.org/abs/1111.4429 BerkeleyGW implementation paper] on arxiv
* More extensive [https://sites.google.com/site/berkeleygw2018/about lecture slides] from a longer tutorial devoted solely to BerkeleyGW in January 2018

##  Instructions  

The first time you log in, execute these lines which will help you see color-coding for what is a link, executable, or directory:

```text
 echo 'alias ls="ls --color"' >> ~/.bashrc.ext
 . ~/.bashrc
```

Each time you log in, you should do this:

```text
 - Load modules
 module load berkeleygw/1.2 octopus/6.0
```

```text
 - Go to the scratch directory, where all runs should happen.
 cd $SCRATCH
```

To begin with the examples,

```text
 - List all examples available
 ls /project/projectdirs/mp149/Benasque2016
```

```text
 - Copy 1-silicon example to your directory
 cp -R /project/projectdirs/mp149/Benasque2016/1-boron_nitride .
```

```text
 - Go to your local folder and follow instructions
 cd 1-boron_nitride
 less README
```

##  Schedule  

* Day 1
** 1-boron_nitride, GW
* Day 2
** 2-benzene, GW and [http://en.wikipedia.org/wiki/Hans_Bethe Bethe]-[http://en.wikipedia.org/wiki/Edwin_Ernest_Salpeter Salpeter]
** 3-xct_LiCl, exciton visualization (copy to your local machine according to blackboard or from Cori, unpack the archive, and follow README)
** 4-silicon, Bethe-Salpeter

##  For historical interest  
* [2014 version of this tutorial](../BerkeleyGW (2014) )

<span class=noprint><hr>
Back to [[Tutorials]]




---------------------------------------------
