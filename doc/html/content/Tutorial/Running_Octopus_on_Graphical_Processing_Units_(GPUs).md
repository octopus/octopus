---
title: "Running Octopus on Graphical Processing Units (GPUs)"
tags: ["Tutorial", "Expert"]
#series: "Tutorial"
---


{{< notice warning >}}
This chapter needs to be updated, as the GPU implementation has changed since it was written!!!
{{< /notice >}}


Recent versions of Octopus support execution on graphical processing
units (GPUs). In this tutorial we explain how the GPU support works
and how it can be used to speed-up calculations.

{{< octopus >}} is in the process of being ported to machines with graphical processing units (GPUs). Not all operations are supported yet, but the number of supported features is growing with time. 

{{< notice note >}}
Note that the code might fall back to CPU operation for unsupported features. 
{{< /notice >}}

###  Consideration before using a GPU  

####   Calculations that will be accelerated by GPUs  

A GPU is a massively parallel processor, to work efficiently it need
large amounts of data to process. This means that using GPUs to
simulate small systems, like molecules of less than 20 atoms or low
dimensionality systems, will probably will be slower than using the CPU.

Not all the calculations you can do with Octopus will be effectively
accelerated by GPUs. Essentially ground state calculations with the
rmmdiis eigensolver (and calculations based on ground state calculations, like geometry optimization) and time propagation with the etrs and aetrs
propagators are the simulations that will use the GPU more
effectively. For other types of calculations you might see some
improvements in perfomance. Moreover, some Octopus features do not
support GPUs at all. Currently these are:

* HGH or relativistic pseudopotentials.
* Non-collinear spin.
* Curvilinear coordinates.

* Periodic systems are partially supported, but there are limitations which might force the code to run on the CPU.

In these cases Octopus will be silently executed on the CPU.

####  Supported GPUs  

Octopus GPU support is based on the OpenCL framework, so it is vendor
independent, currently it has been tested on GPUs from Nvidia and AMD
(ex-ATI). Recently, CUDA support is added using a CUDA wrapping layer.

{{< notice warning >}}
Actually, the non-CUDA version is currently untested, as our development systems are all NVidia CUDA based.
{{< /notice >}}

In order to run Octopus, the GPU must have double precision
support. For AMD/ATI cards this is only available in high-end models:
Radeon 58xx, 69xx, 78xx, 79xx, and Firepro/Firestream cards. All
recent Nvidia cards support double precision (4xx and newer), low- and
middle-end cards however have very limited processing power and are
probably slower running Octopus than a CPU.

The other important factor is GPU memory, you need a GPU with at least
1.5 GiB of RAM to run Octopus, but 3 GiB or more are recommended. 


####  Activating GPU support  

Octopus has GPU support since version 4.0; to activate it you need a version compiled with OpenCL support enabled (see the {{< manual "Building_Octopus_with_GPU_support" "manual" >}} for details). If your version of Octopus was compiled with OpenCL support you should see that 'opencl' among the configuration options:

```text
                            Running octopus
 
 Version                : superciliosus
 Revision               :
 Build time             : Thu Feb 14 22:19:10 EST 2013
 Configuration options  : max-dim=3 openmp '''opencl''' sse2 avx
 Optional libraries     : netcdf '''clamdfft clamdblas'''
```

Ideally your Octopus version should also be compiled with the optional libraries 'clamdfft' and  'clamdblas', that are part of the [http://developer.amd.com/tools/heterogeneous-computing/amd-accelerated-parallel-processing-math-libraries/ AMD APPML] package (these libraries work on hardware from other vendors too, including Nvidia).

###  Starting the tutorial: running without OpenCL  

If you have a version of Octopus compiled with OpenCL support, by setting the {{< variable "DisableAccel" >}} to {{< code "yes" >}} you can tell Octopus not to use OpenCL.

###  Selecting the OpenCL device  

OpenCL is designed to work with multiple devices from different vendors. Most likely you will have a single OpenCL device, but in some cases you need to select the OpenCL device that Octopus should use.

In OpenCL a 'platform' identifies an OpenCL implementation, you can have multiple platforms installed in a system. Octopus will list the available platforms at the beginning of the execution, for example:

```text
 Info: Available CL platforms: 2
     * Platform 0 : NVIDIA CUDA (OpenCL 1.1 CUDA 4.2.1)
       Platform 1 : AMD Accelerated Parallel Processing (OpenCL 1.2 AMD-APP (1084.4))
```

By default Octopus will use the first platform. You can select a different platform with the 
{{< max-version 9 >}}{{< variable "OpenCLPlatform" >}}{{< /max-version >}}{{< min-version 10 >}}{{< variable "AccelPlatform" >}}{{< /min-version >}} 
variable. You can select the variable by name, for example:

{{< min-version "10" >}}
{{< code-block >}}
{{< variable "AccelPlatform" >}} = amd
{{< /code-block >}}
{{< /min-version  >}}

{{< max-version "9" >}}
{{< code-block >}}
{{< variable "OpenCLPlatform" >}} = amd
{{< /code-block >}}
{{< /max-version  >}}

or you can use numeric values that select the platform by its index in the list:

{{< min-version "10" >}}
{{< code-block >}}
{{< variable "AccelPlatform" >}} = 1
{{< /code-block >}}
{{< /min-version >}}

{{< max-version "9" >}}
{{< code-block >}}
{{< variable "OpenCLPlatform" >}} = 1
{{< /code-block >}}
{{< /max-version >}}

The first form should be preferred as it is independent of the order of the platforms.

Each OpenCL platform can contain several devices, the available devices on the selected platform are printed by Octopus, like this:

```text
 Info: Available CL devices: 2
       Device 0 : Tahiti
       Device 1 : Intel(R) Core(TM) i7-3820 CPU @ 3.60GHz
```

Each OpenCL device has a certain type that can be 'cpu', 'gpu', or 'accelerator'. By default Octopus tries to use the first device of type 'gpu', if this fails Octopus tries to use the default device for the system, as defined by OpenCL implementation. To select a different device you can use the 
{{< max-version 9 >}}{{< variable "OpenCLDevice" >}}{{< /max-version >}}{{< min-version 10 >}}{{< variable "AccelDevice" >}}{{< /min-version >}} variable, you can select the device by its index or by its type ('cpu', 'gpu', or 'accelerator').

{{< notice note >}}
the Octopus OpenCL implementation is currently designed to execute efficiently on GPUs. Octopus will run on CPUs using OpenCL, but it will be slower than using the CPU directly.
{{< /notice >}}

After Octopus has selected an OpenCL device it will print some information about it, for example:

```text
 Selected CL device:
      Device vendor          : Advanced Micro Devices, Inc.
      Device name            : Tahiti
      Driver version         : 1084.4 (VM)
      Compute units          : 32
      Clock frequency        : 925 GHz
      Device memory          : 2048 MiB
      Max alloc size         : 512 MiB
      Device cache           : 16 KiB
      Local memory           : 32 KiB
      Constant memory        : 64 KiB
      Max. workgroup size    : 256
      Extension cl_khr_fp64  : yes
      Extension cl_amd_fp64  : yes
```



###  The block-size  

To obtain good peformance on the GPU (and CPUs), Octopus uses blocks of states as the basic data object. The size of these blocks is controlled by the {{< variable "StatesBlockSize" >}} variable. For GPUs the value must be a power of 2 (the default is 32) and for optimal performance must be as large as possible, the catch is that for large values you might run out of GPU memory.

###  StatesPack  

###  The Poisson solver  

If you want to used a GPU-accelerated Poisson solver, you need to use the 'fft' solver and to have the 'clamdfft' library support compiled. To select the solver, just add

```text
 {{< variable "PoissonSolver" >}} = fft
 {{< variable "FFTLibrary" >}} = clfft
```

to your input file.


<span class=noprint><hr>
Back to [[Tutorials]]



---------------------------------------------
