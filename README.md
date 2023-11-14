## Code for Multi-state MASH

## Overview
This directory contains
- A main Python script `mash.py` 
- A script `utils.py` with functions for handling input/output, sampling of initial variables, and debugging
- A Fortran source code in `src/` that runs the trajectories
These allow you to run MASH for various potentials, initial conditions and observables.

## How to compile
Prerequisites:
- gfortran
- f2py for Python3
- lapack

Compile by typing `make [option]` and one of the following options
- `clean` Remove compiled files and restart from clean directory
- `fast` Produces fast parallelized code
- `debug` For more serious debugging i.e. using `valgrind`

## How to run
```mash.py +[args].in```
where `[args].in` is an argument file with one line per argument (allowing commented lines). See `examples` for example input files.
There are plenty of available option flags and you can add more to suit your system.

Hints: 
- Make sure `mash.py` is executable, otherwise `chmod u+x mash.py`
- create a symbolic link to `mash.py` in some place in your PATH, so that you can call it from anywhere.
- You may need to add the present directory to your PYTHONPATH.

## What comes in
Run `mash.py -h` to see the available options. In general you should know about the following:
- `model` string specifying your model system.
- model-specific parameters like `beta`, `Delta` etc.
- `init` Integer that specifies initial state (index starts at 0)
- `dt` Timestep
- `nt` Number of timesteps
- `ntraj` Number of trajectories
- `nucsamp` Nuclear sampling option
- `obstyp` Specify what kind of observables you want to measure (e.g. `pop` for populations)

## What comes out
Depends on what `obstyp` you specified, but for `pop` you should see a file `pop.out` which contains time in the first column and then the state populations in the following columns.

## Add a new potential
The Fortran code contains a few potentials, e.g. `linvib` and `tully`. If you want to create a different potential, add another file with the subroutines `pot` and `grad` specifying the diabatic potential and gradient. Copy the initialization routine from one of the existing potentials, and (if necessary) add a potential-specific init subroutine in the beginning of `f2py.f90`.
If the potential has a particularly simple form that allows for faster computation of the adiabatic gradient, it is also possible to override the generic adiabatic gradient routines - see `frexc.f90` for an example using the Frenkel-exciton model.


***

