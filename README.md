## Code for Multi-state MASH

## Overview
This directory contains
- A main Python script `mash.py` 
- A script `model.py` with functions for setting up the model, handling input/output, sampling initial variables, and running a debug trajectory
- A Fortran source code in `src/` that runs the trajectories
These allow you to run MASH for various potentials, initial conditions and observables.

## How to compile
Prerequisites:
- gfortran (it also compiles and runs with ifort, but you will need to manually update the makefile)
- lapack
- numpy.f2py for Python3 

Compile by typing `make [option]` and one of the following options
- `clean` Remove compiled files and restart from clean directory
- `fast` Produces fast parallelized code
- `debug` Uses plenty of warning flags and should tell at which line the code breaks.

Also check that you have the required Python packages installed. 
One way to ensure this is to (preferably in a virtual environment) run `pip install -r requirements.txt`.

## How to run
```mash.py +[args].in```
where `[args].in` is an argument file with one line per argument (allowing commented lines). See `examples` for example input files.
There are plenty of available option flags and you can add more to suit your system.

Hints: 
- Make sure `mash.py` is executable, otherwise `chmod u+x mash.py`
- Create a symbolic link to `mash.py` in some place in your PATH, so that you can call it from anywhere.
- You may need to add the directory to your PYTHONPATH.

## What comes in
Run `mash.py -h` to see the available options. In general you should know about the following:
- `model` String specifying your model system.
- Model-specific parameters like `beta`, `Delta` etc.
- `init` Integer that specifies initial state (index starts at 0)
- `dt` Timestep
- `nt` Number of timesteps
- `ntraj` Number of trajectories
- `nucsamp` Nuclear sampling option (see "sampling" in model.py)
- `elsamp` Electronic sampling option (see "sampling" in model.py)
- `obstyp` Specify what kind of observables you want to measure (e.g. `pop` for populations)

## What comes out
Depends on what `obstyp` you specified, but for `pop` you should see a file `pop.out` which contains time in the first column and then the state populations in the following columns.

## Add a new potential
The Fortran code contains a few potentials, e.g. `linvib` and `tully`. If you want to create a different potential, add another file with the subroutines `pot` and `grad` specifying the diabatic potential and gradient. Copy the initialization routine from one of the existing potentials, and add a potential-specific init subroutine to `f2py.f90`.


***


## Problems?
Please don't hesitate contacting me if you discover any bugs. (Thanks to Eric Koessler, Rochester; Annina Lieberherr, Oxford; Bokang Huo, Berkeley)
