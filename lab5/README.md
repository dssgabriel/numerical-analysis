# Solving the 1D stationary heat equation

Author: T. Dufaud

---

This directory contains the INCOMPLETE code corresponding to the solution
of Poisson 1D problem by direct method or iterative method.
It is organized in three directories:
- src/ 
- include/
- bin/

`src` contains the source codes, `include` contains the 
header files and `bin` contains the executables. 

The compilation and execution can be done using the Makefile.

Here are the main targets: 
- testenv
- tp2poisson1D_direct
- tp2poisson1D_iter

The commands are:
```sh
# Compile an executable bin/target 
make target

# Compile the executables corresponding to all targets
make all

# Execute ./bin/target
make run_target
```
