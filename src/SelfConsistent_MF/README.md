# Real-Space Self-Consistent Mean-Field Code

### Features

- Flexible with respect to boundary conditions (PBC or OBC) along both directions x and y
- Works for both single and double softwalls on either edges along-x
- Works in both canonical and grand-canonical ensemble
- Carry out Hartree approximation for multi-orbital Hubbard model + NN interaction
- Accelerated convergence with Broyden mixing
- Computes local and non-local charge and magnetic properties
- Calculates topological band structures, spin currents, edge state wave functions
 

### Requirements

- CMake
- g++
- C++11 or newer
- LAPACK and BLAS

### Details and Compilation

To compile the code do:
```bash
make -f Makefile
```
or
```bash
make
```

To run the code do:
```bash
./SelfConsistent input.inp
```

To reproduce the figures from the replace "input.inp" with the desired inputs while running the code. The relevant inputs are added in this repo. 
- For reproducing Figure-2 of the paper run "input_36x36_v1p9.inp". This will give you the spectral function data. 
-For reproducing Supplementary Figure-1 of the paper run "input_36x36_v2p0_e0p0.inp" and "input_36x36_v2p0_e0p3.inp"

Note: if the seed doesn't converge change the seed value or the number of iterations within these inputs. For canonical runs put "Ensemble=CE" and for grand-canonical runs put "Ensemble=GCE" within these inputs.
