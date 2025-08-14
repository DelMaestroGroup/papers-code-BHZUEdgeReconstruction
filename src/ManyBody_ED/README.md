# Many-Body Exact Diagonalization

A small, fast C++ code for lattice many-body **exact diagonalization (ED)** of multi-orbital BHZ model. It entails a 3-site, 2-orbital system that supports full diagonalization at/near half filling in the **grand-canonical** ensemble.

### Features
- Eigenvalues & eigenvectors (full ED)
- Magnetic and charge observables: $S^z$, $S^2$, local densities
- Spin–spin correlations $\langle \mathbf{S}_i \cdot \mathbf{S}_j \rangle$

### Requirements
- CMake
- g++
- C++11 or newer
- LAPACK and BLAS

### Details and Compilation
The system variables are hard-coded in the main.cpp file. To make any changes one needs to modify the main.cpp. System variables includes:
- BHZ parameters: $A, B, M$ 
- Interaction parameters: $U, {J}_H $
- Number of particles: $N$

To compile the code run:
```bash
make -f Makefile
```
