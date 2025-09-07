[![Paper](https://img.shields.io/badge/paper-arXiv%3A2508.10726-B31B1B.svg)](https://arxiv.org/abs/2508.10726)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16875878.svg)](https://zenodo.org/badge/latestdoi/16875878)


# Edge Reconstruction in a Quantum Spin Hall Insulator

Rahul Soni, Matthias Thamm, Gonzalo Alvarez, Bernd Rosenow, and Adrian Del Maestro

[arXiv:2508.10726](https://arxiv.org/abs/2508.10726)

### Abstract
We study interaction-driven edge reconstruction in a quantum spin Hall insulator described by the BHZ model with Kanamori–Hubbard interactions using real-space density matrix renormalization group method in both the grand-canonical and canonical ensembles. 
For a two-dimensional cylinder with one smooth edge, we identify discrete particle-number transitions that lead to spin-polarized edge states stabilized by an emergent ferromagnetic exchange interaction. The reconstruction is orbital-selective, occurring predominantly in the $s$-orbital channel. Our results reveal a fully microscopic mechanism for emergent spin polarization at the edge that could compromise the topological protection of helical edge states by time reversal symmetry.

### Description
This repository includes links, code, scripts, and data to generate the figures in a paper.

### Requirements
The data in this project is generated using three different methods: Exact Diagonalization, Mean-Field and DMRG. Processed data is included in the [data](https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/tree/main/data) directory.


1. The spectral data for the interacting BHZ model was generated via self-consistent real-space mean-field calculations. The code can be found [here](https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/tree/main/codes/SelfConsistent_MF). Detail instructions are provided in this repo regarding compilations, executions and more.
2. The real space charge, magnetic, exact ground-state data for interacting BHZ ladders was generated via DMRG for both canonical and grand-canonical ensemble. For the large DMRG calculations we considered two different codes (benchmarked against each other).

   (a) The in-house DMRG++ software developed by G.A. The documentation for the same is provided [here](https://github.com/g1257/dmrgpp), for compilation follow the steps below:

      **DMRG++ Code and Compilation**
        
      ```bash
      git clone https://code.ornl.gov/gonzalo_3/PsimagLite.git
      cd PsimagLite/lib
      perl configure.pl
      make
      cd ../../
      git clone https://code.ornl.gov/gonzalo_3/dmrgpp.git
      cd dmrgpp/src
      perl configure.pl
      make
      ```
        
      Dependencies include the [BOOST](https://www.boost.org/), [HDF5](https://docs.hdfgroup.org/archive/support/HDF5/doc1.8/cpplus_RM/index.html) and [OpenBLAS](https://www.openblas.net/) libraries
      
      **Running the Code**
      
      This will generate `dmrg` and `observe` executables. Run the dmrg executable first to save the ground state and then use the observe executable to evaluate all the necessary observables. 


    (b) the ITensors DMRG code available [here](https://github.com/DelMaestroGroup/BHZ_DMRG_Julia/tree/main) using the [ITensor.jl](https://docs.itensor.org/ITensors/stable/) package. To run the code, clone the repository
      ```bash
      git clone https://github.com/DelMaestroGroup/BHZ_DMRG_Julia.git
      ```
      Running the code requires Julia 1.11.1 or higher:
      ```bash
      julia ./BHZitensorsDMRG/run_bhz.jl ARGS KWARGS
      ```
      For details, see [README](https://github.com/DelMaestroGroup/BHZ_DMRG_Julia/blob/main/README.md) of the repository.

       



4. For small system we considered the many-body exact diagonalization code. The code can be found [here](https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/tree/main/codes/ManyBody_ED). Detail instructions are provided in this repo regarding compilations, executions and more.

### Support
The creation of these materials was supported in part by the Department of Energy under Award No. [DE-SC0022311]([https://www.nsf.gov/awardsearch/simpleSearchResult?queryText=delmaestro](https://pamspublic.science.energy.gov/WebPAMSExternal/Interface/Common/ViewPublicAbstract.aspx?rv=31bd2b59-7a7a-424c-83cc-fad4b3df485f&rtc=24&PRoleId=10)).

<img width="400px" src="https://science.osti.gov/assets/img/doe-logos/logo.png">


### Figures
#### Figure 1: Lattice Geometry and Potential Profile
<img src="https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/blob/main/figures/Confining_Lattice_new.png" width="400px">

#### Figure 2: Mean-field Orbital and Spin Resolved Spectral Function
<img src="https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/blob/main/figures/OS_resolved_Akyw_for_Vo_9p5_36x36_SSW.png" width="600px">

#### Figure 3: Total particles versus Confining potential
<img src="https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/blob/main/figures/Tot_N_vs_V0_for_12x3_GCE.png" width="400px">

#### Figure 4: GSE for different Sz sectors of $N_0$, $N_0+1$ and $N_0+2$ particle sectors
<img src="https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/blob/main/figures/GSE_vs_Sz_for_diff_Vo_12x3_triple.png" width="500px">

#### Figure 5: Change in particle density between different particle sectors
<img src="https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/blob/main/figures/del_ni_12x3_w4.png" width="400px">

#### Figure 6: Magnetization profiles for the ground-states
<img src="https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/blob/main/figures/szi_12x3_w4.png" width="400px">

#### Supplementary Figure 1: Spectral comparison for different helical potential
<img src="https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/blob/main/figures/Akyw_for_diff_eps_Vo_10p0_36x36_SSW.png" width="600px">

#### Supplementary Figure 2: Finite size scaling of the exchange coupling
<img src="https://github.com/DelMaestroGroup/papers-code-BHZUEdgeReconstruction/blob/main/figures/FS_exchange_vs_1byLx_Nx3_DMRG.png" width="400px">


This figure is relesed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) and can be freely copied, redistributed and remixed.

