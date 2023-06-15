# HierarchicalEOM.jl: An efficient Julia framework for Hierarchical Equations of Motion (HEOM) in open quantum systems

`HierarchicalEOM.jl` is a numerical framework written in [`Julia`](https://julialang.org/). It provides a user-friendly and efficient tool based on hierarchical equations of motion (HEOM) approach to simulate complex open quantum systems, including non-Markovian effects due to non-perturbative interaction with one (or multiple) environment(s). It is inspired by the [Quantum Toolbox in Python (QuTiP)](https://qutip.org).

While integrating many of the features present in other open-source HEOM packages, `HierarchicalEOM.jl` also includes new functionalities, such as the construction of even- and odd-parity [HEOM Liouvillian superoperator (HEOMLS) matrices](@ref doc-HEOMLS-Matrix), the estimation of [importance values](@ref doc-Importance-Value-and-Threshold) for all [auxiliary density operators (ADOs)](@ref doc-ADOs), and the calculation of [spectra](@ref doc-Spectrum) for both bosonic and fermionic systems. 

By wrapping some functions from other [`Julia`](https://julialang.org/) packages ([`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/), [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/) and [`fastExpm.jl`](https://github.com/fmentink/FastExpm.jl)), `HierarchicalEOM.jl` collects different methods and could further optimize the computation for the [stationary state](@ref doc-Stationary-State), and the [time evolution](@ref doc-Time-Evolution) of all [ADOs](@ref doc-ADOs). The required handling of the [ADOs](@ref doc-ADOs) multi-indexes is achieved through a user-friendly interface called [Hierarchy Dictionary](@ref doc-Hierarchy-Dictionary).

![HEOM Ecosystem](assets/heom_ecosystem.jpeg)

We believe that `HierarchicalEOM.jl` will be a valuable tool for researchers working in different fields such as quantum biology, quantum optics, quantum thermodynamics, quantum information, quantum transport, and condensed matter physics.

If you like `HierarchicalEOM.jl` and find the framework useful in your research, we would be grateful if you could cite our publication ( [`arXiv:2306.07522`](https://doi.org/10.48550/arXiv.2306.07522)  ) using the bibtex entry [here](@ref doc-Cite).