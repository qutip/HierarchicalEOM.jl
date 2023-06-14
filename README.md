![Fancy logo](./docs/src/assets/logo-dark.png#gh-dark-mode-only)
![Fancy logo](./docs/src/assets/logo.png#gh-light-mode-only)

# Heom.jl
<!--
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ncku-qfort.github.io/Heom.jl/dev/)
-->
[![](https://img.shields.io/github/release/NCKU-QFort/Heom.jl.svg)](https://github.com/NCKU-QFort/Heom.jl/releases)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ncku-qfort.github.io/Heom.jl/dev/)  
[![Runtests](https://github.com/NCKU-QFort/Heom.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/NCKU-QFort/Heom.jl/actions/workflows/Runtests.yml)
[![codecov](https://codecov.io/gh/NCKU-QFort/Heom.jl/branch/main/graph/badge.svg?token=237Z7F7OOV)](https://codecov.io/gh/NCKU-QFort/Heom.jl)

`Heom.jl` is a numerical framework written in [`Julia`](https://julialang.org/). It provides a user-friendly and efficient tool based on hierarchical equations of motion (HEOM) approach to simulate complex open quantum systems, including non-Markovian effects due to non-perturbative interaction with one (or multiple) environment(s). It is inspired by the [Quantum Toolbox in Python (QuTiP)](https://qutip.org).

## Installation
To install `Heom.jl`, run the following commands inside Julia's interactive session (also known as REPL):
```julia
using Pkg
Pkg.add("Heom")
```
Alternatively, this can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/) by pressing the key `]` in the REPL to use the package mode, and then type the following command:
```julia-REPL
(1.8) pkg> add Heom
```
More information about `Julia`'s package manager can be found at [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/).  
`Heom` requires Julia 1.8 or higher. Installing it on an older version of Julia will result in many errors.

To load the package and check the version information, use the command:
```julia
julia> using Heom
julia> Heom.versioninfo()
```

## Documentation
The documentation can be found in :
- [**STABLE**](https://ncku-qfort.github.io/Heom.jl/stable) : most recently tagged version.
- [**DEVELOP**](https://ncku-qfort.github.io/Heom.jl/dev/) : in-development version.

## Cite `Heom.jl`
If you like `Heom.jl`, we would appreciate it if you starred the repository in order to help us increase its visibility. Furthermore, if you find the framework useful in your research, we would be grateful if you could cite our publication ( [`arXiv:2306.07522`](https://doi.org/10.48550/arXiv.2306.07522)  ) using the following bibtex entry:
```bib
@article{Heom-jl2023,
  title={{Heom.jl}: {A}n efficient {J}ulia framework for hierarchical equations of motion in open quantum systems},
  author={Huang, Yi-Te and Kuo, Po-Chen and Lambert, Neill and Cirio, Mauro and Cross, Simon and Yang, Shen-Liang and Nori, Franco and Chen, Yueh-Nan},
  journal={arXiv preprint arXiv:2306.07522},
  year={2023}
}
```

## License
`Heom.jl` is released under the [Apache 2 license](./LICENSE.md).