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

An efficient julia framework for Hierarchical Equations of Motion (HEOM) in open quantum systems  

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
The documentation and examples can be found in :
- [**STABLE (not ready yet)**](https://ncku-qfort.github.io/Heom.jl/stable) : most recently tagged version.
- [**DEVELOP**](https://ncku-qfort.github.io/Heom.jl/dev/) : in-development version.

## Cite `Heom.jl`
If you like `Heom.jl`, we would appreciate it if you starred the repository in order to help us increase its visibility. Furthermore, if you find the framework useful in your research, we would be grateful if you could cite our publication ([`Update when we published the paper`](https://unknown)  ) using the following bibtex entry:
```bib
@article{Heom2023
  title={},
  author={},
  journal={},
  year={}
}
```

## License
`Heom.jl` is released under the [Apache 2 license](./LICENSE.md).