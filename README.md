![Fancy logo](./docs/src/assets/logo-dark.png#gh-dark-mode-only)
![Fancy logo](./docs/src/assets/logo.png#gh-light-mode-only)

# HierarchicalEOM.jl

| **Release**       | [![Release][release-img]][release-url] [![License][license-img]][license-url] [![Cite][cite-img]][cite-url] |
|:-----------------:|:-------------|
| **Runtests**      | [![Runtests][runtests-img]][runtests-url] [![Coverage][codecov-img]][codecov-url] |
| **Code Quality**  | [![Code Quality][code-quality-img]][code-quality-url] [![Aqua QA][aqua-img]][aqua-url] [![JET][jet-img]][jet-url] |
| **Documentation** | [![Doc-Stable][docs-stable-img]][docs-stable-url] [![Doc-Dev][docs-develop-img]][docs-develop-url] |

[release-img]: https://img.shields.io/github/release/qutip/HierarchicalEOM.jl.svg
[release-url]: https://github.com/qutip/HierarchicalEOM.jl/releases

[license-img]: https://img.shields.io/badge/license-New%20BSD-blue.svg
[license-url]: https://opensource.org/licenses/BSD-3-Clause

[cite-img]: https://img.shields.io/badge/cite-Commun._Phys._6%2C_313_(2023)-blue
[cite-url]: https://doi.org/10.1038/s42005-023-01427-2

[runtests-img]: https://github.com/qutip/HierarchicalEOM.jl/actions/workflows/Runtests.yml/badge.svg
[runtests-url]: https://github.com/qutip/HierarchicalEOM.jl/actions/workflows/Runtests.yml

[codecov-img]: https://codecov.io/gh/qutip/HierarchicalEOM.jl/graph/badge.svg?token=ICFVVNuLHW
[codecov-url]: https://codecov.io/gh/qutip/HierarchicalEOM.jl

[code-quality-img]: https://github.com/qutip/HierarchicalEOM.jl/actions/workflows/Code-Quality.yml/badge.svg 
[code-quality-url]: https://github.com/qutip/HierarchicalEOM.jl/actions/workflows/Code-Quality.yml

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[jet-img]: https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a
[jet-url]: https://github.com/aviatesk/JET.jl

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://qutip.org/HierarchicalEOM.jl/stable/
[docs-develop-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-develop-url]: https://qutip.org/HierarchicalEOM.jl/dev/

`HierarchicalEOM.jl` is a numerical framework written in [`Julia`](https://julialang.org/). It provides a user-friendly and efficient tool based on hierarchical equations of motion (HEOM) approach to simulate complex open quantum systems, including non-Markovian effects due to non-perturbative interaction with one (or multiple) environment(s). It is built upon [`QuantumToolbox.jl`](https://github.com/qutip/QuantumToolbox.jl).

![](./docs/src/assets/heom_ecosystem.jpeg)

## Installation

> [!NOTE]
> `HierarchicalEOM.jl` requires `Julia 1.10+`.

To install `HierarchicalEOM.jl`, run the following commands inside Julia's interactive session (also known as REPL):
```julia
using Pkg
Pkg.add("HierarchicalEOM")
```
Alternatively, this can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/) by pressing the key `]` in the REPL to use the package mode, and then type the following command:
```julia-REPL
(1.10) pkg> add HierarchicalEOM
```
More information about `Julia`'s package manager can be found at [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/).  

To load the package and check the version information, use either `HierarchicalEOM.versioninfo()` or `HierarchicalEOM.about()`, namely
```julia
using HierarchicalEOM
HierarchicalEOM.versioninfo()
HierarchicalEOM.about()
```

## Documentation
The documentation can be found in :
- [**STABLE**](https://qutip.org/HierarchicalEOM.jl/stable) : most recently tagged version.
- [**DEVELOP**](https://qutip.org/HierarchicalEOM.jl/dev/) : in-development version.

## Cite `HierarchicalEOM.jl`
If you like `HierarchicalEOM.jl`, we would appreciate it if you starred the repository in order to help us increase its visibility. Furthermore, if you find the framework useful in your research, we would be grateful if you could cite our publication [ [Communications Physics 6, 313 (2023)](https://doi.org/10.1038/s42005-023-01427-2)  ] using the following bibtex entry:
```bib
@article{HierarchicalEOM-jl2023,
  doi = {10.1038/s42005-023-01427-2},
  url = {https://doi.org/10.1038/s42005-023-01427-2},
  year = {2023},
  month = {Oct},
  publisher = {Nature Portfolio},
  volume = {6},
  number = {1},
  pages = {313},
  author = {Huang, Yi-Te and Kuo, Po-Chen and Lambert, Neill and Cirio, Mauro and Cross, Simon and Yang, Shen-Liang and Nori, Franco and Chen, Yueh-Nan},
  title = {An efficient {J}ulia framework for hierarchical equations of motion in open quantum systems},
  journal = {Communications Physics}
}
```

## License
`HierarchicalEOM.jl` is released under the [BSD 3-Clause License](./LICENSE.md).

## Contributing to HierarchicalEOM.jl

You are most welcome to contribute to `HierarchicalEOM.jl` development by forking this repository and sending pull requests (PRs), or filing bug reports at the issues page. You can also help out with users' questions, or discuss proposed changes in the [QuTiP discussion group](https://groups.google.com/g/qutip).

For more information about contribution, including technical advice, please see the [Contributing to QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/stable/resources/contributing) section of the `QuantumToolbox.jl` documentation.

