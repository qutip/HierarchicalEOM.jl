# Misc.

```@docs
Heom.versioninfo()
```
The outputs will be something like the following:
```
julia> Heom.versioninfo()

                                   __
                                  /  \
 __     __                     __ \__/ __
|  |   |  |                   /  \    /  \
|  |   |  | ______   ______   \__/_  _\__/
|  |___|  |/  __  \ /  __  \ / '   \/     \
|   ___   |  |__)  |  /  \  |    _     _   |
|  |   |  |   ____/| (    ) |   / \   / \  |
|  |   |  |  |____ |  \__/  |  |   | |   | |
|__|   |__|\______) \______/|__|   |_|   |_|

Julia framework for Hierarchical Equations of Motion
≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
Copyright © NCKU-QFORT 2022 and later.
Lead  developer : Yi-Te Huang
Other developers:
    Simon Cross, Neill Lambert, Po-Chen Kuo and Shen-Liang Yang

Package information:
====================================
Heom            Version: 0.1.0
JLD2            Version: 0.4.31
LinearSolve     Version: 1.42.0
OrdinaryDiffEq  Version: 6.52.0
ProgressMeter   Version: 1.7.2
PrecompileTools Version: 1.1.2

System information:
====================================
Julia Version: 1.8.0
OS       : macOS (x86_64-apple-darwin21.4.0)
CPU      : 12 × Intel(R) Core(TM) i7-8700B CPU @ 3.20GHz
Memory   : 16.0 GB
WORD_SIZE: 64
LIBM     : libopenlibm
LLVM     : libLLVM-13.0.1 (ORCJIT, skylake)
BLAS     : libopenblas64_.0.3.20.dylib (ilp64)
```

```@docs
Heom.print_logo(io::IO=stdout)
```
The output will be something like the following:
```
julia> Heom.print_logo()
                                   __
                                  /  \
 __     __                     __ \__/ __
|  |   |  |                   /  \    /  \
|  |   |  | ______   ______   \__/_  _\__/
|  |___|  |/  __  \ /  __  \ / '   \/     \
|   ___   |  |__)  |  /  \  |    _     _   |
|  |   |  |   ____/| (    ) |   / \   / \  |
|  |   |  |  |____ |  \__/  |  |   | |   | |
|__|   |__|\______) \______/|__|   |_|   |_|
```