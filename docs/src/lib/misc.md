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
    Neill Lambert, Po-Chen Kuo and Shen-Liang Yang

Package informations:
===================================
Heom              Version: 0.1.0
DistributedArrays Version: 0.6.6
JLD2              Version: 0.4.25
LinearSolve       Version: 1.26.0
OrdinaryDiffEq    Version: 6.27.2
ProgressMeter     Version: 1.7.2
SnoopPrecompile   Version: 1.0.1

System informations:
===================================
Julia Version: 1.8.0
OS       : macOS (x86_64-apple-darwin21.4.0)
CPU      : 12 × Intel(R) Core(TM) i7-8700B CPU @ 3.20GHz
Memory   : 16.0 GB
WORD_SIZE: 64
LIBM     : libopenlibm
LLVM     : libLLVM-13.0.1 (ORCJIT, skylake)
BLAS     : libblastrampoline (f2c_capable)
```

```@docs
Heom.print_logo(io::IO=stdout)
```
The output will be something like the following:
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
```