# Documentation

First, run the following command to create the examples in markdown format ***if there is any updates in [notebooks](./src/examples/notebooks)***.  
(Note that in order to execute the notebooks locally, remember to add a cell `import Pkg; Pkg.activate("../../../")` in the beginning of the notebooks.)
```sh
> julia nb2md.jl
```

Use the following command to build the documentation locally, and remember to install `Documenter.jl` before running it.

```sh
> julia make.jl local
```