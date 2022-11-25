# Documentation

First, run the following command to create the examples in markdown format ***if there is any updates in [notebooks](./src/examples/notebooks)***.
```sh
> julia nb2md.jl
```

Use the following command to build the documentation locally, and remember to install `Documenter.jl` before running it.

```sh
> julia make.jl local
```