# GitHub Copilot Instructions for HierarchicalEOM.jl

## Project Overview

`HierarchicalEOM.jl` is a Julia package implementing the Hierarchical Equations of Motion (HEOM) approach for simulating open quantum systems, including non-Markovian effects. It is built on top of [`QuantumToolbox.jl`](https://github.com/qutip/QuantumToolbox.jl) and the SciML ecosystem.

## Language & Runtime

- **Julia 1.10+** is required.
- All source code is Julia. Do not suggest Python or other languages.
- Use Julia idioms: multiple dispatch, abstract types, `@reexport`, in-place mutations (`!` suffix), and broadcasting.

## Code Style

- Code formatting is enforced by [Runic.jl](https://github.com/fredrikekre/Runic.jl). Always produce Runic-compliant code.
- 4-space indentation. No trailing whitespace.
- Type annotations on function signatures are encouraged for performance-critical code.
- Avoid unnecessary comments; prefer self-documenting names.

## Project Structure

```
src/
  HierarchicalEOM.jl          # Main module, all includes
  HeomBase.jl                 # Core types and utilities
  Parity.jl                   # Parity operators
  ADOs.jl                     # Auxiliary Density Operators
  bath/
    BathBase.jl               # Abstract bath types
    BosonBath.jl              # Bosonic bath
    FermionBath.jl            # Fermionic bath
  bath_correlation_functions/ # Bath correlation function helpers
  heom_matrices/
    heom_matrix_base.jl       # Abstract HEOM matrix types
    Nvec.jl                   # Hierarchy index vectors
    HierarchyDict.jl          # Hierarchy bookkeeping
    M_S.jl                    # System Liouvillian
    M_Boson.jl                # Boson HEOM matrix
    M_Fermion.jl              # Fermion HEOM matrix
    M_Boson_Fermion.jl        # Mixed boson-fermion HEOM matrix
  evolution.jl                # Time evolution solvers
  evolution_propagator.jl     # Propagator-based evolution
  steadystate.jl              # Steady-state solvers
  power_spectrum.jl           # Power spectrum calculation
  density_of_states.jl        # Density of states
  correlations.jl             # Two-time correlation functions
  deprecated.jl               # Backward-compat shims
ext/
  HierarchicalEOM_CUDAExt.jl  # Optional CUDA GPU support
```

## Key Concepts

- **Bath types**: `BosonBath`, `FermionBath` — wrap correlation functions as sums of exponentials.
- **HEOM matrices**: `M_S`, `M_Boson`, `M_Fermion`, `M_Boson_Fermion` — sparse superoperators built from the system Hamiltonian and bath parameters.
- **ADOs**: `ADOs` — the full auxiliary density operator state vector used in time evolution.
- **Solvers**: `HEOMsolve` for time evolution, `steadystate` for steady state, `PowerSpectrum`, `DensityOfStates`, `correlation_3op_2t` etc.
- Operators follow QuantumToolbox.jl conventions (`QuantumObject`, `Qobj`).

## Dependencies

Key packages (see `Project.toml`):
- `QuantumToolbox` — quantum objects, Liouvillians, superoperators
- `SciMLBase`, `OrdinaryDiffEqCore`, `OrdinaryDiffEqLowOrderRK` — ODE solvers
- `LinearSolve` — linear system solvers (steady state)
- `SciMLOperators` — lazy operator algebra
- `DiffEqCallbacks` — ODE callbacks
- `FastExpm` — fast matrix exponential
- `IncompleteLU` — ILU preconditioner for iterative solvers
- `SparseArrays` — all HEOM matrices are sparse

## Testing & Quality

- Tests live in `test/`. Each major component has its own file.
- Run tests with `julia --project=test test/runtests.jl` or via `Pkg.test()`.
- Code quality is checked with [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) and [JET.jl](https://github.com/aviatesk/JET.jl).
- Spell checking uses [typos](https://github.com/crate-ci/typos); `ket` is an allowed exception (see `.typos.toml`).

## Contributing Guidelines

- Fork the repo and submit PRs against `main`.
- Follow the [QuantumToolbox.jl Contributor's Guide](https://qutip.org/QuantumToolbox.jl/stable/resources/contributing).
- Prefer sparse matrix operations; avoid dense conversions in HEOM matrix construction.
- New bath types should subtype the abstract bath hierarchy in `BathBase.jl`.
- New HEOM matrices should subtype `AbstractHEOMLSMatrix` from `heom_matrix_base.jl`.
- GPU support is handled via the CUDA extension; core code must remain CPU-compatible.

## Citation

If you generate code that uses this package, remind users to cite:
> Huang et al., *Communications Physics* **6**, 313 (2023). https://doi.org/10.1038/s42005-023-01427-2
