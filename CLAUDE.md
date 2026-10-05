# HierarchicalEOM.jl

See [`.github/copilot-instructions.md`](.github/copilot-instructions.md) for full project context, domain concepts, structure, and conventions.

## Common Commands

```bash
# Install Runic formatter (one-time setup)
make setup

# Format code in-place with Runic
make format

# Run tests (4 threads, GROUP=All by default)
make test

# Build documentation
make docs
```

## Notes

- `make format` uses Runic.jl (not JuliaFormatter). The Makefile comment is stale — the actual tool is Runic.
- `make test` runs `Pkg.update()` before testing; use `julia -t 4 --project -e 'using Pkg; Pkg.test()'` directly if you want to skip the update.
- Docs require a separate `--project=docs` environment; `make docs` handles the `Pkg.develop` step automatically.
- CUDA extension tests run on Buildkite CI (not in local `make test`).
