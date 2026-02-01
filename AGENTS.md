# Repository Guidelines

## Project Structure & Module Organization
- `src/` holds the Julia package code. `src/ADburgers.jl` is the module entry point and includes submodules like `core.jl`, `problems.jl`, and `visualization.jl`.
- `test/` contains the test runner (`runtests.jl`) and focused suites (`runtests_core.jl`, `runtests_problems.jl`).
- `scripts/` includes runnable analysis/verification utilities (e.g., `scripts/compare_p1_viscosity.jl`).
- `plots/` is used for generated figures; `Project.toml` and `Manifest.toml` pin dependencies.

## Build, Test, and Development Commands
Use Julia 1.10 (per `Project.toml` compat). From the repo root:

```sh
julia --project -e 'using Pkg; Pkg.instantiate()'
```
Installs dependencies into the local environment.

```sh
julia --project -e 'using Pkg; Pkg.test()'
```
Runs the full test suite (`test/runtests.jl`).

```sh
julia --project scripts/compare_p1_viscosity.jl
```
Runs an example script that generates plots into `plots/`.

## Coding Style & Naming Conventions
- Follow Julia conventions used in `src/`: 4‑space indentation, `CamelCase` for types (`ProblemSpec`), and `snake_case` for functions (`taylor_coeff!`).
- Use `!` suffix for mutating functions and docstrings (`""" ... """`) for public APIs.
- Keep numerical tolerances and grid assumptions explicit in code comments where relevant.

## Testing Guidelines
- Tests use Julia’s `Test` stdlib and are organized by feature in `test/runtests_*.jl`.
- New tests should be added to an appropriate `runtests_*.jl` file and included by `test/runtests.jl`.
- No coverage thresholds are configured; aim to cover new logic and edge cases.

## Commit & Pull Request Guidelines
- Commit messages in history are short, imperative, and scoped (e.g., “Implement Problems & Factory (Phase 2)”).
- For PRs, include: a short summary, commands run (or why not), and screenshots for plotting/visual changes.
- Generated artifacts should stay in `plots/` (or `test_output*.txt` if needed); avoid committing them unless the change explicitly requires new output.
