# Eigenvalue tracking for a quantum soliton ODE

Finds eigenvalues λ of a second-order linear ODE depending on parameters `k` (wavenumber), `n` (Matsubara mode), and `λ` (eigenvalue), and tracks them as `k` varies.

## ODE

```
P₂(z) F''(z) + P₁(z) F'(z) + P₀(z; λ, n, k) F(z) = 0
```

with Dirichlet boundary condition `F(z₁) = 0` at the right endpoint. The coefficients `P₀, P₁, P₂` are polynomials in `z` derived from the physical model. `λ` enters linearly in `P₀`, and `n` enters up to quadratically in `P₀`.

The integration domain is `z ∈ [z₀, z₁]` where `z₀ ≈ 10⁻⁵` and `z₁ = 1 − offset`.

## Methodology

1. **Frobenius series** around the singular point `z = 0` to obtain approximate boundary conditions at `z₀`.
2. **Numerical integration** (`NDSolve` with `"StiffnessSwitching"`) from `z₀` to `z₁`.
3. **Interpolation-based root-finding** on `F(z₁) = 0` to locate eigenvalues. The residual `F(z₁)` is sampled on a λ-grid, interpolated, and the root of the interpolant is found — no initial guess needed.
4. **Two-stage refinement**: a coarser scan (medium stage) brackets the root, then a finer scan (fine stage) with a narrowed window refines it.
5. **Eigenvalue tracking** across `k` using scan-based continuation: at each new `k`, the previous λ defines the window center.

## File structure

```
solver/
├── config.wl              — Stage presets, λ-windows, and output config
├── eigenvalues.wl         — Main orchestrator script
├── plot.wl                — Post-processing: load .m files and plot λ(k) [Not yet wired]
├── tools/
│   ├── equation.wl        — ODE coefficients P₀,P₁,P₂, Equ, makeFrobeniusSolver
│   ├── export.wl          — setupRunDirectory (auto-increment numbered run dirs)
│   ├── modes.wl           — FindFirstLamb, RefineFindFirstLamb, makeKGrid, trackMode
│   └── solver.wl          — SOL, SOL3 (NDSolve wrapper)
```

Output goes to `output/<NNN>/` — each run gets its own numbered directory.
The `output/` directory is gitignored.

### `config.wl`

Defines the stage presets, λ-search windows, and output configuration:

| Preset | `prec` | `orden` | `z1Offset` | Use |
|--------|--------|---------|------------|-----|
| `mediumStage` | 35 | 5 | 10⁻⁶ | Coarse scan (fast) |
| `fineStage` | 38 | 7 | 10⁻⁷ | Fine scan (accurate) |

The `lambWindow` association controls the λ search:
- `"center"` — center of the search interval
- `"width"` — full width of the search interval
- `"Npoints"` — number of λ points to sample

The `outputConfig` association controls the output directory and file names:
- `"baseDir"` — parent directory (default `"output"`)
- `"evolutionsFile"`, `"fitsFile"`, `"logFile"` — file names inside each run directory

Each run creates an auto-incremented numbered subdirectory under `baseDir`
(e.g., `output/001/`, `output/002/`, …). Output files are gitignored.

### `tools/equation.wl`

- `P2[z, λ, n, k]`, `P1[z, λ, n, k]`, `P0[z, λ, n, k]` — ODE coefficient polynomials
- `Equ[λ, n, k]` — the ODE operator applied to `F[z]`
- `makeFrobeniusSolver[orden]` — computes the Frobenius series expansion around `z = 0` up to order `orden`. Returns an association `{ "solu" -> fn, "Dsolu" -> fn }` with functions that evaluate the series and its derivative at a given `z₀`.

### `tools/solver.wl`

- `SOL[λ, n, k, frob, stage]` — integrates the ODE with `NDSolve` using the given Frobenius solver and stage settings. Returns the NDSolve rule list.
- `SOL3[λ, n, k, frob, stage]` — evaluates the solution at `z₁`, returning `F(z₁)`. This is the residual function used for root-finding.

### `tools/export.wl`

- `setupRunDirectory[baseFolder, config]` — scans `output/` for existing numbered
  directories, creates the next one (zero-padded to 3 digits: `001`, `002`, …),
  and returns its path as `<|"dir" -> path|>`. Called from `eigenvalues.wl` at
  the start of each run.

### `tools/modes.wl`

- `FindFirstLamb[n, k, lambWindow, stage, frob]` — scans λ on a uniform grid, interpolates the residual `F(z₁)`, and finds the root via `FindRoot` on the interpolant. Returns `<|"lambda" -> λ, "grid" -> ...|>` or `Failure[...]`.
- `RefineFindFirstLamb[n, k, lambWindow, stage1, stage2, frob1, frob2]` — two-stage root finding: coarse scan with `stage1`/`frob1`, then narrow the window and refine with `stage2`/`frob2`.
- `makeKGrid[kStart, stepSize, nPoints]` — creates a decreasing k-grid from `kStart` downward.
- `trackMode[n, kGrid, lambWindow, stage, lambdaseeds, frob]` — scans through a k-grid, using `FindFirstLamb` at each step with the previous λ as the window center.

### `eigenvalues.wl`

Main orchestrator. Run with:

```bash
wolframscript -file eigenvalues.wl
```

Execution flow:
1. Loads all modules and config.
2. Calls `setupRunDirectory` to create `output/<NNN>/` with auto-increment.
3. Copies `config.wl` into the run directory as a snapshot for reproducibility.
4. Builds one Frobenius solver at the fine-stage order (shared across all n).
5. **Stage 1** — finds seed eigenvalues at `kInitial = 0.01` for `n = 1..2` using `RefineFindFirstLamb`.
6. Builds the k-grid: `kInitial` → `kInitial - step*nPoints` (downward).
7. **Stage 2** — tracks each mode across the k-grid using `trackMode` with scan-based continuation.
8. Exports `evolutions.wl` and `fits.wl` to the run directory. All `Print`
   output is also teed to `run.log`.

