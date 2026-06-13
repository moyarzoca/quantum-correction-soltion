# Eigenvalue tracking for a quantum soliton ODE

Finds eigenvalues of a second-order linear ODE depending on parameters `k` (wavenumber), `n` (quantum number), and `λ` (eigenvalue), and tracks them as `k` varies.

## Methodology

1. **Frobenius series** around the singular point `z=0` to obtain boundary conditions at `z₀ = 10⁻⁶`
2. **Numerical integration** (`NDSolve`) from `z₀` to `z₁ = 1-10⁻⁶`
3. **Root-finding** (`FindRoot`) on `F(z₁) = 0` to locate eigenvalues
4. **Eigenvalue tracking** across `k` using previous `λ` as initial guess

## Files

### `eigenvalues_finder.wl`
Main computation script. Run with:

```bash
wolframscript -file eigenvalues_finder.wl
```

#### Parameters (line 2)

| Variable | Value | Meaning |
|---|---|---|
| `orden` | 10 | Series expansion order |
| `prec` | 30 | Working precision (digits) |
| `z0` | 10⁻⁶ | Left boundary (near singular point z=0) |
| `z1` | 1-10⁻⁶ | Right boundary (near singular point z=1) |

#### Functions

| Function | Inputs | Output | Purpose |
|---|---|---|---|
| `P2[z]`, `P1[z]`, `P0[z]` | `z` | Polynomial | ODE coefficients: `P₂(z) F'' + P₁(z) F' + P₀(z) F = 0` |
| `bcsolu[λ,n,k]` | λ, n, k | Number | Series solution at z₀ (boundary condition for F) |
| `Dbcsolu[λ,n,k]` | λ, n, k | Number | Derivative of series at z₀ (BC for F') |
| `Equ[λ,n,k]` | λ, n, k | Expression | The ODE operator applied to F[z] |
| `SOL[λ,n,k]` | λ, n, k | NDSolve rules | Integrates the ODE, returns `{F -> InterpolatingFunction}` |
| `SOL3[λ,n,k]` | λ, n, k | Number | Evaluates the solution at z₁: `F(z₁)` |
| `BuscarModos[n,k,λmin,λmax,nPuntos]` | n, k, interval | List of λ | Scans for eigenvalues via sign changes + FindRoot |
| `trackModo[n,kStart,kEnd,dk,λSeed]` | n, k range, seed | List of `{k, λ}` | Tracks an eigenvalue continuously as k varies |

#### Execution flow

1. Build power series around `z=0`, solve linear system for coefficients `a[i]`
2. Use the series to define boundary conditions at `z₀`
3. Find seed eigenvalues at `kᵢ = 0.001` for `n = 1..4` via coarse scan (`BuscarModos`)
4. Track each eigenvalue from `k=0.001` to `k=0.01` in steps of `1/3000` (`trackModo`)
5. Save results to `evolucion_n{1..4}.m`

### `plot.wl`
Post-processing script. Loads the saved `.m` files and plots `λ(k)` for each `n` with linear fits.

## ODE

The differential equation is:

```
P₂(z) F''(z) + P₁(z) F'(z) + P₀(z; λ, n, k) F(z) = 0
```

with boundary condition `F(z₁) = 0` (Dirichlet at right boundary). The functions `P₀, P₁, P₂` are polynomials in `z` derived from the physical model, with `λ` entering linearly in `P₀`.
