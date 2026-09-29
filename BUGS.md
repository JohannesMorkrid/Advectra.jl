# Known bugs

Open bugs found while writing the tests in `TESTS_TODO.md`. Bugs 1–9 were found in
section D: each one was reproduced, and a fix was written and verified at the time. Those
changes were later reverted, so they are still open. Bug 10 was found in section E. Bugs 11
and 12 were found by running the GPU scripts in `test/gpu/` on a CUDA machine. All of them
already existed before this round of test work (commit `ac05e99`). Line numbers refer to the
current code.

Ordered roughly by impact.

---

## 1. Potential and enstrophy dissipation integrals have the wrong sign

**Where:** `src/diagnostics/energy_integrals.jl`, `potential_dissipation_integral` (l. 312)
and `enstrophy_dissipation_integral` (l. 391).

**What happens:** `energy_evolution_integral` and `enstrophy_evolution_integral` don't
match the actual rate of change of the energy and enstrophy in a simulation.

**Why:** With `hyper_laplacian` = ∇⁶ ↦ −k⁶, a hyper-diffusion term ν∇⁶n damps. The
potential dissipation is computed as `ν∫n∇⁶n`, which is **negative** for damping. The energy
evolution then subtracts it (`dE/dt = Γ − Γc − D^E`), so hyper-diffusion *adds* energy. The
kinetic dissipation `μ∫ϕ∇⁶Ω = μ∫(∇²Ω)² ≥ 0` and the resistive dissipation `C∫(n − ϕ)² ≥ 0`
are already positive, so only these two are inconsistent. `enstrophy_dissipation_integral`
has the same problem: `∫(n − Ω)(ν∇⁶n − μ∇⁶Ω)` is the (negative) rate, but it is subtracted.

**Evidence:** Hasegawa–Wakatani run (32², `C = 0.05, ν = 2e-2, μ = 3e-2`, damping `+ν∇⁶n`,
`+μ∇⁶Ω`, flux term `−∂ϕ/∂y`). The central-difference time derivatives are compared with the
integrals:

| | dE/dt | dW/dt, W = ½∫(n − Ω)² |
|---|---|---|
| Finite difference | −0.0167424 | 0.0178967 |
| Current code | −0.0161634 (3.5% off, exactly 2·D^E_N) | 0.0205551 (15% off) |
| With the sign fixed | −0.0167424 | 0.0178967 |

**Fix:** Return `−ν∫n∇⁶n` and `−∫(n − Ω)(ν∇⁶n − μ∇⁶Ω)`, so all dissipation terms are
positive and the evolution formulas hold as written.

**Related:** `enstrophy_energy_integral` computes ½∫Ω², while `enstrophy_evolution_integral`
is the evolution of the Hasegawa–Wakatani enstrophy ½∫(n − Ω)². The two are different quantities.
One of them should be renamed or changed.

---

## 2. Spectra from real transforms are wrong

**Where:** `src/diagnostics/spectral.jl`, `_spectral_sum` (l. 82) and `energy_spectrum`
(l. 121–150).

**What happens:** With `real_transform=true` (the default), the spectra from the spectrum
diagnostics (potential, kinetic, flux, enstrophy, electrostatic potential) are wrong:

- **Poloidal spectrum:** sums to about half of the corresponding integral, e.g. 0.0376 instead of
  0.0526 for the potential energy. The complex transform gives the right value.
- **Wavenumber spectrum:** differs from the complex-transform result by up to 20%.
- **Radial spectrum:** correct in total, but the individual E(kx) values differ from the
  complex-transform result.

**Why:** A real transform only stores ky ≥ 0. Every row except ky = 0 (and the Nyquist
row) also stands for the −ky modes (Hermitian symmetry Â(−k) = conj(Â(k))).
- Poloidal (sum over kx per ky) and wavenumber spectra don't account for this. They
  weight each stored mode once.
- The radial spectrum uses `spectral_sum(...; dims=1)`, which doubles the ky > 0 rows *within
  the same column*. But the −ky modes of column kx are stored in column **−kx**, not kx. So each
  E(kx) mixes kx and −kx. Summing over kx gives the right total, but the individual values are wrong.

**Fix:**
- Poloidal and wavenumber: weight the ky rows by 1 (ky = 0 and Nyquist) or 2 (all others) before
  summing or binning. For the wavenumber spectrum this means weighting the bin masks.
- Radial: in `_spectral_sum` with `dims=1` (over all kx), average each column with its mirror,
  `(S .+ circshift(reverse(S; dims=2), (0, 1))) ./ 2`.

With these, the radial and wavenumber spectra are identical for both transforms, and all
spectra sum to their integrals.

---

## 3. `probe_radial_velocity` always crashes

**Where:** `src/diagnostics/probe.jl:310`.

**What happens:** `UndefVarError: state_hat not defined`, whenever the diagnostic is called.

**Why:** The argument is named `state`, but the body uses `state_hat`.

**Fix:** Rename the argument to `state_hat`.

---

## 4. `probe_all` overwrites the state

**Where:** `src/diagnostics/probe.jl:360–366`.

**What happens:** With a real transform, calling `probe_all` changes the spectral state
passed to it. Inside a simulation, the diagnostic corrupts the solution.

**Why:** `mul!(cache, bwd_plan(domain), n_hat)` applies FFTW's inverse real transform
(c2r) in place. FFTW's multi-dimensional c2r transforms overwrite their input, and FFTW
can't be told to preserve it. `bwd_plan(domain) * n_hat` copies first, `mul!` doesn't.

**Fix:** Use `bwd_plan(domain) * n_hat` (and the same for Ω, ϕ, v_x) instead of `mul!`.

See also bug 8, which is the same problem in other operators.

---

## 5. Probe interpolation can't be used

**Where:** `src/diagnostics/probe.jl:39` and `:88`.

**What happens:** Passing `interpolation=linear_interpolation` (or anything else) to a
probe diagnostic gives a `MethodError` when the probe is evaluated.

**Why:** The method requires `interpolation::AbstractInterpolation` (an interpolant
*object*), but then calls it as `interpolation((domain.y, domain.x), field)`, i.e. as a
*constructor* such as Interpolations.jl's `linear_interpolation`. No value can satisfy both.

**Fix:** Type the argument as `interpolation::Function` (e.g. `linear_interpolation`,
`cubic_spline_interpolation`). With that change, linear interpolation of a linear field is exact.

---

## 6. The CFL number ignores negative velocities

**Where:** `src/diagnostics/CFL.jl:56`, `:61`, `:65–66` (`compute_cfl` for `:x`, `:y`, `:both`).

**What happens:** The reported maximum CFL number is the maximum of the *signed* velocity
times dt/dx. If the largest speed is in the negative direction, the CFL number is too small,
or even negative. `:magnitude` is not affected (it uses `hypot`).

**Fix:** Use `abs.(first(velocities))` and `abs.(last(velocities))`.

---

## 7. The last time step can be skipped

**Where:** `src/spectralSolve.jl:16` (`floor`) and `src/outputer.jl:466` (`ceil`).

**What happens:** For `tspan = [0, 0.3]` and `dt = 0.1`, the solver stops at t = 0.2 and
the last reserved sample is never written. Reading `"Density/t"` returns `[0.0, 0.1, 0.2, 0.0]`.

**Why:** `0.3 / 0.1 == 2.9999999999999996`. The solver uses `floor` (2 steps), while the
output uses `ceil` (3 steps) to reserve space for the samples. Common values such as
`0.1 / 0.001` happen to be exact, which is why this went unnoticed.

**Fix:** Use one function for both. It should round when `T/dt` is within floating-point
error of an integer, and round up otherwise:

```julia
N = (last(tspan) - first(tspan)) / dt
isapprox(N, round(N)) ? round(Int, N) : ceil(Int, N)
```

---

## 8. Without dealiasing, operators overwrite their inputs

**Where:** Every `mul!(U, bwd(transforms), padded ? pad!(up, u, ...) : u)`. That includes
`quadratic_term` (`src/operators/quadraticTerm.jl:62`), the spectral functions
(`src/operators/spectralFunctions.jl:48`), `integral_of_quadratic_term`
(`src/diagnostics/energy_integrals.jl:47`), and the non-Boussinesq `solve_phi`
(`src/operators/solvePhi.jl`, e.g. l. 151–152).

**What happens:** With `dealiased=false` and a real transform, these operators overwrite the
spectral arrays passed to them. For example, `solve_phi(n_hat, ϖ_hat)` changes `n_hat` and
`ϖ_hat`, and `quadratic_term(n_hat, Ω_hat)` would corrupt the state. With dealiasing (the
default), the input is first copied into a padded buffer, so it is safe.

**Why:** Same as bug 4. FFTW's c2r transforms overwrite their input and can't be told not to.

**Fix:** When not padded, copy the input into a scratch buffer before the inverse transform. For
example, allocate `up`/`vp` without padding too and always copy into them. This costs one
copy per transform, only in the non-dealiased case.

---

## 9. Two-argument `parseval_integral` returns a wrong imaginary part

**Where:** `src/diagnostics/energy_integrals.jl:36`.

**What happens:** For real transforms, `parseval_integral(a_hat, b_hat, domain)` returns a
complex number. The real part is correct, but the imaginary part is **not** round-off: 13.3 for
random fields. The flux (already fixed), dissipation and evolution integrals all go through this
function. As a result, `energy_evolution_integral` etc. return `Complex` values.

**Why:** `spectral_sum` doubles the ky > 0 rows to account for the missing −ky modes.
For a product `â·conj(b̂)`, the mirror mode contributes the complex conjugate, so the correct
factor is `2·Re(…)`. Doubling the complex value keeps a spurious imaginary part.

**Fix:** Return `real(integral)` when the physical fields are real
(`physical_eltype(domain) <: Real`, i.e. real transforms).

---

## 10. MSS3 is only second order

**Where:** `src/schemes.jl`, the first start-up step of `perform_step!` for `MSS3Cache` and
`MSS3ConstantCache`.

**What happens:** MSS3 converges at second order, not third. For du/dt = ν∇²u + au (a single
Fourier mode, dt = 0.1 → 0.05 → 0.025), the errors are about the same as MSS2's:

| | dt = 0.1 | dt = 0.05 | dt = 0.025 | order |
|---|---|---|---|---|
| MSS2 | 2.5e-3 | 6.1e-4 | 1.5e-4 | 2.0 |
| MSS3 | 2.6e-3 | 6.9e-4 | 1.8e-4 | 1.9–2.0 |
| MSS3, exact start-up values | 3.8e-5 | 5.4e-6 | 7.1e-7 | 2.8–2.9 |

**Why:** MSS3 needs two previous values, so it starts with one MSS1 step and one MSS2 step.
The MSS1 step (backward-Euler-like) has a local error of O(dt²), and that error stays in the
solution for the rest of the run. So the whole run is only O(dt²) accurate, even though every
later step is third order.

**Fix:** Make the first start-up value accurate to O(dt³). The MSS2 step after it is already
accurate enough. Replace the single MSS1 step with Richardson extrapolation:
u₁ = 2·(two MSS1 steps of dt/2) − (one MSS1 step of dt). This costs two extra cheap steps at
the start only, about 25–30 lines:
- `get_cache(prob, MSS1())` accepts an optional `dt`, since `c = (1 − D·dt)⁻¹` depends on it.
- A shared helper does the three MSS1 steps and combines them.
- The MSS3 caches (in-place and out-of-place) call it in their first start-up step.

Sub-stepping the start-up with plain MSS1 would need about 1/dt sub-steps, so Richardson
extrapolation is the better option.

**Test:** `test/schemes_tests.jl` currently checks that MSS3 is second order in a normal run,
and third order with exact start-up values. Once this is fixed, the normal-run check should
be changed to third order.

---

## 11. `memory_type` clashes with CUDA.jl

**Where:** `src/Advectra.jl` (export list) and `src/domains/domain.jl` (`memory_type`).

**What happens:** After `using Advectra, CUDA`, calling `memory_type` fails with
`UndefVarError: memory_type not defined`. Julia adds the hint that two modules export
different bindings with this name. The usual way to move data to the GPU,
`u0 |> memory_type(domain, Physical())`, therefore breaks exactly when CUDA is loaded. This is
also why 6 of the 7 GPU scripts fail at their first `memory_type` call.

**Why:** CUDA.jl exports its own `memory_type` function. If two modules loaded with `using`
export the same name, Julia refuses to pick one, so the name can only be used qualified
(`Advectra.memory_type`). The export was added in #121 ("Domain refactoring"). The same PR
changed the GPU scripts from `Advectra.memory_type(domain)`, which works alongside CUDA, to the
unqualified name.

**Workaround:** Write `Advectra.memory_type(...)`, or `import Advectra: memory_type`.

**Fix (API decision):** Either rename the function (e.g. `field_type` or `array_type`) or
stop exporting it and document the qualified form. Both change the public API.

---

## 12. The GPU scripts in `test/gpu/` are broken

**Where:** `test/gpu/*.jl`, run with `ADVECTRA_TEST_GPU=true`.

**What happens:** Plain run on an RTX 3050 (`CUDA.functional()` is true): 6 of 7 scripts fail
at `memory_type` (bug 11). With `memory_type` imported explicitly, 4 still fail:

| Script | Failure |
|---|---|
| `cfl_tests.jl:30` | `Unknown component: :something`: a deliberately invalid component, not wrapped in `@test_throws` |
| `spectral_tests.jl:35` | `axis has to be either :kx, :ky, :both or :diag`: a deliberately invalid `axis=:all`, not wrapped in `@test_throws` |
| `energy_integrals_tests.jl:69`, `:87` | `MethodError`: the removed `compute_density` keyword of `parsevals_theorem` (now `average`) |
| `probe_tests.jl:17` | `FieldError: type Tuple has no field solve_phi`: the script passes `operators=()`, but `probe_all` needs `solve_phi` and `diff_y` |

`COM_tests.jl`, `fluxes_tests.jl` and `vorticity_tests.jl` run through. Everything before the
failing lines ran on the GPU without errors.

**Why:** These are the old scratch scripts, moved unchanged in item B6 of `TESTS_TODO.md`.
They contain no `@test`s, so even when they run, they only show that nothing throws, not that
the GPU results are correct. CUDA.jl and Plots.jl are also not test dependencies, so
`ADVECTRA_TEST_GPU=true` fails at `using CUDA` unless they are added by hand.

**Fix:** Replace the scripts with a small GPU suite: rerun a selection of the CPU tests with
`MemoryType=CuArray` and check that the results match the CPU ones. Give it its own
`test/gpu/Project.toml` with CUDA.jl, so it can be run with a single command.
