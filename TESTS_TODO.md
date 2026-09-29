# Test suite TODO

Status at time of writing: 289 tests pass, but only 5 of 15 test files are included in
`runtests.jl` and line coverage is 34.6%. Time estimates are rough and assume familiarity
with the codebase; "h" = focused working hours.

Suggested order: A → B → C → D → E → F → G.

## A. Source bugs found while reviewing (fix before writing tests for these areas)

- [x] **A1. Spectrum defaults use `Val{:radial}` (a type) instead of `Val(:radial)`.**
  `kinetic_energy_spectrum`, `flux_spectrum`, `enstrophy_spectrum` and
  `electrostatic_potential_spectrum` throw a `MethodError` when called without a spectrum
  (`src/diagnostics/spectral.jl:218`, `:257`, `:295`, `:327`). — **15 min**
- [x] **A2. Profiles slice the wrong dimension.** All functions in
  `src/diagnostics/profiles.jl` use `ndims(prob.domain)` instead of `ndims(prob.domain) + 1`
  (`radial_density_profile` returns a 2×1 matrix). `radial_flux_profile` calls the
  non-existent `vExB`. None of them has a `build_diagnostic`, so they can't be used in
  `@diagnostics`. — **1–2 h**
- [x] **A3. CFL never requests its operators.** `requires_operator(::Val{cfl}; velocity_method)`
  should be `requires_operator(::Val{:cfl}; velocity=:ExB, kwargs...)`
  (`src/diagnostics/CFL.jl:155`). — **15 min**
- [x] **A4. ComponentArrays extension has the wrong argument order.**
  `_spectral_transform!(du, u::ComponentArray, p)` should be `(du, p, u::ComponentArray)`,
  including the inner call (`ext/AdvectraComponentArraysExt.jl`). — **30 min** (+ test in F6)
- [x] **A5. Delete dead `src/rhs.jl`.** It isn't included anywhere and doesn't compile. —
  **5 min**
- [x] **A6. `Domain` has no explicit size check.** `Domain(-64)` only throws via
  `LinRange`. Add `Nx > 0 && Ny > 0 || throw(ArgumentError(...))`. — **10 min**
- [x] **A7. `SpectralODEProblem` without `p` can't be written to HDF5.** The default
  `NullParameters` has no `keys`, so `Output` throws a `MethodError` when writing the
  parameter attributes. Workaround in `test/profile_tests.jl` (passes a dummy `p`). —
  **30 min**
- [x] **A8. `Output(prob; store_hdf=false)` crashes.** `setup_hdf5_storage` calls
  `rm(simulation.file.filename)` even when `simulation` is `nothing`
  (`src/outputer.jl:146`). — **15 min**
- [x] **A9. `radial_flux`/`poloidal_flux` return a complex number** (tiny imaginary part from
  `parseval_integral` of `a .* conj(b)`). Should return `real(...)`. — **15 min**

## B. Test infrastructure

- [x] **B1. Restructure `runtests.jl`.** Include every test file, each in its own
  `@testset`/module (or use `SafeTestsets.jl`) so files can't depend on each other's
  imports. — **1 h**
- [x] **B2. Clean up test dependencies.** Use only `test/Project.toml` (drop `[extras]` from
  the main `Project.toml`) and add HDF5, ComponentArrays and SMTPClient. — **30 min**
- [x] **B3. Run integration tests with `debug=true`.** `spectral_solve` swallows exceptions
  by default, so a crashing solve can still pass. — **10 min**
- [x] **B4. Write test output to a temp dir.** Use `mktempdir()` (or `store_hdf=false`)
  instead of writing `.h5` files into `test/output/`. — **20 min**
- [x] **B5. Run plot tests headless.** Set `ENV["GKSwstype"] = "100"` at the top of
  `runtests.jl`, before Plots/Advectra is loaded, so the 21 plots in `display_tests.jl`
  render off-screen instead of opening windows or failing on CI. — **10 min**
- [x] **B6. Move GPU tests to `test/gpu/`.** Run them only when `CUDA.functional()` is true
  or an env flag is set, and keep them out of the default CI run. — **1 h**
- [x] **B7. Rename `progressbar_test.jl` → `progress_tests.jl`** for consistent naming. —
  **5 min**
- [x] **B8. Delete `test/testutilities.jl`.** It is completely outdated (old `Output` API,
  `domain.SC`, Roots and PlotlyJS). Port the two convergence helpers into C1 if they are
  useful. — **15 min**

## C. Fix the tests that currently run

- [x] **C1. `integration_tests.jl`.** Compare linear diffusion against the exact solution
  `û(t) = û₀·exp(ν k² t)` instead of hard-coded matrices. Keep the non-linear case as a
  regression test, but drop the redundant second assertion and the `println`. — **2 h**
- [x] **C2. `progress_tests.jl`.** Fix the loop time (`i*dt` → `i`), remove the `sleep`s
  and `@test true`, and import `build_diagnostic` explicitly. — **20 min**
- [x] **C3. `operator_tests.jl`.** Remove the `try`/`catch` around GradDotGrad and use a case
  with a nonzero expected result. — **45 min**
- [x] **C4. `domain_tests.jl`.** Add the missing `@test` on line 59, check the
  `lengths`/`differential_elements` order with Lx ≠ Ly, and test the explicit size check
  from A6. — **30 min**

## D. Rewrite the scratch test files that never run (CPU, real assertions)

The original GPU scripts now live in `test/gpu/` (run with `ADVECTRA_TEST_GPU=true`). Write
the CPU versions as new files in `test/` and add them to `TEST_FILES` in `runtests.jl`.

- [ ] **D1. `cfl_tests.jl`.** Cover every component mode × ExB/burger velocity, compare
  against a hand-computed CFL for a known velocity field, and check the `silent` flag and
  the metadata. — **2 h**
- [ ] **D2. `COM_tests.jl`.** The COM of a shifted Gaussian equals its centre, the velocity
  equals Δx/Δt across two calls, and there is no division by zero at equal times. —
  **1 h**
- [ ] **D3. `fluxes_tests.jl`.** Radial and poloidal flux for an analytic `n` and `ϕ`, plus
  `flux_magnitude`. — **1.5 h**
- [ ] **D4. `energy_integrals_tests.jl`.** Fix the outdated `compute_density` → `average`
  keyword. Parseval gives the same result for real and complex transforms, and a constant
  field averages to 1. Check the kinetic/potential/total/enstrophy integrals against
  physical-space sums, check all dissipation/evolution integrals with custom coefficient
  symbols, and cover `integral_of_quadratic_term`. — **3–4 h**
- [ ] **D5. `probe_tests.jl`.** Errors for a wrong tuple length and for out-of-domain
  positions, type promotion, the value at a grid point, interpolation, all 5 probe
  types, and `(Ny, Nx)` vs `(Ny, Nx, 1)` shapes. — **2–3 h**
- [ ] **D6. `spectral_tests.jl`.** Cover `get_modes` for every axis (single- and
  multi-field), `get_log_modes`, and every spectrum type for each spectrum diagnostic.
  Invalid options must throw. Turn the `spectral_sum` exploration (lines 68–155) into
  real assertions. — **3 h**
- [ ] **D7. `vorticity_tests.jl`.** Rewrite on CPU with a small grid and no plots: the
  Boussinesq `solve_phi` inverts the Laplacian, and the non-Boussinesq version satisfies
  `∇·(N∇ϕ) ≈ ϖ` for both `density=:linear` and `density=:log`. — **2–3 h**
- [ ] **D8. `output_tests.jl`.** Rewrite against `parse_storage_limit`,
  `determine_sampling_strategy`, `validate_stride`, `recommend_stride`,
  `nearest_divisor`/`next_divisor` and `format_bytes`, including the edge cases listed
  in the file. — **2–3 h**

## E. New tests for untested core functionality

- [ ] **E1. Time-stepper convergence.** Run MSS1, MSS2 and MSS3 on linear diffusion, halve
  `dt`, and check that the error slopes are ≈ 1/2/3. Cover in-place and out-of-place
  right-hand sides (the out-of-place caches are untested). — **3–4 h**
- [ ] **E2. Checkpoint/resume round trip.** Solve to T/2, resume to T, and compare with
  solving straight to T. `resume=true` must reject a mismatched domain
  (`validate_resume_attributes`). — **2–3 h**
- [ ] **E3. HDF5 output layout.** Check the groups, attributes, strides, `storage_limit`,
  `simulation_name` handling and `physical_transform`. — **2 h**
- [ ] **E4. Operators.** Poisson bracket on analytic fields, `quadratic_term` vs a
  physical-space product with and without dealiasing, the `reciprocal`/`spectral_exp`/
  `spectral_expm1`/`spectral_log` functions, and `SpectralConstant`/`Source`
  (`sources.jl`, 0% covered). — **3 h**
- [ ] **E5. `@diagnostics` macro and `required_operators`.** Cover vector, block and single
  forms, check that the alias syntax errors, and that operators are pulled in
  (regression test for A3). — **1 h**
- [ ] **E6. Profiles.** Profile values for an analytic field, plus construction through
  `@diagnostics` (after A2). — **1 h**
- [ ] **E7. Sample diagnostics.** Cover density, vorticity, temperature and potential
  (`sample.jl`, 14% covered). — **45 min**
- [ ] **E8. `SpectralODEProblem`.** Operator recipes (`:default`, `:all`, custom),
  `convert_parameters` precision conversion, `isinplace` detection and `show`. — **2 h**

## F. Utilities and extensions

- [ ] **F1. Initial conditions.** Each `initial_condition` function gives the right shape
  and element type, including `isolated_blob`/`isolated_temperature_blob` `:lin`/`:log`
  and `@nobroadcast`. — **1.5 h**
- [ ] **F2. Mode removal.** `remove_zonal/streamer/asymmetric/nyquist_modes!` zero the
  correct entries for real and complex transforms. — **1 h**
- [ ] **F3. `add_constant(!)`, `logspace`, `spectral_sum` edge cases.** — **30 min**
- [ ] **F4. Float32 coverage.** Parameterise the core operator and solver tests over
  `Float32` and `Float64`. — **1 h**
- [ ] **F5. SMTPClient extension.** Test with the network send stubbed out. — **1 h**
- [ ] **F6. ComponentArrays extension.** A ComponentArray state through
  `SpectralODEProblem` and `spectral_solve` (regression test for A4). — **1 h**

## G. CI and tooling

- [ ] **G1. Add `Aqua.jl`.** It checks for method ambiguities, stale dependencies and
  missing compat bounds, and will likely flag the duplicate 3-argument
  `potential_energy_spectrum` methods. — **1 h** (+ fixing what it finds)
- [ ] **G2. Upload coverage to Codecov** (`julia-processcoverage` + `codecov-action`) and add
  a badge to the README. — **45 min**
- [ ] **G3. Run doctests in CI** (`doctest(Advectra)` in the test suite or a docs job). —
  **30 min** (+ fixing failing doctests)
- [ ] **G4. Optional: `JET.jl`** type-stability/error analysis. — **1–2 h**
- [ ] **G5. Optional: move Plots to a package extension** (existing `# TODO make ext`),
  which speeds up loading and tests. — **3–4 h**

## Total estimate

| Section | Estimate |
|---|---|
| A. Source bugs | ~3.5–4.5 h |
| B. Infrastructure | ~3.5 h |
| C. Fix running tests | ~3.5 h |
| D. Rewrite scratch tests | ~17–21 h |
| E. New core tests | ~15–17 h |
| F. Utilities & extensions | ~6 h |
| G. CI & tooling | ~2.5 h required, +4–6 h optional |
| **Total** | **~51–58 h required, plus ~4–6 h optional** |
