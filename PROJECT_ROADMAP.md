# 🗺️ Project Roadmap — Three-Body Problem

> Tracking document for polishing this project to research-grade quality.
> Target: demonstrate physics + computational competency for OIST / MEXT Masters application.
>
> **Status legend:** ✅ Done · 🔄 In Progress · ⬜ Not Started

---

## Phase 1 — Quick Fixes
> ~30 minutes. Things that make the project look unfinished to any reviewer.

| # | Task | Status | Notes |
|---|------|--------|-------|
| 1.1 | Replace `[Your Name]` in README.md license block | ✅ | Bottom of README |
| 1.2 | Replace `yourusername` in README.md GitHub badge/link | ✅ | Top of README |
| 1.3 | Remove or add missing configs: README lists `euler_collinear` and `broucke_henon` but they don't exist in code | ⬜ | Either add them (Phase 3) or remove from README table |
| 1.4 | Document or remove `frontend/` — currently undocumented in README | ⬜ | Add a section explaining what it does and how to run it |
| 1.5 | Explain `simulation/index.html` — standalone HTML simulation, not mentioned anywhere | ⬜ | Add note in README or move to `frontend/` |
| 1.6 | Notebook 06 mismatch — `06_3d_visualization.ipynb` exists but physics engine is 2D | ⬜ | Either add 3D engine (Phase 5) or rename/scope the notebook correctly |

---

## Phase 2 — Visual Output
> ~1-2 hours. A professor cloning the repo should see results without running anything.

| # | Task | Status | Notes |
|---|------|--------|-------|
| 2.1 | Generate and commit reference plots for all 6 configs → `docs/figures/` | ⬜ | figure-eight, pythagorean, lagrange, solar_system, sun_earth_moon, perturbed |
| 2.2 | Embed 2-3 key images in README.md | ⬜ | figure-eight trajectory + Lyapunov convergence at minimum |
| 2.3 | Run all notebooks and save with output cells populated | ⬜ | Notebooks should show inline plots, not blank cells |
| 2.4 | Generate and commit Lyapunov convergence plot for Pythagorean | ⬜ | Core result of the project — must be visible |

---

## Phase 3 — Physics Completeness
> ~2-3 hours. Fixes gaps a physicist would notice.

| # | Task | Status | Notes |
|---|------|--------|-------|
| 3.1 | Add `euler_collinear` initial condition | ⬜ | Euler (1767) collinear solution — one of the three known exact families. Referenced in docs but missing from code. |
| 3.2 | Fix `angular_momentum_error()` division by zero | ⬜ | When L₀≈0 (figure-eight), function returns ~10¹⁹. Guard: if `|L₀| < 1e-10`, return absolute error. File: `src/analysis/analysis.py` |
| 3.3 | Fix Leapfrog docstring — "no drift" claim is misleading | ⬜ | Correct claim: *bounded oscillation*, not zero drift. Close encounters cause transient spikes. File: `src/core/integrator.py` |
| 3.4 | Fix `sun_earth_moon` Moon velocity inconsistency | ⬜ | Formula uses `3e-6` for Earth mass but config value is `3.003e-6`. File: `src/core/initial_conditions.py` |
| 3.5 | Add `broucke_henon` initial condition (optional) | ⬜ | Broucke (1975) stable family. Low priority — only if time permits. |

---

## Phase 4 — Original Experiments
> ~3-4 hours. This is what separates "I built a tool" from "I did science with it."
> **Highest impact for the application.**

| # | Task | Status | Notes |
|---|------|--------|-------|
| 4.1 | **Experiment: Lyapunov vs perturbation size** | ⬜ | Compute λ_L for `figure_eight_perturbed` across δ = [1e-4, 1e-3, 1e-2, 0.1] AU. Does λ_L depend on δ or converge to same value? 3-sentence conclusion. Add to notebook 03. |
| 4.2 | **Experiment: Integrator accuracy comparison** | ⬜ | Plot energy error vs time for all 4 integrators on same axes for figure-eight orbit. DOP853 flat, RK4 drifting, Leapfrog oscillating. Written interpretation. New section in notebook 04. |
| 4.3 | **Experiment: Chaos onset in Pythagorean** | ⬜ | Plot body separations over time. Mark ejection event. Show Lyapunov converges *before* ejection. Conclusion: chaos precedes ejection. Add to notebook 02. |
| 4.4 | **Notebook 05: quantitative bridge** | ⬜ | Add numerical comparison: your computed λ_L (Pythagorean) vs Schwarzschild photon sphere γ=π. Discuss why they differ in magnitude. This is the research payoff. |

---

## Phase 5 — Research-Grade Features
> ~4-6 hours. Makes this look like preparation for actual research.

| # | Task | Status | Notes |
|---|------|--------|-------|
| 5.1 | **3D physics engine** | ⬜ | Extend `equations.py` and `integrator.py` for 3D (positions shape `(3,)` instead of `(2,)`). Makes `06_3d_visualization.ipynb` legitimate. |
| 5.2 | **Poincaré section** | ⬜ | `plot_poincare_section()` is declared in `visualize.py` but not implemented. Standard tool for visualizing chaos in Hamiltonian systems. |
| 5.3 | **HDF5 save/load for SimulationResult** | ⬜ | h5py is in requirements.txt but unused. Add `result.save(path)` and `SimulationResult.load(path)`. Useful for long runs. |
| 5.4 | **Notebook 05: code-level bridge to Kerr** | ⬜ | End notebook 05 with a concrete code sketch of the Kerr geodesic integrator using the same `solve_ivp` call signature. Show the parallel in code, not just prose. |
| 5.5 | **Virial theorem verification** | ⬜ | For stable orbits, time-averaged `<KE> = -½<PE>`. Compute this and verify numerically. Adds another physics validation layer. |

---

## Phase 6 — Presentation Layer
> ~1 hour. How it looks on GitHub matters for an application reviewer.

| # | Task | Status | Notes |
|---|------|--------|-------|
| 6.1 | Add GitHub repo description + topics | ⬜ | Topics: `gravitational-dynamics`, `chaos`, `numerical-methods`, `python`, `astrophysics`, `three-body-problem` |
| 6.2 | Set repo social preview image | ⬜ | Use the figure-eight trajectory plot |
| 6.3 | Create `RESULTS.md` | ⬜ | One page: what each experiment showed, with inline images. Link this from CV and README. |
| 6.4 | Add GitHub Actions CI | ⬜ | Run `pytest tests/` on every push. Adds a green ✓ badge to README — signals working, tested code. |

---

## Already Completed ✅

| Task | What was done |
|------|--------------|
| Fix figure-eight initial conditions | Velocities were in G=1 units; scaled by 2π to G=4π² units. Energy error: 4.49e-3 → 5e-9 |
| Fix angular momentum test | L₀=0 for figure-eight causes ÷0. Test now uses absolute error for zero-L systems |
| Fix Leapfrog energy test | Replaced Lagrange triangle (close encounter spike) with figure-eight for fair comparison |
| Fix periodicity test threshold | 0.1 AU → 0.01 AU now that ICs are correct |
| Fix `.gitignore` | Added `frontend/node_modules/`, `*.egg-info/` |
| Polish `setup.py` | Replaced placeholder author/email, added classifiers, pinned deps |
| All 16 tests passing | `pytest tests/test_physics.py` → 16 passed in 12.5s |

---

## Current Test Status

```
pytest tests/test_physics.py -v
16 passed in 12.45s ✅
```

---

## Notes for Interview / Application

- Be ready to explain **why G = 4π²** (AU/yr/M☉ unit system makes G exact, no floating point factors)
- Be ready to explain **what DOP853 does** (adaptive 8th-order Dormand-Prince with 5th and 3rd order error estimates for step control)
- Be ready to explain **the state vector** (converts 3 second-order ODEs → 12 first-order ODEs for solve_ivp)
- Be ready to draw the **Lyapunov → photon ring** connection on a whiteboard
- Know the **three known exact solution families**: Lagrange (equilateral triangle), Euler (collinear), Figure-eight (Moore/Chenciner-Montgomery)
