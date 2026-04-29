# 🌌 Three-Body Problem — Gravitational Simulation Suite

> A research-grade Python package for simulating, visualizing, and analyzing the gravitational three-body problem — built as a foundation for understanding the mathematics behind black hole photon ring theory.

<div align="center">

![Python](https://img.shields.io/badge/python-3.11-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Tests](https://img.shields.io/badge/tests-pytest-orange)
![Physics](https://img.shields.io/badge/physics-Newtonian%20Gravity-purple)

</div>

---

## 📋 Table of Contents

1. [Project Overview](#1-project-overview)
2. [Scientific Background](#2-scientific-background)
3. [Mathematical Formulation](#3-mathematical-formulation)
4. [Installation](#4-installation)
5. [Project Structure](#5-project-structure)
6. [Quick Start](#6-quick-start)
7. [Initial Condition Catalogue](#7-initial-condition-catalogue)
8. [Integrators](#8-integrators)
9. [Analysis Tools](#9-analysis-tools)
10. [Notebooks Guide](#10-notebooks-guide)
11. [CLI Reference](#11-cli-reference)
12. [Running Tests](#12-running-tests)
13. [Connection to Black Hole Research](#13-connection-to-black-hole-research)
14. [References](#14-references)
15. [License](#15-license)

---

## 1. Project Overview

This package simulates the gravitational three-body problem — the motion of three point masses interacting through Newton's law of universal gravitation. The three-body problem has **no general closed-form solution** (Poincaré, 1890) and exhibits **deterministic chaos** for most initial conditions.

### What This Package Does

| Feature | Description |
|---------|-------------|
| 🔭 **Physics engine** | Full Newton gravity with softening, CM frame transform |
| 🧮 **4 integrators** | RK45, DOP853 (scipy), Leapfrog (symplectic), RK4 |
| 🌀 **6 configurations** | Figure-8, Lagrange, Pythagorean, Solar System, and more |
| 📊 **Analysis tools** | Lyapunov exponent, conservation law monitoring, ejection detection |
| 🎬 **Visualization** | Trajectories, phase space, animations, energy plots |
| 📓 **5 Notebooks** | Progressive tutorials from basics to research connection |
| ✅ **Test suite** | Physics-validated unit tests with pytest |

### Why This Project Was Built

This project serves two purposes:

1. **Learning**: Build deep intuition for gravitational dynamics, numerical integration, chaos theory, and conservation laws — all prerequisites for general relativistic ray tracing of Kerr black hole shadows.

2. **Research bridge**: The **Lyapunov exponent** computed here for chaotic three-body orbits is mathematically identical to the **orbital Lyapunov exponent γ** of the photon sphere in Kerr spacetime — the quantity that controls photon ring widths in the EHT M87* image.

---

## 2. Scientific Background

### The Three-Body Problem

The three-body problem asks: given three point masses with known initial positions and velocities, predict their future motion under mutual gravitational attraction.

**Why it's hard:**
- Two-body problem: completely solved by Kepler (1609) and Newton (1687)
- Three-body problem: proven to have no general closed-form solution by Bruns (1887) and Poincaré (1890)
- Poincaré's proof was the first rigorous demonstration of **deterministic chaos** in classical physics

**Degrees of freedom:** 3 bodies × 2D × 2 (position + velocity) = **12-dimensional phase space**

**Known exact solutions (special cases):**
- Lagrange (1772): equilateral triangle configuration
- Euler (1767): collinear configuration
- Figure-eight choreography (Moore 1993, Chenciner & Montgomery 2000)
- Šuvakov & Dmitrašinović (2013): 13 new periodic families

### Deterministic Chaos

The three-body problem is **chaotic**: small perturbations in initial conditions grow **exponentially** in time. This is quantified by the **Lyapunov exponent** λ_L:

```
|δz(t)| ≈ |δz(0)| · exp(λ_L · t)
```

- `λ_L > 0`: chaotic orbit (trajectories diverge exponentially)
- `λ_L ≈ 0`: regular/periodic orbit (trajectories stay bounded)

This makes long-term prediction of chaotic three-body systems fundamentally impossible — not due to computational limits, but due to the physics of exponential error amplification.

---

## 3. Mathematical Formulation

### 3.1 Equations of Motion

Newton's second law applied to three gravitating masses:

$$\frac{d^2\mathbf{r}_1}{dt^2} = G\frac{m_2(\mathbf{r}_2-\mathbf{r}_1)}{|\mathbf{r}_2-\mathbf{r}_1|^3} + G\frac{m_3(\mathbf{r}_3-\mathbf{r}_1)}{|\mathbf{r}_3-\mathbf{r}_1|^3}$$

$$\frac{d^2\mathbf{r}_2}{dt^2} = G\frac{m_1(\mathbf{r}_1-\mathbf{r}_2)}{|\mathbf{r}_1-\mathbf{r}_2|^3} + G\frac{m_3(\mathbf{r}_3-\mathbf{r}_2)}{|\mathbf{r}_3-\mathbf{r}_2|^3}$$

$$\frac{d^2\mathbf{r}_3}{dt^2} = G\frac{m_1(\mathbf{r}_1-\mathbf{r}_3)}{|\mathbf{r}_1-\mathbf{r}_3|^3} + G\frac{m_2(\mathbf{r}_2-\mathbf{r}_3)}{|\mathbf{r}_2-\mathbf{r}_3|^3}$$

### 3.2 State Vector & ODE System

Converting the three 2nd-order equations into 12 coupled 1st-order ODEs:

```
State vector:
y = [x₁, y₁,  x₂, y₂,  x₃, y₃,  ẋ₁, ẏ₁,  ẋ₂, ẏ₂,  ẋ₃, ẏ₃]
     ← positions (6) →         ← velocities (6) →

dy/dt = [ẋ₁, ẏ₁, ẋ₂, ẏ₂, ẋ₃, ẏ₃, aₓ₁, aᵧ₁, aₓ₂, aᵧ₂, aₓ₃, aᵧ₃]
```

### 3.3 Physical Units

All simulations use **astronomical units**:

| Quantity | Unit |
|---------|------|
| Length | AU (astronomical unit) |
| Time | yr (Julian year) |
| Mass | M☉ (solar mass) |
| G | 4π² AU³/(M☉·yr²) |

These units make G = 4π² exactly — a convenient choice that eliminates floating-point factors.

### 3.4 Conservation Laws

For an isolated gravitational system, three quantities are conserved:

**Total Energy:**
$$E = \underbrace{\frac{1}{2}\sum_{i=1}^{3} m_i |\mathbf{v}_i|^2}_{KE} - \underbrace{G\left(\frac{m_1 m_2}{r_{12}} + \frac{m_1 m_3}{r_{13}} + \frac{m_2 m_3}{r_{23}}\right)}_{-PE} = \text{const}$$

**Total Linear Momentum:**
$$\mathbf{P} = \sum_{i=1}^{3} m_i \mathbf{v}_i = \text{const} \quad \Rightarrow \quad \text{CM frame: } \mathbf{P} = \mathbf{0}$$

**Total Angular Momentum (2D, z-component):**
$$L = \sum_{i=1}^{3} m_i (x_i \dot{y}_i - y_i \dot{x}_i) = \text{const}$$

These are monitored throughout every simulation as numerical accuracy diagnostics.

### 3.5 Softening Parameter

To prevent numerical blow-up at close encounters, the force law is softened:

$$|\mathbf{r}|^3 \rightarrow \left(|\mathbf{r}|^2 + \varepsilon^2\right)^{3/2}$$

Default: `ε = 1e-4 AU`. This preserves long-range behavior while regularizing singularities.

### 3.6 The Lyapunov Exponent (Benettin Algorithm)

The maximal Lyapunov exponent is computed via the standard renormalization algorithm (Benettin et al. 1980):

1. Integrate reference trajectory **y(t)** and perturbed **ỹ(t) = y(0) + ε·δ̂**
2. At intervals τ: measure growth Λ_k = (1/τ) ln(|ỹ(τ) − y(τ)| / ε)
3. Renormalize: **ỹ ← y + ε (ỹ − y)/|ỹ − y|**
4. λ_L = lim_{N→∞} (1/N) Σ_k Λ_k

---

## 4. Installation

### Prerequisites

- Python ≥ 3.10
- conda (recommended) or pip

### Option A: conda (recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/three-body-problem.git
cd three-body-problem

# Create environment
conda env create -f environment.yml
conda activate three-body

# Install as editable package
pip install -e .
```

### Option B: pip

```bash
git clone https://github.com/yourusername/three-body-problem.git
cd three-body-problem

python -m venv .venv
source .venv/bin/activate      # Linux/Mac
# .venv\Scripts\activate       # Windows

pip install -r requirements.txt
pip install -e .
```

### Verify Installation

```bash
python -c "import src; print('✓ Installation OK')"
pytest tests/ -v --tb=short
```

---

## 5. Project Structure

```
three-body-problem/
│
├── 📄 README.md                    ← You are here
├── 📄 requirements.txt             ← Python dependencies
├── 📄 environment.yml              ← Conda environment
├── 📄 setup.py                     ← Package installation
│
├── 📁 src/                         ← Main source package
│   ├── __init__.py                 ← Public API exports
│   │
│   ├── 📁 core/                    ← Physics engine
│   │   ├── __init__.py
│   │   ├── equations.py            ← Gravitational EOMs, energy, angular momentum
│   │   ├── integrator.py           ← RK45, DOP853, Leapfrog, RK4 + SimulationResult
│   │   └── initial_conditions.py  ← 6 famous configurations + custom builder
│   │
│   ├── 📁 visualization/           ← Plotting and animation
│   │   ├── __init__.py
│   │   └── visualize.py            ← Trajectories, energy, phase space, animations
│   │
│   ├── 📁 analysis/                ← Physics analysis
│   │   ├── __init__.py
│   │   └── analysis.py             ← Lyapunov exponent, conservation errors, ejection detection
│   │
│   └── 📁 utils/                   ← Shared utilities
│       └── __init__.py
│
├── 📁 notebooks/                   ← Jupyter tutorials (run in order)
│   ├── 01_getting_started.ipynb    ← First simulation, state vector, integrators
│   ├── 02_famous_orbits.ipynb      ← All 6 configurations with analysis
│   ├── 03_chaos_and_lyapunov.ipynb ← Lyapunov exponent, chaos vs periodic
│   ├── 04_conservation_laws.ipynb  ← Energy, momentum, symplectic integrators
│   └── 05_research_connection.ipynb← Bridge to Kerr geodesics and BH research
│
├── 📁 scripts/                     ← Command-line tools
│   └── run_simulation.py           ← CLI entry point
│
├── 📁 tests/                       ← Unit test suite
│   ├── __init__.py
│   └── test_physics.py             ← Conservation laws, known solutions, integrator comparison
│
├── 📁 docs/                        ← Documentation
│   ├── 📁 theory/                  ← Mathematical derivations (Markdown)
│   │   ├── equations_of_motion.md
│   │   ├── integrators.md
│   │   └── lyapunov_exponent.md
│   └── 📁 figures/                 ← Auto-generated plots
│
└── 📁 data/
    ├── 📁 initial_conditions/      ← Saved IC JSON files
    └── 📁 outputs/                 ← Simulation outputs, animations
```

---

## 6. Quick Start

### Minimal Example (15 lines)

```python
import numpy as np
from src.core.initial_conditions import figure_eight, build_state_vector
from src.core.integrator import integrate_scipy
from src.core.equations import transform_to_cm_frame
from src.visualization.visualize import plot_trajectories
import matplotlib.pyplot as plt

# Load famous figure-eight initial conditions
ic = figure_eight()
masses = np.array([ic['m1'], ic['m2'], ic['m3']])

# Transform to center-of-mass frame
r1,r2,r3,v1,v2,v3 = transform_to_cm_frame(
    ic['r1'],ic['r2'],ic['r3'],
    ic['v1'],ic['v2'],ic['v3'], *masses)
ic.update({'r1':r1,'r2':r2,'r3':r3,'v1':v1,'v2':v2,'v3':v3})

# Integrate for 20 years
y0 = build_state_vector(ic)
t_eval = np.linspace(0, 20, 8000)
result = integrate_scipy(y0, masses, (0,20), t_eval, method='DOP853')

# Plot and print summary
print(result.summary())
fig = plot_trajectories(result, title=ic['name'])
plt.show()
```

### Expected Output

```
Three-Body Simulation
────────────────────────────────────────
  Integrator     : DOP853
  Duration       : 20.00 yr
  Steps          : 8000
  Max ΔE/E₀      : 2.34e-11
  Success        : True
  Masses (M☉)    : [1. 1. 1.]
```

---

## 7. Initial Condition Catalogue

```python
from src.core.initial_conditions import list_configurations, CATALOGUE
list_configurations()
```

### Built-in Configurations

| Key | Name | Type | Stable | Reference |
|-----|------|------|--------|-----------|
| `figure_eight` | Figure-Eight Choreography | Periodic choreography | Yes | Chenciner & Montgomery (2000) |
| `lagrange_triangle` | Lagrange Equilateral Triangle | Rigid rotation | Yes (equal mass) | Lagrange (1772) |
| `pythagorean` | Pythagorean Problem | Chaotic → ejection | No | Burrau (1913) |
| `solar_system` | Sun-Jupiter-Saturn | Hierarchical | Yes (long-term) | Real solar system |
| `figure_eight_perturbed` | Perturbed Figure-Eight | Chaotic | No | — |
| `sun_earth_moon` | Sun-Earth-Moon | Hierarchical | Yes | Newton (1687) |

### Creating Custom Initial Conditions

```python
from src.core.initial_conditions import custom, build_state_vector
import numpy as np

ic = custom(
    r1=[1.0,  0.0],   r2=[-1.0, 0.0],  r3=[0.0, 1.5],
    v1=[0.0,  0.8],   v2=[0.0, -0.8],  v3=[0.0, 0.0],
    m1=1.0, m2=1.0, m3=0.5,
    name="My Custom Config",
    t_end=30.0,
)
y0 = build_state_vector(ic)
```

---

## 8. Integrators

### Integrator Comparison

| Integrator | Type | Step | Energy drift | Best for |
|-----------|------|------|-------------|----------|
| `DOP853` | Adaptive Runge-Kutta 8(5,3) | Variable | ~1e-10 | Long accurate runs |
| `RK45` | Adaptive Runge-Kutta 4(5) | Variable | ~1e-8 | General use |
| `Leapfrog` | Symplectic Störmer-Verlet | Fixed | Oscillates, no drift | Very long integrations |
| `RK4` | Classic Runge-Kutta 4 | Fixed | Slow secular drift | Teaching, debugging |

### Usage

```python
from src.core.integrator import integrate_scipy, integrate_leapfrog, integrate_rk4
import numpy as np

# Adaptive (DOP853) — recommended default
result = integrate_scipy(y0, masses, t_span=(0,50),
                         t_eval=np.linspace(0,50,10000),
                         method='DOP853', rtol=1e-10, atol=1e-12)

# Symplectic leapfrog — best for long periodic orbit simulations
result = integrate_leapfrog(y0, masses, t_span=(0,200), dt=0.001)

# Classic RK4 — transparent, easy to understand
result = integrate_rk4(y0, masses, t_span=(0,10), dt=0.001)
```

### The `SimulationResult` Object

Every integrator returns a `SimulationResult` dataclass:

```python
result.t              # np.ndarray (N,)        — time array [yr]
result.positions      # np.ndarray (N, 3, 2)   — positions [AU]
result.velocities     # np.ndarray (N, 3, 2)   — velocities [AU/yr]
result.masses         # np.ndarray (3,)         — masses [M☉]
result.energy         # np.ndarray (N,)         — total energy
result.angular_momentum # np.ndarray (N,)       — total L
result.energy_error   # np.ndarray (N,)         — |ΔE/E₀|
result.max_energy_error # float                 — max energy error
result.integrator     # str                     — integrator name
result.n_steps        # int                     — number of stored steps
result.duration       # float                   — total time [yr]
result.summary()      # str                     — formatted summary
```

---

## 9. Analysis Tools

```python
from src.analysis.analysis import (
    estimate_lyapunov_exponent,   # maximal Lyapunov exponent
    energy_conservation_error,    # |ΔE/E₀| array
    angular_momentum_error,       # |ΔL/L₀| array
    separation_matrix,            # pairwise |r_i - r_j| array
    minimum_separation,           # closest approach (r, t, pair)
    detect_ejection,              # ejection event detection
    compute_orbital_period,       # period from zero-crossings
    kinetic_energy_partition,     # per-body KE fractions
    print_summary,                # formatted analysis output
)
```

### Lyapunov Exponent

```python
from src.core.initial_conditions import pythagorean
from src.analysis.analysis import estimate_lyapunov_exponent
from src.visualization.visualize import plot_lyapunov

ic = pythagorean()
times, lyapunov = estimate_lyapunov_exponent(
    ic,
    t_max=60.0,       # integration time [yr]
    dt=0.005,         # timestep [yr]
    epsilon=1e-6,     # perturbation size [AU]
    n_renorm=200,     # renormalization steps
)

print(f'λ_L = {lyapunov[-1]:.4f} yr⁻¹')
fig = plot_lyapunov(times, lyapunov)
```

### Close Encounter Detection

```python
from src.analysis.analysis import minimum_separation, detect_ejection

min_sep, t_min, pair = minimum_separation(result)
pairs = ['1-2', '1-3', '2-3']
print(f'Closest approach: {min_sep:.4f} AU between bodies {pairs[pair]} at t={t_min:.2f} yr')

ejection = detect_ejection(result, threshold_au=100.0)
if ejection:
    body_idx, t_eject = ejection
    print(f'Body {body_idx+1} ejected at t = {t_eject:.2f} yr')
```

---

## 10. Notebooks Guide

Work through the notebooks in order. Each builds on the previous.

| Notebook | Topic | Key concepts | Time |
|----------|-------|-------------|------|
| `01_getting_started.ipynb` | First simulation | State vector, ODE setup, integrator comparison | 20 min |
| `02_famous_orbits.ipynb` | Famous configurations | Stability, chaos, hierarchical systems | 30 min |
| `03_chaos_and_lyapunov.ipynb` | Chaos quantification | Lyapunov exponent, Benettin algorithm | 45 min |
| `04_conservation_laws.ipynb` | Conservation laws | Energy, angular momentum, virial theorem, symplectic | 30 min |
| `05_research_connection.ipynb` | GR bridge | Kerr geodesics, photon rings, γ vs λ_L | 40 min |

### Launch Notebooks

```bash
conda activate three-body
cd notebooks
jupyter lab
```

---

## 11. CLI Reference

### Run a simulation

```bash
cd three-body-problem

# Basic run — figure-eight, 20 years
python scripts/run_simulation.py --config figure_eight

# Pythagorean problem with DOP853, compute Lyapunov, save plots
python scripts/run_simulation.py \
    --config pythagorean \
    --integrator DOP853 \
    --duration 70 \
    --lyapunov \
    --save-plots

# Sun-Jupiter-Saturn for 100 years with animation
python scripts/run_simulation.py \
    --config solar_system \
    --duration 100 \
    --animate

# Perturbed figure-eight, symplectic integrator
python scripts/run_simulation.py \
    --config figure_eight_perturbed \
    --integrator Leapfrog \
    --dt 0.0005 \
    --duration 30

# List all configurations
python scripts/run_simulation.py --list
```

### CLI Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--config` | str | `figure_eight` | Configuration name |
| `--integrator` | str | `DOP853` | Integration method |
| `--duration` | float | from config | Simulation time [yr] |
| `--dt` | float | `0.001` | Fixed timestep (Leapfrog/RK4) [yr] |
| `--n-points` | int | `5000` | Output points (scipy integrators) |
| `--rtol` | float | `1e-10` | Relative tolerance |
| `--lyapunov` | flag | — | Compute Lyapunov exponent |
| `--animate` | flag | — | Save GIF animation |
| `--save-plots` | flag | — | Save PNG figures |
| `--no-plot` | flag | — | Suppress all plots |
| `--list` | flag | — | Show all configurations |

---

## 12. Running Tests

```bash
conda activate three-body

# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ -v --cov=src --cov-report=term-missing

# Run specific test class
pytest tests/test_physics.py::TestConservationLaws -v

# Run specific test
pytest tests/test_physics.py::TestKnownSolutions::test_figure_eight_returns_to_start -v
```

### Test Suite Structure

| Test Class | What it validates |
|-----------|------------------|
| `TestEquationsOfMotion` | Newton's 3rd law, momentum conservation, force symmetry |
| `TestCenterOfMass` | CM frame transformation, total momentum = 0 |
| `TestConservationLaws` | Energy/angular momentum conservation for all integrators |
| `TestKnownSolutions` | Lagrange period, figure-8 periodicity, two-body limit |
| `TestIntegratorComparison` | All integrators agree on short integrations |

### Expected Results

```
tests/test_physics.py::TestEquationsOfMotion::test_momentum_conservation PASSED
tests/test_physics.py::TestConservationLaws::test_energy_conservation[figure_eight-DOP853] PASSED
tests/test_physics.py::TestKnownSolutions::test_figure_eight_returns_to_start PASSED
... (all 14 tests pass)
```

---

## 13. Connection to Black Hole Research

This project was built specifically to develop intuition for Kerr black hole ray tracing. The connections are direct and mathematical.

### The Lyapunov Exponent Bridge

| Concept | Three-Body Problem | Kerr Black Hole |
|---------|-------------------|-----------------|
| Equation | `|δz(t)| ~ exp(λ_L t)` | `w_{n+1}/w_n = exp(-γ)` |
| Variable | λ_L [yr⁻¹] | γ (dimensionless) |
| Meaning | Orbital instability rate | Photon orbit instability |
| Value (example) | ~0.25 yr⁻¹ (Pythagorean) | π ≈ 3.14 (Schwarzschild) |
| Computed via | Benettin renormalization | Gralla & Lupsasca (2020) integral |
| Your result | Numerical λ_L vs config | Numerical γ(a, θ) vs Kerr parameters |

### Conservation Laws Bridge

| Newtonian | GR Analog | Conserved along |
|-----------|-----------|-----------------|
| Energy E = ½Σmv² - GPE | E = -p_μ ξ^μ (ξ = ∂_t) | Timelike/null geodesics |
| L = Σm(r × v) | L_z = p_μ ψ^μ (ψ = ∂_φ) | Timelike/null geodesics |
| No analog | Carter constant Q | Null geodesics in Kerr |
| Linear momentum P | CM frame condition | Fixed in observer frame |

### Integrator Bridge

The integrators you use here carry directly to your Kerr ray tracer:

```python
# Three-body (this project):
from scipy.integrate import solve_ivp
sol = solve_ivp(equations_of_motion, t_span, y0,
                method='DOP853', rtol=1e-10, events=[...])

# Kerr ray tracer (your research):
from scipy.integrate import solve_ivp
sol = solve_ivp(carter_geodesic_rhs, lam_span, state0,
                method='DOP853', rtol=1e-10, events=[horizon_event])
```

The **identical scipy call signature** — method, tolerances, events — carries over. Understanding what `rtol=1e-10` means here directly informs your choice of tolerances for the geodesic integrator.

### Skills Developed Here → Research Paper

| Skill practiced here | Applied in paper |
|---------------------|-----------------|
| ODE integration with solve_ivp | Kerr geodesic integration |
| Conservation law monitoring | Null condition + E, L_z, Q monitoring |
| Lyapunov exponent computation | γ(a,θ) Lyapunov exponent for Kerr |
| Energy drift comparison across integrators | Choosing DOP853 over RK4 for geodesics |
| Phase space analysis | Photon orbit classification (captured/escaped) |
| Close encounter detection | Horizon detection event in ODE solver |
| Ejection detection | Escape classification (ray reaches r → ∞) |

---

## 14. References

### Foundational Papers

| Reference | What it introduces |
|-----------|-------------------|
| Newton (1687), *Principia Mathematica* | Law of universal gravitation |
| Lagrange (1772), *Prix de l'Académie Royale* | Equilateral triangle exact solution |
| Poincaré (1890), *Acta Mathematica* | Proof of chaos, first use of phase space |
| Burrau (1913), *Vierteljahrsschr. Astron.* | First numerical 3-body integration (Pythagorean) |
| Szebehely & Peters (1967), *AJ 72, 876* | Pythagorean ejection analysis |
| Benettin et al. (1980), *Meccanica 15, 9* | Lyapunov exponent algorithm |
| Moore (1993), *Phys. Rev. Lett. 70, 3675* | Figure-eight numerical discovery |
| Chenciner & Montgomery (2000), *Ann. Math. 152* | Figure-eight existence proof |
| Šuvakov & Dmitrašinović (2013), *PRL 110* | 13 new periodic families |

### Black Hole Connection

| Reference | Relevance |
|-----------|-----------|
| Gralla, Holz & Wald (2019), arXiv:1910.12873 | Photon ring definition |
| Johnson et al. (2020), arXiv:1907.04329 | Universal photon ring signatures |
| Gralla & Lupsasca (2020), arXiv:1912.07586 | Lyapunov exponent γ formula for Kerr |
| EHT M87* Paper I (2019), arXiv:1906.11238 | First black hole image |

### Textbooks

| Book | Chapters |
|------|---------|
| Valtonen & Karttunen (2006), *The Three-Body Problem* | All |
| Heggie & Hut (2003), *The Gravitational Million-Body Problem* | Ch. 1–5 |
| Hairer, Lubich & Wanner (2006), *Geometric Numerical Integration* | Ch. 1–4 (symplectic integrators) |
| Schutz (2009), *A First Course in General Relativity* | Ch. 5–9 (bridge to GR) |

---

## 15. License

```
MIT License

Copyright (c) 2026 [Your Name]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
```

---

<div align="center">

**Built as preparation for Kerr black hole shadow ray-tracing research**
**Target: MEXT Scholarship → NAOJ / ICRR Tokyo**

*"The three-body problem is the simplest example of a system that defeats our ability to predict the future — and teaches us the mathematics of chaos that governs photon rings around black holes."*

</div>