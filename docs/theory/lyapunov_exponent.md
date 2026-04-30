# The Lyapunov Exponent: Three-Body Chaos & Black Hole Photon Rings

## 1. What is the Lyapunov Exponent?

The **maximal Lyapunov exponent** λ_L quantifies the rate of exponential divergence of nearby trajectories in phase space.

Given two trajectories starting at $\mathbf{y}_0$ and $\mathbf{y}_0 + \boldsymbol{\delta}_0$ (with $|\boldsymbol{\delta}_0| = \varepsilon \ll 1$):

$$|\boldsymbol{\delta}(t)| \approx \varepsilon \cdot e^{\lambda_L t}$$

- **λ_L > 0**: Chaotic orbit — nearby trajectories diverge exponentially
- **λ_L = 0**: Regular/periodic orbit — trajectories remain bounded
- **λ_L < 0**: Attracting orbit — trajectories converge (rare in Hamiltonian systems)

## 2. The Benettin Algorithm

Direct measurement of $e^{\lambda_L t}$ fails because once $|\boldsymbol{\delta}(t)|$ becomes large, the linear approximation breaks down. The solution: **periodic renormalization**.

**Algorithm (Benettin et al. 1980):**

Divide total time $T$ into $N$ intervals of length $\tau = T/N$.

For $k = 1, \ldots, N$:
1. Integrate reference **y** and perturbed **ỹ** for time $\tau$
2. Measure growth: 
$$\Lambda_k = \frac{1}{\tau} \ln\frac{|\tilde{\mathbf{y}}(k\tau) - \mathbf{y}(k\tau)|}{\varepsilon}$$
3. Renormalize:
$$\tilde{\mathbf{y}} \leftarrow \mathbf{y}(k\tau) + \varepsilon \cdot \frac{\tilde{\mathbf{y}} - \mathbf{y}}{|\tilde{\mathbf{y}} - \mathbf{y}|}$$

**Result:** 
$$\lambda_L = \lim_{N \to \infty} \frac{1}{N} \sum_{k=1}^N \Lambda_k$$

In practice, this converges after $\sim 50$–$200$ renormalization steps for three-body orbits.

## 3. Example Values

| Configuration | λ_L (yr⁻¹) | Type | e-folding time |
|--------------|------------|------|---------------|
| Figure-Eight (exact) | ~0.001 | Regular | ~1000 yr |
| Lagrange Triangle | ~0.000 | Periodic | ∞ |
| Figure-Eight (perturbed) | ~0.3–0.5 | Chaotic | 2–3 yr |
| Pythagorean | ~0.2–0.4 | Chaotic | 3–5 yr |
| Sun-Jupiter-Saturn | ~0.001 | Quasi-periodic | ~1000 yr |

## 4. The Black Hole Connection

This is the central insight connecting this project to your research:

**The Lyapunov exponent of the photon sphere** in Kerr spacetime is defined by the same equation — it measures how fast nearby photon orbits diverge.

### For Schwarzschild (non-rotating black hole)

Photons on nearly-circular orbits near the photon sphere ($r = 3M$) diverge as:

$$|\delta r(\lambda)| \sim |\delta r_0| \cdot e^{\gamma \cdot n_\text{orbits}}$$

where $\gamma$ is the orbital Lyapunov exponent. For Schwarzschild:

$$\boxed{\gamma_\text{Schwarzschild} = \pi \approx 3.14159}$$

This is a **universal constant** — independent of the photon's angular momentum!

### For Kerr (rotating black hole, spin $a$)

The Kerr Lyapunov exponent depends on spin $a$ and observer inclination $\theta$:

$$\gamma(a, \theta) = \frac{\pi}{\int_{r_-}^{r_+} \frac{dr}{\sqrt{\tilde{R}(r)}}} \times (\text{normalization})$$

where $\tilde{R}(r)$ is the radial potential evaluated on the critical orbit. The exact formula is given in Gralla & Lupsasca (2020b), arXiv:1912.07586.

**Your research contribution:** Numerically compute $\gamma(a, \theta)$ across the full Kerr parameter space and compare to the Gralla & Lupsasca analytical predictions.

## 5. Photon Ring Width Formula

Each successive photon ring sub-image is exponentially narrower:

$$w_{n+1} = w_n \cdot e^{-\gamma}$$

and exponentially fainter:

$$F_{n+1} = F_n \cdot e^{-\gamma}$$

**Numerical values:**

| BH | Spin $a/M$ | γ | $e^{-\gamma}$ | n=1 ring width |
|----|-----------|---|--------------|----------------|
| Schwarzschild | 0 | π ≈ 3.14 | 0.043 | 23× thinner than n=0 |
| Kerr (M87*) | ~0.9 | ~2.76 | 0.063 | 16× thinner |
| Extremal Kerr | 1.0 | ~2.2 | 0.11 | 9× thinner |

## 6. Side-by-Side Comparison

```
Three-Body Chaos              Black Hole Photon Sphere
─────────────────────────     ────────────────────────────────
Phase space: (r, v) ∈ ℝ¹²    Phase space: (r, θ, p_r, p_θ) ∈ ℝ⁴

Unstable orbit:               Unstable orbit:
  circular orbit r(t)=const     photon sphere r=3M (Schwarzschild)

Perturbation grows as:        Perturbation grows as:
  |δz(t)| ~ exp(λ_L t)         |δb| ~ exp(γ · n_orbits)

Benettin measurement:         Gralla & Lupsasca formula:
  λ_L from numerical sims       γ from analytical integral over Kerr

Physical consequence:         Physical consequence:
  chaotic orbit = unpredictable  photon ring n+1 is e^{-γ} narrower
  regular orbit = predictable    larger γ = harder to detect n≥1 rings
```

## 7. Why This Matters for ngEHT

The ngEHT (next-generation Event Horizon Telescope) aims to detect the n=1 photon ring of M87*. Whether this is feasible depends on γ:

- **Large γ** → ring is very thin → requires very high baseline sensitivity
- **Small γ** → ring is less suppressed → easier to detect

For M87* with spin $a \approx 0.9M$ and inclination $\theta \approx 17°$:
$$\gamma \approx 2.76, \quad e^{-\gamma} \approx 0.063$$

The n=1 ring is ~16× thinner than the direct image. ngEHT's target baseline sensitivity is designed around this number.

**Your contribution:** A systematic numerical study of $\gamma(a, \theta)$ across the full Kerr parameter space, validated against Gralla & Lupsasca, gives observational predictions for a range of BH spin and observer geometries — directly useful for ngEHT mission planning.

## 8. References

- Benettin et al. (1980), *Meccanica* **15**, 9 — Lyapunov algorithm
- Gralla & Lupsasca (2020a), arXiv:1912.07585 — photon ring theory
- Gralla & Lupsasca (2020b), arXiv:1912.07586 — Lyapunov exponent formula
- Johnson et al. (2020), arXiv:1907.04329 — universal photon ring signatures
- Cardoso et al. (2009), arXiv:0812.1806 — Lyapunov in BH physics context