# Numerical Integrators: Theory & Implementation

## 1. The ODE Problem

The three-body equations of motion form a system of 12 coupled first-order ODEs:

$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}), \qquad \mathbf{y}(0) = \mathbf{y}_0$$

where $\mathbf{y} \in \mathbb{R}^{12}$ is the state vector.

**Goal:** Given $\mathbf{y}_0$, find $\mathbf{y}(t)$ for $t \in [0, T]$.

## 2. Classic Runge-Kutta 4 (RK4)

Fixed-step, 4th-order method. Error per step: $O(h^5)$. Global error: $O(h^4)$.

$$\mathbf{k}_1 = \mathbf{f}(t_n, \mathbf{y}_n)$$
$$\mathbf{k}_2 = \mathbf{f}(t_n + h/2, \mathbf{y}_n + \frac{h}{2}\mathbf{k}_1)$$
$$\mathbf{k}_3 = \mathbf{f}(t_n + h/2, \mathbf{y}_n + \frac{h}{2}\mathbf{k}_2)$$
$$\mathbf{k}_4 = \mathbf{f}(t_n + h, \mathbf{y}_n + h\mathbf{k}_3)$$
$$\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{h}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)$$

**Pros:** Simple, well-understood, easy to implement
**Cons:** Fixed step size, energy drifts secularly for long integrations

## 3. Dormand-Prince 8(5,3) — DOP853

Adaptive step size, 8th-order method. The workhorse of scientific ODE solving.

- Computes 13 function evaluations per step
- Estimates local error by comparing 8th and 5th order solutions
- Automatically adjusts step size to maintain `rtol` and `atol`
- Global error: $O(h^8)$

**Step size control:**
$$h_\text{new} = h \cdot \min\left(f_\text{max}, \max\left(f_\text{min}, f_\text{safe} \cdot \left(\frac{\varepsilon_\text{tol}}{\varepsilon_\text{local}}\right)^{1/8}\right)\right)$$

**Pros:** Highly accurate, adaptive, handles stiff regions automatically
**Cons:** Not symplectic — energy can drift slightly over very long runs

**Tolerance guide for three-body problem:**

| Use case | rtol | atol |
|----------|------|------|
| Quick exploration | `1e-6` | `1e-8` |
| Publication figures | `1e-10` | `1e-12` |
| Lyapunov computation | `1e-12` | `1e-14` |

## 4. Symplectic Leapfrog (Störmer-Verlet)

Fixed-step, 2nd-order **symplectic** integrator. The key property: it exactly conserves a modified Hamiltonian $\tilde{H} = H + O(h^2)$, preventing secular energy drift.

**Algorithm (velocity Verlet form):**

$$\mathbf{a}(t) = \mathbf{F}(\mathbf{r}(t)) / m$$

$$\mathbf{r}(t+h) = \mathbf{r}(t) + \mathbf{v}(t)h + \frac{1}{2}\mathbf{a}(t)h^2$$

$$\mathbf{a}(t+h) = \mathbf{F}(\mathbf{r}(t+h)) / m$$

$$\mathbf{v}(t+h) = \mathbf{v}(t) + \frac{1}{2}\left[\mathbf{a}(t) + \mathbf{a}(t+h)\right]h$$

**Why symplectic matters:**

Non-symplectic methods (RK4, DOP853) allow energy to drift slowly — $\Delta E \propto t$ for long integrations. Symplectic methods bound energy error — $\Delta E$ oscillates but does not grow:

| Integrator | Energy behavior | Long-run stability |
|-----------|----------------|-------------------|
| RK4 (non-symplectic) | Secular drift $\propto t$ | Poor |
| DOP853 (adaptive, non-symplectic) | Very small drift | Good for moderate $T$ |
| Leapfrog (symplectic) | Oscillates, no secular drift | Excellent for large $T$ |

**Recommended:** Use Leapfrog for periodic orbits over many periods. Use DOP853 for general trajectories.

## 5. Adaptive vs. Fixed Step Size

**Fixed step:** Same $h$ at every step. Fast per step, but requires small $h$ everywhere (even where solution is smooth).

**Adaptive step:** Step size adjusted based on local error estimate. Large steps where solution is smooth, small steps near close encounters.

**For three-body problem:** Adaptive (DOP853) is usually better because close encounters require very small steps — a fixed-step integrator wastes computation on smooth regions.

## 6. Event Detection

`scipy.integrate.solve_ivp` supports stopping the integration when a function crosses zero — used for horizon and escape detection in the Kerr ray tracer.

**Horizon event (stops integration when r < r_+):**
```python
def horizon_event(t, y, *args):
    r = y[1]            # radial coordinate
    return r - r_plus   # zero when r = r_+

horizon_event.terminal  = True    # stop integration
horizon_event.direction = -1      # only when r decreasing
```

This exact pattern is used in your Kerr geodesic integrator.

## 7. Error Analysis and Validation

### Conservation Law Monitoring

At each stored timestep, compute:
$$\delta E(t) = \frac{|E(t) - E(0)|}{|E(0)|}$$

This is your primary validation diagnostic. If $\delta E$ grows, either:
1. Tolerances are too loose (tighten `rtol`/`atol`)
2. Step size is too large (reduce `dt` for fixed-step)
3. Close encounter is occurring (add softening or event detection)

### Cross-Validation

Compare results from two different integrators:
$$\|\mathbf{y}_\text{RK45}(T) - \mathbf{y}_\text{DOP853}(T)\| < \text{tolerance}$$

If they agree, you can trust the result. If they disagree, at least one is wrong — the higher-order method is usually more reliable.

## 8. Connection to Kerr Ray Tracer

The ODE solver call in your Kerr ray tracer is structurally identical:

```python
# Three-body (this project)
sol = solve_ivp(
    equations_of_motion,          # f(t, y) — Newtonian acceleration
    t_span=(0, T),
    y0=initial_state,             # [r1, r2, r3, v1, v2, v3]
    method='DOP853',
    rtol=1e-10, atol=1e-12,
    events=[horizon_event],
)

# Kerr ray tracer (your research)
sol = solve_ivp(
    carter_geodesic_rhs,          # f(λ, y) — Carter equations
    t_span=(0, lam_max),
    y0=photon_state,              # [t, r, θ, φ, p_t, p_r, p_θ, p_φ]
    method='DOP853',
    rtol=1e-10, atol=1e-12,
    events=[horizon_event, escape_event],
)
```

The structure is the same. The `equations_of_motion` function becomes `carter_geodesic_rhs` — but the numerical machinery is identical.