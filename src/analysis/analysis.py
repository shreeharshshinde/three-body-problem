"""
analysis.py
===========
Physics analysis tools for three-body simulation results.

Modules
-------
1. Conservation law monitoring  (energy, momentum, angular momentum)
2. Lyapunov exponent estimation (chaos quantification)
3. Orbital element extraction   (period, eccentricity, semi-major axis)
4. Close encounter detection    (minimum separation tracking)
5. Escape/ejection detection    (energy partition, Hill radius)

The Lyapunov Exponent Connection to Your Research
--------------------------------------------------
The Lyapunov exponent λ_L quantifies how fast nearby trajectories
diverge in phase space:

    |δz(t)| ~ |δz(0)| * e^(λ_L * t)

For chaotic three-body orbits:  λ_L > 0 (trajectories diverge)
For stable periodic orbits:     λ_L = 0 (trajectories stay close)

In your black hole research, the SAME mathematical object governs
photon ring widths:

    w_{n+1} / w_n = e^{-γ}

where γ = πγ_Kerr is the orbital Lyapunov exponent of the photon
sphere. Larger γ → more unstable photon orbits → narrower rings.

Computing λ_L for the three-body problem gives you direct intuition
for what γ means in the Kerr spacetime context.

References
----------
- Benettin et al. (1980), Meccanica 15, 9  (Lyapunov exponent algorithm)
- Heggie & Hut (2003), "The Gravitational Million-Body Problem"
- Valtonen & Karttunen (2006), "The Three-Body Problem"
"""

import numpy as np
from typing import Tuple, List, Optional
from .integrator import SimulationResult
from .equations import equations_of_motion, G_UNITS


def energy_conservation_error(result: SimulationResult) -> np.ndarray:
    """
    Compute relative energy error over time.

    ΔE/E₀ = |E(t) - E(0)| / |E(0)|

    A well-behaved integration should have ΔE/E₀ < 1e-6.
    The Leapfrog integrator keeps this bounded; RK45 can show
    slow drift for very long integrations.

    Returns
    -------
    rel_error : np.ndarray, shape (N,)
        Relative energy error at each timestep
    """
    return result.energy_error


def angular_momentum_error(result: SimulationResult) -> np.ndarray:
    """
    Relative angular momentum conservation error.

    ΔL/L₀ = |L(t) - L(0)| / |L(0)|

    Angular momentum is exactly conserved in closed gravitational
    systems. Any drift indicates numerical error.
    """
    L0 = result.angular_momentum[0]
    return np.abs((result.angular_momentum - L0) / (np.abs(L0) + 1e-30))


def separation_matrix(result: SimulationResult) -> np.ndarray:
    """
    Compute pairwise separations |r_i - r_j| at all timesteps.

    Returns
    -------
    sep : np.ndarray, shape (N, 3)
        sep[:, 0] = |r1 - r2|  (body 1–2 separation)
        sep[:, 1] = |r1 - r3|  (body 1–3 separation)
        sep[:, 2] = |r2 - r3|  (body 2–3 separation)
    """
    pos = result.positions
    sep = np.zeros((result.n_steps, 3))
    sep[:, 0] = np.linalg.norm(pos[:, 0] - pos[:, 1], axis=1)
    sep[:, 1] = np.linalg.norm(pos[:, 0] - pos[:, 2], axis=1)
    sep[:, 2] = np.linalg.norm(pos[:, 1] - pos[:, 2], axis=1)
    return sep


def minimum_separation(result: SimulationResult) -> Tuple[float, float, int]:
    """
    Find the closest approach between any two bodies.

    Returns
    -------
    min_sep : float
        Minimum separation [AU]
    time : float
        Time of closest approach [yr]
    pair : int
        Which pair: 0=(1,2), 1=(1,3), 2=(2,3)
    """
    seps = separation_matrix(result)
    min_idx = np.unravel_index(np.argmin(seps), seps.shape)
    t_idx, pair = min_idx
    return float(seps[t_idx, pair]), float(result.t[t_idx]), pair


def detect_ejection(
    result: SimulationResult,
    threshold_au: float = 100.0,
) -> Optional[Tuple[int, float]]:
    """
    Detect if any body escapes the system.

    A body is considered ejected when its distance from the center of
    mass exceeds threshold_au. This is a simplified criterion; a more
    rigorous one uses the Hill radius and energy partition.

    Returns
    -------
    (body_index, time) if ejection detected, else None
        body_index : 0, 1, or 2
        time : time of ejection [yr]
    """
    masses = result.masses
    M = masses.sum()
    pos = result.positions

    for i, t in enumerate(result.t):
        r_cm = (masses[0]*pos[i,0] + masses[1]*pos[i,1] + masses[2]*pos[i,2]) / M
        for body in range(3):
            dist = np.linalg.norm(pos[i, body] - r_cm)
            if dist > threshold_au:
                return (body, float(t))
    return None


def estimate_lyapunov_exponent(
    ic: dict,
    t_max: float = 50.0,
    dt: float = 0.01,
    epsilon: float = 1e-6,
    n_renorm: int = 100,
    G: float = G_UNITS,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Estimate the maximal Lyapunov exponent using the standard algorithm.

    Algorithm (Benettin et al. 1980)
    ---------------------------------
    1. Integrate reference trajectory y(t) from y0
    2. Integrate perturbed trajectory ỹ(t) from y0 + ε*δ
       where δ is a random unit vector in phase space
    3. At regular intervals τ:
       a. Compute growth: Λ_k = (1/τ) * ln(|ỹ - y| / ε)
       b. Rescale: reset ỹ = y + ε * (ỹ - y) / |ỹ - y|
    4. The Lyapunov exponent: λ_L = lim_{N→∞} (1/N) Σ_k Λ_k

    Connection to your BH research:
    The Lyapunov exponent γ of the photon sphere is defined identically —
    it measures how fast nearby photon orbits diverge in phase space.
    For Schwarzschild: γ = π (universal).
    Your code measures γ(a, θ) for Kerr numerically.

    Parameters
    ----------
    ic : dict
        Initial condition dictionary (from initial_conditions.py)
    t_max : float
        Total integration time [yr]
    dt : float
        Integration timestep [yr]
    epsilon : float
        Initial perturbation size [AU]
    n_renorm : int
        Renormalization steps (more = smoother estimate)
    G : float
        Gravitational constant

    Returns
    -------
    times : np.ndarray
        Times at renormalization steps [yr]
    lyapunov : np.ndarray
        Running estimate of λ_L at each step [1/yr]
    """
    from .integrator import integrate_rk4
    from .initial_conditions import build_state_vector

    masses = np.array([ic['m1'], ic['m2'], ic['m3']])
    y0 = build_state_vector(ic)

    # Random perturbation in phase space
    rng = np.random.default_rng(0)
    delta = rng.standard_normal(12)
    delta /= np.linalg.norm(delta)
    y0_pert = y0 + epsilon * delta

    tau = t_max / n_renorm   # time between renormalizations
    lyapunov_sum = 0.0
    times = []
    lyapunov_running = []

    t_current = 0.0
    y_ref  = y0.copy()
    y_pert = y0_pert.copy()

    for k in range(n_renorm):
        # Integrate both trajectories for time tau
        res_ref  = integrate_rk4(y_ref,  masses, (t_current, t_current + tau), dt=dt, store_every=1, G=G)
        res_pert = integrate_rk4(y_pert, masses, (t_current, t_current + tau), dt=dt, store_every=1, G=G)

        y_ref  = np.concatenate([res_ref.positions[-1].flatten(),  res_ref.velocities[-1].flatten()])
        y_pert = np.concatenate([res_pert.positions[-1].flatten(), res_pert.velocities[-1].flatten()])

        # Separation in phase space
        delta_new = y_pert - y_ref
        d_new = np.linalg.norm(delta_new)

        if d_new > 0:
            # Growth factor
            lyapunov_sum += np.log(d_new / epsilon)

            # Renormalize
            y_pert = y_ref + epsilon * delta_new / d_new

        t_current += tau
        times.append(t_current)
        lyapunov_running.append(lyapunov_sum / (t_current))

    return np.array(times), np.array(lyapunov_running)


def compute_orbital_period(
    result: SimulationResult,
    body: int = 0,
    coordinate: int = 0,
) -> Optional[float]:
    """
    Estimate orbital period from zero-crossing analysis.

    Finds successive zero-crossings of x(t) (or y(t)) for a given body,
    then averages the half-period intervals.

    Returns
    -------
    period : float or None
        Estimated period [yr], or None if not enough crossings found
    """
    signal = result.positions[:, body, coordinate]
    t = result.t

    # Find zero crossings
    crossings = []
    for i in range(1, len(signal)):
        if signal[i-1] * signal[i] < 0:
            # Linear interpolation for crossing time
            t_cross = t[i-1] - signal[i-1] * (t[i] - t[i-1]) / (signal[i] - signal[i-1])
            crossings.append(t_cross)

    if len(crossings) < 2:
        return None

    # Period from successive same-direction crossings
    half_periods = np.diff(crossings)
    period = 2 * np.mean(half_periods)
    return float(period)


def compute_specific_angular_momentum(result: SimulationResult, body: int) -> np.ndarray:
    """
    Compute specific angular momentum for a single body: l = r × v (z-component in 2D).

    Parameters
    ----------
    body : int
        Body index (0, 1, or 2)

    Returns
    -------
    l : np.ndarray, shape (N,)
        Specific angular momentum [AU^2/yr]
    """
    pos = result.positions[:, body, :]
    vel = result.velocities[:, body, :]
    return pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0]


def kinetic_energy_partition(result: SimulationResult) -> np.ndarray:
    """
    Compute kinetic energy fraction for each body.

    Returns
    -------
    ke_fraction : np.ndarray, shape (N, 3)
        KE of each body as fraction of total KE
    """
    masses = result.masses
    vel = result.velocities
    ke = np.zeros((result.n_steps, 3))

    for i in range(3):
        ke[:, i] = 0.5 * masses[i] * np.sum(vel[:, i, :]**2, axis=1)

    ke_total = ke.sum(axis=1, keepdims=True)
    return ke / (ke_total + 1e-30)


def print_summary(result: SimulationResult) -> None:
    """Print a comprehensive analysis summary."""
    print(result.summary())
    print(f"  Max ΔE/E₀      : {result.max_energy_error:.2e}")
    print(f"  Max ΔL/L₀      : {angular_momentum_error(result).max():.2e}")

    seps = separation_matrix(result)
    min_sep, t_min, pair = minimum_separation(result)
    pairs = ['1-2', '1-3', '2-3']
    print(f"  Min separation : {min_sep:.4f} AU (bodies {pairs[pair]} at t={t_min:.2f} yr)")

    ejection = detect_ejection(result)
    if ejection:
        print(f"  EJECTION: body {ejection[0]+1} at t = {ejection[1]:.2f} yr")
    else:
        print(f"  Ejection       : None detected")
    print()