"""
integrator.py
=============
Numerical integration of the three-body ODEs.

Available Integrators
---------------------
1. RK45   — Runge-Kutta 4(5), adaptive step, scipy (default)
2. DOP853 — Dormand-Prince 8(5,3), high-accuracy adaptive, scipy
3. Leapfrog — Symplectic integrator, fixed step, energy-conserving
4. RK4    — Classic 4th-order Runge-Kutta, fixed step

Integrator Choice Guide
-----------------------
| Integrator | Energy conservation | Speed | Best for             |
|------------|--------------------|----|----------------------|
| RK45       | Good (adaptive)    | Fast  | General exploration  |
| DOP853     | Excellent          | Medium| Long integrations    |
| Leapfrog   | Excellent (symplectic)| Fast| Periodic orbits      |
| RK4        | Good (fixed step)  | Fast  | Teaching/debugging   |

For research-grade simulations of periodic orbits (figure-8, Lagrange),
use the Leapfrog (Störmer-Verlet) integrator — it exactly conserves a
modified Hamiltonian, preventing secular energy drift over long times.

References
----------
- Hairer, Lubich & Wanner (2006), "Geometric Numerical Integration"
- Leimkuhler & Reich (2004), "Simulating Hamiltonian Dynamics"
- Dormand & Prince (1980), J. Comput. Appl. Math. 6, 19-26
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
from typing import Optional, Callable
from .equations import equations_of_motion, total_energy, total_angular_momentum, G_UNITS


@dataclass
class SimulationResult:
    """
    Container for a completed three-body simulation.

    Attributes
    ----------
    t : np.ndarray, shape (N,)
        Time array [yr]
    positions : np.ndarray, shape (N, 3, 2)
        Position of each body at each timestep [AU]
        positions[:, 0, :] = body 1 positions (x,y)
        positions[:, 1, :] = body 2 positions (x,y)
        positions[:, 2, :] = body 3 positions (x,y)
    velocities : np.ndarray, shape (N, 3, 2)
        Velocity of each body [AU/yr]
    masses : np.ndarray, shape (3,)
        Masses [M_sun]
    energy : np.ndarray, shape (N,)
        Total energy at each timestep [M_sun AU^2/yr^2]
    angular_momentum : np.ndarray, shape (N,)
        Total angular momentum [M_sun AU^2/yr]
    energy_error : np.ndarray, shape (N,)
        Relative energy error |(E - E0) / E0|
    integrator : str
        Name of integrator used
    success : bool
        Whether integration completed without errors
    message : str
        Status message from integrator
    """
    t: np.ndarray
    positions: np.ndarray
    velocities: np.ndarray
    masses: np.ndarray
    energy: np.ndarray
    angular_momentum: np.ndarray
    energy_error: np.ndarray
    integrator: str = "RK45"
    success: bool = True
    message: str = ""

    @property
    def n_steps(self) -> int:
        return len(self.t)

    @property
    def duration(self) -> float:
        return self.t[-1] - self.t[0]

    @property
    def max_energy_error(self) -> float:
        return float(np.max(self.energy_error))

    def summary(self) -> str:
        return (
            f"Three-Body Simulation\n"
            f"{'─'*40}\n"
            f"  Integrator     : {self.integrator}\n"
            f"  Duration       : {self.duration:.2f} yr\n"
            f"  Steps          : {self.n_steps}\n"
            f"  Max ΔE/E₀      : {self.max_energy_error:.2e}\n"
            f"  Success        : {self.success}\n"
            f"  Masses (M☉)    : {self.masses}\n"
        )


def integrate_scipy(
    y0: np.ndarray,
    masses: np.ndarray,
    t_span: tuple,
    t_eval: np.ndarray,
    method: str = "RK45",
    rtol: float = 1e-10,
    atol: float = 1e-12,
    G: float = G_UNITS,
) -> SimulationResult:
    """
    Integrate using scipy.integrate.solve_ivp (RK45 or DOP853).

    Parameters
    ----------
    y0 : np.ndarray, shape (12,)
        Initial state vector [x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3]
    masses : np.ndarray, shape (3,)
        [m1, m2, m3] in solar masses
    t_span : tuple (t_start, t_end)
        Integration interval [yr]
    t_eval : np.ndarray
        Times at which to record solution
    method : str
        'RK45' or 'DOP853'
    rtol, atol : float
        Relative and absolute tolerances
    G : float
        Gravitational constant

    Returns
    -------
    SimulationResult
    """
    sol = solve_ivp(
        equations_of_motion,
        t_span,
        y0,
        method=method,
        t_eval=t_eval,
        args=(masses, G),
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    return _build_result(sol.t, sol.y.T, masses, method, sol.success, sol.message)


def integrate_leapfrog(
    y0: np.ndarray,
    masses: np.ndarray,
    t_span: tuple,
    dt: float = 0.001,
    store_every: int = 10,
    G: float = G_UNITS,
) -> SimulationResult:
    """
    Symplectic Störmer-Verlet (Leapfrog) integrator.

    Algorithm (velocity Verlet form):
        a(t)       = F(r(t)) / m
        r(t+dt)    = r(t) + v(t)*dt + 0.5*a(t)*dt²
        a(t+dt)    = F(r(t+dt)) / m
        v(t+dt)    = v(t) + 0.5*(a(t) + a(t+dt))*dt

    This is a symplectic integrator — it exactly conserves a modified
    Hamiltonian H̃ = H + O(dt²), preventing secular energy drift.
    Critical for long-duration periodic orbit simulations.

    Parameters
    ----------
    y0 : np.ndarray, shape (12,)
        Initial state vector
    masses : np.ndarray, shape (3,)
        Body masses [M_sun]
    t_span : tuple
        (t_start, t_end) [yr]
    dt : float
        Fixed timestep [yr]. Default 0.001 yr ≈ 9 hours
    store_every : int
        Store state every N steps (reduces memory for long runs)
    G : float
        Gravitational constant

    Returns
    -------
    SimulationResult
    """
    from .equations import acceleration as _accel

    t_start, t_end = t_span
    n_steps = int((t_end - t_start) / dt)
    m1, m2, m3 = masses

    # Unpack initial state
    r1 = y0[0:2].copy()
    r2 = y0[2:4].copy()
    r3 = y0[4:6].copy()
    v1 = y0[6:8].copy()
    v2 = y0[8:10].copy()
    v3 = y0[10:12].copy()

    # Storage
    store_n = n_steps // store_every + 1
    t_arr   = np.zeros(store_n)
    y_arr   = np.zeros((store_n, 12))
    t_arr[0] = t_start
    y_arr[0] = y0.copy()
    store_idx = 1

    # Initial acceleration
    a1, a2, a3 = _accel(r1, r2, r3, m1, m2, m3, G=G)

    t = t_start
    for step in range(1, n_steps + 1):
        # Update positions
        r1 += v1 * dt + 0.5 * a1 * dt**2
        r2 += v2 * dt + 0.5 * a2 * dt**2
        r3 += v3 * dt + 0.5 * a3 * dt**2

        # New accelerations
        a1_new, a2_new, a3_new = _accel(r1, r2, r3, m1, m2, m3, G=G)

        # Update velocities
        v1 += 0.5 * (a1 + a1_new) * dt
        v2 += 0.5 * (a2 + a2_new) * dt
        v3 += 0.5 * (a3 + a3_new) * dt

        a1, a2, a3 = a1_new, a2_new, a3_new
        t += dt

        if step % store_every == 0:
            t_arr[store_idx] = t
            y_arr[store_idx] = np.concatenate([r1, r2, r3, v1, v2, v3])
            store_idx += 1

    # Trim unused storage
    t_arr = t_arr[:store_idx]
    y_arr = y_arr[:store_idx]

    return _build_result(t_arr, y_arr, masses, "Leapfrog", True, "OK")


def integrate_rk4(
    y0: np.ndarray,
    masses: np.ndarray,
    t_span: tuple,
    dt: float = 0.001,
    store_every: int = 10,
    G: float = G_UNITS,
) -> SimulationResult:
    """
    Classic 4th-order Runge-Kutta integrator (fixed step).

    The RK4 method:
        k1 = f(t,       y)
        k2 = f(t+dt/2,  y + dt/2 * k1)
        k3 = f(t+dt/2,  y + dt/2 * k2)
        k4 = f(t+dt,    y + dt   * k3)
        y(t+dt) = y(t) + dt/6 * (k1 + 2k2 + 2k3 + k4)

    Error: O(dt^5) per step, O(dt^4) globally.
    Not symplectic — energy drifts slowly over long integrations.
    Good for short-duration simulations and teaching.
    """
    t_start, t_end = t_span
    n_steps = int((t_end - t_start) / dt)

    store_n = n_steps // store_every + 1
    t_arr = np.zeros(store_n)
    y_arr = np.zeros((store_n, 12))
    t_arr[0] = t_start
    y_arr[0] = y0.copy()
    store_idx = 1

    y = y0.copy()
    t = t_start

    for step in range(1, n_steps + 1):
        k1 = equations_of_motion(t,          y,              masses, G)
        k2 = equations_of_motion(t + dt/2,   y + dt/2 * k1, masses, G)
        k3 = equations_of_motion(t + dt/2,   y + dt/2 * k2, masses, G)
        k4 = equations_of_motion(t + dt,     y + dt   * k3, masses, G)

        y = y + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t += dt

        if step % store_every == 0:
            t_arr[store_idx] = t
            y_arr[store_idx] = y.copy()
            store_idx += 1

    t_arr = t_arr[:store_idx]
    y_arr = y_arr[:store_idx]

    return _build_result(t_arr, y_arr, masses, "RK4", True, "OK")


def _build_result(t, y, masses, method, success, message):
    """Build SimulationResult from raw integration output."""
    m1, m2, m3 = masses
    N = len(t)

    positions  = np.zeros((N, 3, 2))
    velocities = np.zeros((N, 3, 2))
    energy     = np.zeros(N)
    ang_mom    = np.zeros(N)

    for i in range(N):
        positions[i, 0] = y[i, 0:2]
        positions[i, 1] = y[i, 2:4]
        positions[i, 2] = y[i, 4:6]
        velocities[i, 0] = y[i, 6:8]
        velocities[i, 1] = y[i, 8:10]
        velocities[i, 2] = y[i, 10:12]

        energy[i] = total_energy(
            positions[i,0], positions[i,1], positions[i,2],
            velocities[i,0], velocities[i,1], velocities[i,2],
            m1, m2, m3,
        )
        ang_mom[i] = total_angular_momentum(
            positions[i,0], positions[i,1], positions[i,2],
            velocities[i,0], velocities[i,1], velocities[i,2],
            m1, m2, m3,
        )

    E0 = energy[0]
    energy_error = np.abs((energy - E0) / (E0 + 1e-30))

    return SimulationResult(
        t=t, positions=positions, velocities=velocities,
        masses=np.array(masses), energy=energy,
        angular_momentum=ang_mom, energy_error=energy_error,
        integrator=method, success=success, message=str(message),
    )