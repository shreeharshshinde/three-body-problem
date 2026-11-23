"""
equations.py
============
Newton's gravitational equations of motion for the three-body problem.

Mathematical Foundation
-----------------------
The gravitational force on body i due to body j is:

    F_ij = G * m_i * m_j * (r_j - r_i) / |r_j - r_i|^3

Summing over all other bodies gives the acceleration of body i:

    a_i = d²r_i/dt² = G * Σ_{j≠i} m_j * (r_j - r_i) / |r_j - r_i|^3

For three bodies (i = 1, 2, 3), the full system is:

    d²r_1/dt² = G*m_2*(r_2-r_1)/|r_2-r_1|^3 + G*m_3*(r_3-r_1)/|r_3-r_1|^3
    d²r_2/dt² = G*m_1*(r_1-r_2)/|r_1-r_2|^3 + G*m_3*(r_3-r_2)/|r_3-r_2|^3
    d²r_3/dt² = G*m_1*(r_1-r_3)/|r_1-r_3|^3 + G*m_2*(r_2-r_3)/|r_2-r_3|^3

This gives 12 coupled ODEs in 2D (or 18 in 3D):
  - 6 (or 9) position components
  - 6 (or 9) velocity components

State vector (2D, n=12):
    y = [x1,y1, x2,y2, x3,y3, vx1,vy1, vx2,vy2, vx3,vy3]

State vector (3D, n=18):
    y = [x1,y1,z1, x2,y2,z2, x3,y3,z3, vx1,vy1,vz1, vx2,vy2,vz2, vx3,vy3,vz3]

References
----------
- Valtonen & Karttunen (2006), "The Three-Body Problem", Cambridge
- Musielak & Quarles (2014), Rep. Prog. Phys. 77, 065901
- Burrau (1913), first numerical integration of three-body problem
"""

import numpy as np
from typing import Tuple

# Gravitational constant in units where:
#   mass  → solar masses (M_sun)
#   length → AU (astronomical units)
#   time  → years
# G = 4π² AU³ / (M_sun · yr²)
G_UNITS = 4.0 * np.pi**2   # AU^3 / (M_sun * yr^2)

# Softening parameter to avoid singularities at close approaches
SOFTENING = 1e-4  # AU


def acceleration(
    r1: np.ndarray,
    r2: np.ndarray,
    r3: np.ndarray,
    m1: float,
    m2: float,
    m3: float,
    G: float = G_UNITS,
    softening: float = SOFTENING,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute gravitational accelerations for all three bodies.

    Parameters
    ----------
    r1, r2, r3 : np.ndarray, shape (2,) or (3,)
        Position vectors of each body [AU]
    m1, m2, m3 : float
        Masses of each body [M_sun]
    G : float
        Gravitational constant [AU^3 / (M_sun * yr^2)]
    softening : float
        Softening length to regularize close encounters [AU]

    Returns
    -------
    a1, a2, a3 : np.ndarray
        Acceleration vectors [AU / yr^2]

    Notes
    -----
    The softened force law replaces |r|^3 with (|r|^2 + ε^2)^(3/2),
    preventing numerical blow-up at close encounters while preserving
    long-range behavior. This is standard in N-body simulations.
    Works for both 2D and 3D position vectors — np.dot is dimension-agnostic.
    """
    r12 = r2 - r1
    r13 = r3 - r1
    r23 = r3 - r2

    d12 = (np.dot(r12, r12) + softening**2) ** 1.5
    d13 = (np.dot(r13, r13) + softening**2) ** 1.5
    d23 = (np.dot(r23, r23) + softening**2) ** 1.5

    a1 = G * (m2 * r12 / d12 + m3 * r13 / d13)
    a2 = G * (m1 * (-r12) / d12 + m3 * r23 / d23)
    a3 = G * (m1 * (-r13) / d13 + m2 * (-r23) / d23)

    return a1, a2, a3


def equations_of_motion(t: float, y: np.ndarray, masses: np.ndarray, G: float = G_UNITS) -> np.ndarray:
    """
    Full system of ODEs for scipy.integrate.solve_ivp.

    Supports both 2D (12-component) and 3D (18-component) state vectors.
    Dimensionality is inferred automatically from len(y).

    State vector layout (2D, n=12):
        y[0:2]   = r1 (x1, y1)
        y[2:4]   = r2 (x2, y2)
        y[4:6]   = r3 (x3, y3)
        y[6:8]   = v1 (vx1, vy1)
        y[8:10]  = v2 (vx2, vy2)
        y[10:12] = v3 (vx3, vy3)

    State vector layout (3D, n=18):
        y[0:3]   = r1 (x1, y1, z1)
        y[3:6]   = r2 (x2, y2, z2)
        y[6:9]   = r3 (x3, y3, z3)
        y[9:12]  = v1 (vx1, vy1, vz1)
        y[12:15] = v2 (vx2, vy2, vz2)
        y[15:18] = v3 (vx3, vy3, vz3)

    Parameters
    ----------
    t : float
        Current time [yr] (required by solve_ivp, not used explicitly)
    y : np.ndarray, shape (12,) or (18,)
        State vector
    masses : np.ndarray, shape (3,)
        [m1, m2, m3] in solar masses
    G : float
        Gravitational constant

    Returns
    -------
    dydt : np.ndarray, same shape as y
        Time derivative of state vector
    """
    n = len(y)
    if n == 12:
        d = 2
    elif n == 18:
        d = 3
    else:
        raise ValueError(f"State vector must have length 12 (2D) or 18 (3D), got {n}")

    r1, r2, r3 = y[0:d], y[d:2*d], y[2*d:3*d]
    v1, v2, v3 = y[3*d:4*d], y[4*d:5*d], y[5*d:6*d]
    m1, m2, m3 = masses

    a1, a2, a3 = acceleration(r1, r2, r3, m1, m2, m3, G=G)

    return np.concatenate([v1, v2, v3, a1, a2, a3])


def total_energy(
    r1: np.ndarray, r2: np.ndarray, r3: np.ndarray,
    v1: np.ndarray, v2: np.ndarray, v3: np.ndarray,
    m1: float, m2: float, m3: float,
    G: float = G_UNITS,
) -> float:
    """
    Compute total mechanical energy E = KE + PE.

    E = (1/2)(m1|v1|² + m2|v2|² + m3|v3|²)
      - G(m1*m2/r12 + m1*m3/r13 + m2*m3/r23)

    Works for both 2D and 3D vectors. np.dot and np.linalg.norm are
    dimension-agnostic.

    Returns
    -------
    E : float
        Total energy [M_sun * AU^2 / yr^2]
    """
    KE = 0.5 * (m1 * np.dot(v1, v1) + m2 * np.dot(v2, v2) + m3 * np.dot(v3, v3))

    r12 = np.linalg.norm(r2 - r1)
    r13 = np.linalg.norm(r3 - r1)
    r23 = np.linalg.norm(r3 - r2)

    PE = -G * (m1 * m2 / r12 + m1 * m3 / r13 + m2 * m3 / r23)

    return KE + PE


def total_angular_momentum(
    r1: np.ndarray, r2: np.ndarray, r3: np.ndarray,
    v1: np.ndarray, v2: np.ndarray, v3: np.ndarray,
    m1: float, m2: float, m3: float,
) -> np.ndarray:
    """
    Compute total angular momentum L = Σ m_i (r_i × v_i).

    Supports both 2D and 3D position/velocity vectors.

    - 2D input: returns a 1-element array [Lz] (z-component of the cross product)
    - 3D input: returns a 3-element array [Lx, Ly, Lz]

    The return value is always a NumPy array. For callers that need a scalar,
    use ``angular_momentum_magnitude()`` instead.

    Returns
    -------
    L : np.ndarray, shape (1,) for 2D or shape (3,) for 3D
        Total angular momentum [M_sun * AU^2 / yr]
    """
    dim = len(r1)
    if dim == 2:
        def cross(r, v):
            return np.array([r[0] * v[1] - r[1] * v[0]])
    else:
        def cross(r, v):
            return np.cross(r, v)

    return m1 * cross(r1, v1) + m2 * cross(r2, v2) + m3 * cross(r3, v3)


def angular_momentum_magnitude(
    r1: np.ndarray, r2: np.ndarray, r3: np.ndarray,
    v1: np.ndarray, v2: np.ndarray, v3: np.ndarray,
    m1: float, m2: float, m3: float,
) -> float:
    """
    Scalar magnitude of the total angular momentum vector.

    Works for both 2D and 3D. Used by the integrators to store a
    single scalar conservation diagnostic in SimulationResult.

    Returns
    -------
    |L| : float  [M_sun * AU^2 / yr]
    """
    L = total_angular_momentum(r1, r2, r3, v1, v2, v3, m1, m2, m3)
    return float(np.linalg.norm(L))


def center_of_mass(
    r1: np.ndarray, r2: np.ndarray, r3: np.ndarray,
    m1: float, m2: float, m3: float,
) -> np.ndarray:
    """
    Compute center of mass position.

    R_cm = (m1*r1 + m2*r2 + m3*r3) / (m1 + m2 + m3)

    Works for both 2D and 3D vectors.
    """
    M = m1 + m2 + m3
    return (m1 * r1 + m2 * r2 + m3 * r3) / M


def transform_to_cm_frame(
    r1, r2, r3, v1, v2, v3, m1, m2, m3
):
    """
    Shift positions and velocities to the center-of-mass frame.

    In the CM frame: Σ m_i r_i = 0, Σ m_i v_i = 0
    This eliminates linear drift of the system.
    Works for both 2D and 3D vectors.
    """
    M = m1 + m2 + m3
    r_cm = center_of_mass(r1, r2, r3, m1, m2, m3)
    v_cm = (m1 * v1 + m2 * v2 + m3 * v3) / M

    return (
        r1 - r_cm, r2 - r_cm, r3 - r_cm,
        v1 - v_cm, v2 - v_cm, v3 - v_cm,
    )
