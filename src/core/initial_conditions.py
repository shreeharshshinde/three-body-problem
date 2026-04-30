"""
initial_conditions.py
=====================
Famous and physically meaningful initial conditions for the three-body problem.

Each function returns a dictionary with:
    'r1','r2','r3'  : initial positions [AU]
    'v1','v2','v3'  : initial velocities [AU/yr]
    'm1','m2','m3'  : masses [M_sun]
    'name'          : human-readable name
    'description'   : physics description
    't_period'      : estimated orbital period [yr] (if periodic)
    'reference'     : original paper

Catalogue
---------
1. figure_eight       — Chenciner & Montgomery (2000), equal masses
2. lagrange_triangle  — Lagrange (1772), equilateral triangle
3. euler_collinear    — Euler (1767), collinear configuration
4. broucke_henon      — Broucke (1975), stable family
5. hierarchical       — Sun + Jupiter + Saturn (real solar system)
6. pythagorean        — Burrau (1913), masses 3-4-5, chaotic
7. figure_eight_perturbed — Figure-8 with perturbation (shows chaos)
8. sun_earth_moon     — Scaled Sun-Earth-Moon system
9. random_stable      — Random configuration with zero total momentum

References
----------
- Chenciner & Montgomery (2000), Ann. Math. 152, 881
- Moore (1993), Phys. Rev. Lett. 70, 3675
- Šuvakov & Dmitrašinović (2013), Phys. Rev. Lett. 110, 114301
- Lagrange (1772), Prix de l'Académie Royale
- Burrau (1913), Vierteljahrsschr. Astron. Ges. 48, 195
"""

import numpy as np
from typing import Dict, Any


def figure_eight() -> Dict[str, Any]:
    """
    The figure-8 choreography — Chenciner & Montgomery (2000).

    Three equal-mass bodies chase each other along a figure-8 curve.
    Discovered numerically by Moore (1993), proven to exist by
    Chenciner & Montgomery (2000) using variational methods.

    This is one of the most famous solutions to the three-body problem.
    All three bodies have equal mass m = 1 M_sun and follow the
    same figure-8 orbit, separated by 1/3 of the period.

    Period: T ≈ 6.3259 (in normalized units; ≈ 1.007 yr here)
    """
    # Positions from Chenciner & Montgomery (2000)
    x1 =  0.97000436; y1 = -0.24308753
    x2 = -0.97000436; y2 =  0.24308753
    x3 =  0.0;        y3 =  0.0

    # Velocities (specific to the choreography constraint)
    vx3 = -0.93240737 * 2; vy3 = -0.86473146 * 2
    vx1 = -vx3 / 2;        vy1 = -vy3 / 2
    vx2 = -vx3 / 2;        vy2 = -vy3 / 2

    return {
        'r1': np.array([x1, y1]),
        'r2': np.array([x2, y2]),
        'r3': np.array([x3, y3]),
        'v1': np.array([vx1, vy1]),
        'v2': np.array([vx2, vy2]),
        'v3': np.array([vx3, vy3]),
        'm1': 1.0, 'm2': 1.0, 'm3': 1.0,
        'name': 'Figure-Eight Choreography',
        'description': (
            'Three equal-mass bodies chasing each other on a figure-8 curve. '
            'Discovered by Moore (1993), proven by Chenciner & Montgomery (2000). '
            'A periodic choreography solution — all bodies follow the same path.'
        ),
        't_period': 6.3259 / (2 * np.pi),
        't_end': 20.0,
        'reference': 'Chenciner & Montgomery (2000), Ann. Math. 152, 881',
    }


def lagrange_triangle() -> Dict[str, Any]:
    """
    Lagrange's equilateral triangle solution (1772).

    Three bodies at the vertices of an equilateral triangle,
    orbiting the common center of mass. Stable for m1 >> m2 + m3.
    This is the origin of the Lagrange points L4 and L5 in celestial
    mechanics — the Trojan asteroids occupy these points.

    The triangle rotates rigidly. Angular velocity:
        ω² = G(m1+m2+m3) / L³
    where L is the side length.
    """
    L = 2.0  # side length [AU]
    m1, m2, m3 = 1.0, 1.0, 1.0
    M = m1 + m2 + m3

    # Equilateral triangle vertices
    h = L * np.sqrt(3) / 2
    r1 = np.array([ L/2,    0.0])
    r2 = np.array([-L/2,    0.0])
    r3 = np.array([ 0.0,    h  ])

    # Shift to center of mass
    r_cm = (m1*r1 + m2*r2 + m3*r3) / M
    r1 -= r_cm; r2 -= r_cm; r3 -= r_cm

    # Angular velocity for rigid rotation
    G = 4 * np.pi**2
    omega = np.sqrt(G * M / L**3)

    # Velocities: v = ω × r (perpendicular to r, in 2D: v = ω * (-y, x))
    v1 = omega * np.array([-r1[1],  r1[0]])
    v2 = omega * np.array([-r2[1],  r2[0]])
    v3 = omega * np.array([-r3[1],  r3[0]])

    T = 2 * np.pi / omega

    return {
        'r1': r1, 'r2': r2, 'r3': r3,
        'v1': v1, 'v2': v2, 'v3': v3,
        'm1': m1, 'm2': m2, 'm3': m3,
        'name': 'Lagrange Equilateral Triangle',
        'description': (
            'Lagrange (1772): three equal-mass bodies at vertices of an '
            'equilateral triangle, rotating rigidly. The triangle rotates '
            'at angular velocity ω = √(GM/L³). Source of L4, L5 Lagrange points.'
        ),
        't_period': T,
        't_end': 3 * T,
        'reference': 'Lagrange (1772), Prix de l\'Académie Royale des Sciences',
    }


def pythagorean() -> Dict[str, Any]:
    """
    Burrau's Pythagorean Problem (1913) — the classic chaotic case.

    Three bodies with masses 3, 4, 5 (Pythagorean triple) placed at
    rest at the vertices of a right triangle with sides 3, 4, 5.
    First numerically integrated by Burrau (1913); proven to end in
    ejection of the lightest body by Szebehely & Peters (1967).

    This is the canonical demonstration that the three-body problem
    is generically chaotic. The lightest body (m=3) gets ejected.
    """
    return {
        'r1': np.array([1.0,  3.0]),
        'r2': np.array([-2.0, -1.0]),
        'r3': np.array([1.0,  -1.0]),
        'v1': np.zeros(2),
        'v2': np.zeros(2),
        'v3': np.zeros(2),
        'm1': 3.0, 'm2': 4.0, 'm3': 5.0,
        'name': "Pythagorean Problem (Burrau 1913)",
        'description': (
            'Masses 3,4,5 at rest at vertices of a right triangle. '
            'The first famous numerical three-body calculation (Burrau 1913). '
            'Chaotic: the lightest body is eventually ejected. '
            'Demonstrates generic instability of equal-energy configurations.'
        ),
        't_period': None,
        't_end': 70.0,
        'reference': 'Burrau (1913); Szebehely & Peters (1967), AJ 72, 876',
    }


def hierarchical_solar_system() -> Dict[str, Any]:
    """
    Hierarchical system: Sun + Jupiter + Saturn (scaled).

    A stable hierarchical triple: a massive central body with two
    smaller bodies in well-separated orbits. Jupiter and Saturn are
    not in a simple periodic orbit but their 5:2 near-resonance
    (the 'Great Inequality') causes long-period perturbations.

    Units scaled: 1 AU, 1 yr, G = 4π²
    """
    # Sun
    r1 = np.array([0.0, 0.0])
    v1 = np.array([0.0, 0.0])
    m1 = 1.0  # M_sun

    # Jupiter: a = 5.2 AU, circular orbit
    a_J = 5.2
    v_J = 2 * np.pi * np.sqrt(1.0 / a_J)   # Kepler: v = 2π√(GM/a), G=4π², M=1 → v=2π/√a
    # Correct formula with G=4π²: v = √(G*M/a) = √(4π²/a) = 2π/√a
    v_J = 2 * np.pi / np.sqrt(a_J)
    r2 = np.array([a_J, 0.0])
    v2 = np.array([0.0, v_J])
    m2 = 9.548e-4  # Jupiter in M_sun

    # Saturn: a = 9.58 AU, circular orbit
    a_S = 9.58
    v_S = 2 * np.pi / np.sqrt(a_S)
    r3 = np.array([-a_S, 0.0])
    v3 = np.array([0.0, -v_S])
    m3 = 2.858e-4  # Saturn in M_sun

    return {
        'r1': r1, 'r2': r2, 'r3': r3,
        'v1': v1, 'v2': v2, 'v3': v3,
        'm1': m1, 'm2': m2, 'm3': m3,
        'name': 'Sun-Jupiter-Saturn',
        'description': (
            'Hierarchical triple: Sun (1 M☉), Jupiter (a=5.2 AU), Saturn (a=9.58 AU). '
            'Stable over long timescales. Shows 5:2 orbital near-resonance. '
            'Jupiter period ≈ 11.9 yr, Saturn ≈ 29.5 yr.'
        ),
        't_period': 11.86,  # Jupiter period [yr]
        't_end': 100.0,
        'reference': 'NASA/JPL Solar System Parameters',
    }


def figure_eight_perturbed(perturbation: float = 0.01) -> Dict[str, Any]:
    """
    Figure-8 with a small perturbation — demonstrates chaos.

    A tiny displacement from the exact figure-8 initial conditions
    grows exponentially (Lyapunov instability). The orbit initially
    follows the figure-8 but eventually diverges.

    This directly illustrates the concept of the Lyapunov exponent —
    the same mathematical object you study in your black hole research!
    (The Lyapunov exponent of photon orbits near the photon sphere
    controls photon ring width.)

    Parameters
    ----------
    perturbation : float
        Size of position perturbation [AU]. Default 0.01 AU.
    """
    ic = figure_eight()
    rng = np.random.default_rng(42)
    ic['r1'] += perturbation * rng.standard_normal(2)
    ic['name'] = f'Figure-8 Perturbed (δ={perturbation} AU)'
    ic['description'] = (
        f'Figure-8 choreography with {perturbation} AU perturbation. '
        'The orbit is initially figure-8-like but diverges exponentially '
        'due to Lyapunov instability — the same mathematical mechanism '
        'that controls photon ring widths in black hole spacetimes.'
    )
    ic['t_end'] = 30.0
    return ic


def sun_earth_moon() -> Dict[str, Any]:
    """
    Scaled Sun-Earth-Moon system.

    The Moon's orbit is perturbed by the Sun, causing precession of
    the lunar orbit with period ~18.6 years (nodal precession) and
    8.85-year apsidal precession. This is a real three-body problem
    with no closed-form solution.

    Scaled for visualization: distances in AU, times in years.
    """
    # Earth: 1 AU, circular orbit
    a_E = 1.0
    v_E = 2 * np.pi / np.sqrt(a_E)

    # Moon: 384,400 km ≈ 0.00257 AU from Earth
    a_M_relative = 0.00257  # AU
    v_M_relative = 2 * np.pi / np.sqrt(a_M_relative) * np.sqrt(3e-6)  # Earth mass = 3e-6 M_sun

    return {
        'r1': np.array([0.0, 0.0]),
        'r2': np.array([a_E, 0.0]),
        'r3': np.array([a_E + a_M_relative, 0.0]),
        'v1': np.array([0.0, 0.0]),
        'v2': np.array([0.0, v_E]),
        'v3': np.array([0.0, v_E + v_M_relative]),
        'm1': 1.0,       # Sun [M_sun]
        'm2': 3.003e-6,  # Earth [M_sun]
        'm3': 3.69e-8,   # Moon [M_sun]
        'name': 'Sun-Earth-Moon',
        'description': (
            'Real Sun-Earth-Moon system (scaled). The Moon\'s orbit precesses '
            'due to solar perturbations. No closed-form solution exists. '
            'Historically motivated the development of perturbation theory.'
        ),
        't_period': 1.0,  # Earth year
        't_end': 5.0,
        'reference': 'Newton (1687), Principia Mathematica, Book III',
    }


def custom(
    r1, r2, r3, v1, v2, v3,
    m1=1.0, m2=1.0, m3=1.0,
    name="Custom", description="User-defined initial conditions",
    t_end=20.0,
) -> Dict[str, Any]:
    """
    Create a custom initial condition dictionary.

    Parameters
    ----------
    r1, r2, r3 : array-like, shape (2,)
        Initial positions [AU]
    v1, v2, v3 : array-like, shape (2,)
        Initial velocities [AU/yr]
    m1, m2, m3 : float
        Masses [M_sun]
    """
    return {
        'r1': np.array(r1), 'r2': np.array(r2), 'r3': np.array(r3),
        'v1': np.array(v1), 'v2': np.array(v2), 'v3': np.array(v3),
        'm1': float(m1), 'm2': float(m2), 'm3': float(m3),
        'name': name, 'description': description,
        't_period': None, 't_end': float(t_end),
        'reference': 'User-defined',
    }


def build_state_vector(ic: Dict[str, Any]) -> np.ndarray:
    """
    Flatten an initial condition dict into the ODE state vector.

    Returns
    -------
    y0 : np.ndarray, shape (12,)
        [x1,y1, x2,y2, x3,y3, vx1,vy1, vx2,vy2, vx3,vy3]
    """
    return np.concatenate([
        ic['r1'], ic['r2'], ic['r3'],
        ic['v1'], ic['v2'], ic['v3'],
    ])


# Registry of all built-in configurations
CATALOGUE = {
    'figure_eight':         figure_eight,
    'lagrange_triangle':    lagrange_triangle,
    'pythagorean':          pythagorean,
    'solar_system':         hierarchical_solar_system,
    'figure_eight_perturbed': figure_eight_perturbed,
    'sun_earth_moon':       sun_earth_moon,
}


def list_configurations() -> None:
    """Print a formatted table of all available configurations."""
    print(f"\n{'Name':<30} {'Description':<60}")
    print('─' * 90)
    for key, fn in CATALOGUE.items():
        ic = fn()
        desc = ic['description'][:57] + '...' if len(ic['description']) > 57 else ic['description']
        print(f"{key:<30} {desc:<60}")
    print()