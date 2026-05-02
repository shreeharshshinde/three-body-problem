"""
Three-Body Problem Simulation Suite
=====================================
A research-grade Python package for simulating the gravitational
three-body problem with multiple integrators, analysis tools,
and publication-quality visualization.
"""

from .core.equations import (
    acceleration,
    equations_of_motion,
    total_energy,
    total_angular_momentum,
    center_of_mass,
    transform_to_cm_frame,
    G_UNITS,
)
from .core.integrator import (
    integrate_scipy,
    integrate_leapfrog,
    integrate_rk4,
    SimulationResult,
)
from .core.initial_conditions import (
    figure_eight,
    lagrange_triangle,
    pythagorean,
    hierarchical_solar_system,
    figure_eight_perturbed,
    sun_earth_moon,
    custom,
    build_state_vector,
    CATALOGUE,
    list_configurations,
)

__version__ = "1.0.0"
__author__  = "Three-Body Research Project"
__all__ = [
    "acceleration", "equations_of_motion", "total_energy",
    "total_angular_momentum", "center_of_mass", "transform_to_cm_frame",
    "integrate_scipy", "integrate_leapfrog", "integrate_rk4",
    "SimulationResult",
    "figure_eight", "lagrange_triangle", "pythagorean",
    "hierarchical_solar_system", "figure_eight_perturbed",
    "sun_earth_moon", "custom", "build_state_vector",
    "CATALOGUE", "list_configurations", "G_UNITS",
]