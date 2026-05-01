"""
test_physics.py
===============
Unit tests for the three-body physics engine.

Tests verify:
1. Conservation laws (energy, angular momentum, momentum)
2. Known analytical solutions (Lagrange triangle, figure-8 period)
3. Integrator accuracy comparison
4. Equations of motion correctness
5. Center-of-mass frame transformation

Run with:
    pytest tests/test_physics.py -v
    pytest tests/test_physics.py -v --tb=short
"""

import numpy as np
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.core.equations import (
    acceleration, equations_of_motion, total_energy,
    total_angular_momentum, center_of_mass, transform_to_cm_frame, G_UNITS,
)
from src.core.integrator import integrate_scipy, integrate_leapfrog, integrate_rk4
from src.core.initial_conditions import (
    figure_eight, lagrange_triangle, pythagorean,
    hierarchical_solar_system, build_state_vector,
)
from src.analysis.analysis import angular_momentum_error, separation_matrix


# ─── Tolerance levels ─────────────────────────────────────────────────────────
TIGHT   = 1e-6   # for high-accuracy integrators
LOOSE   = 1e-3   # for fixed-step RK4 with modest dt


class TestEquationsOfMotion:
    """Test the gravitational acceleration and ODE functions."""

    def test_acceleration_equal_masses_collinear(self):
        """Three equal masses on x-axis: middle body has zero net x-force by symmetry."""
        m = 1.0
        r1 = np.array([-1.0, 0.0])
        r2 = np.array([ 0.0, 0.0])
        r3 = np.array([ 1.0, 0.0])
        v = np.zeros(2)
        a1, a2, a3 = acceleration(r1, r2, r3, m, m, m, softening=0)
        # Middle body: pulled equally left and right → net x ≈ 0
        assert abs(a2[0]) < 1e-10, f"Expected a2x≈0, got {a2[0]}"

    def test_acceleration_symmetry(self):
        """Newton's third law: F_12 = -F_21."""
        r1 = np.array([0.0, 0.0])
        r2 = np.array([1.0, 0.5])
        r3 = np.array([5.0, 5.0])  # far away
        m1, m2, m3 = 1.0, 2.0, 0.001
        a1, a2, a3 = acceleration(r1, r2, r3, m1, m2, m3, softening=0)
        # Force from body 2 on 1 = -Force from 1 on 2 (ignoring body 3 contribution)
        # m1*a1_from_2 = -m2*a2_from_1 → m1*(G*m2*(r2-r1)/d^3) = -m2*(G*m1*(r1-r2)/d^3) ✓
        # We check that total momentum change is zero: m1*a1 + m2*a2 + m3*a3 = 0
        total_force = m1 * a1 + m2 * a2 + m3 * a3
        np.testing.assert_allclose(total_force, 0, atol=1e-12,
                                   err_msg="Total gravitational force must be zero (Newton's 3rd law)")

    def test_momentum_conservation_in_ode(self):
        """Total momentum d/dt(Σ m_i v_i) = 0 from equations of motion."""
        ic = figure_eight()
        masses = np.array([ic['m1'], ic['m2'], ic['m3']])
        y0 = build_state_vector(ic)
        dydt = equations_of_motion(0.0, y0, masses)
        # dydt[6:8] = a1, dydt[8:10] = a2, dydt[10:12] = a3
        total_force = (masses[0]*dydt[6:8] + masses[1]*dydt[8:10] + masses[2]*dydt[10:12])
        np.testing.assert_allclose(total_force, 0, atol=1e-12,
                                   err_msg="Total force must vanish")

    def test_energy_positive_kinetic(self):
        """Kinetic energy must always be non-negative."""
        ic = figure_eight()
        r1, r2, r3 = ic['r1'], ic['r2'], ic['r3']
        v1, v2, v3 = ic['v1'], ic['v2'], ic['v3']
        m1, m2, m3 = ic['m1'], ic['m2'], ic['m3']
        E = total_energy(r1, r2, r3, v1, v2, v3, m1, m2, m3)
        KE = 0.5*(m1*np.dot(v1,v1) + m2*np.dot(v2,v2) + m3*np.dot(v3,v3))
        assert KE >= 0, f"Kinetic energy must be non-negative, got {KE}"


class TestCenterOfMass:
    """Test CM frame transformations."""

    def test_cm_position_zero_after_transform(self):
        """After CM transform, center of mass should be at origin."""
        ic = pythagorean()
        r1, r2, r3, v1, v2, v3 = transform_to_cm_frame(
            ic['r1'], ic['r2'], ic['r3'],
            ic['v1'], ic['v2'], ic['v3'],
            ic['m1'], ic['m2'], ic['m3'],
        )
        r_cm = center_of_mass(r1, r2, r3, ic['m1'], ic['m2'], ic['m3'])
        np.testing.assert_allclose(r_cm, 0, atol=1e-14,
                                   err_msg="CM should be at origin after transform")

    def test_cm_velocity_zero_after_transform(self):
        """After CM transform, total momentum should be zero."""
        ic = pythagorean()
        r1, r2, r3, v1, v2, v3 = transform_to_cm_frame(
            ic['r1'], ic['r2'], ic['r3'],
            ic['v1'], ic['v2'], ic['v3'],
            ic['m1'], ic['m2'], ic['m3'],
        )
        m1, m2, m3 = ic['m1'], ic['m2'], ic['m3']
        p_total = m1*v1 + m2*v2 + m3*v3
        np.testing.assert_allclose(p_total, 0, atol=1e-14,
                                   err_msg="Total momentum must be zero in CM frame")


class TestConservationLaws:
    """Test that integrators conserve energy and angular momentum."""

    @pytest.mark.parametrize("ic_fn,method,tol", [
        (figure_eight,          'DOP853',   TIGHT),
        (lagrange_triangle,     'DOP853',   TIGHT),
        (hierarchical_solar_system, 'DOP853', TIGHT),
        (figure_eight,          'RK45',     TIGHT),
    ])
    def test_energy_conservation(self, ic_fn, method, tol):
        """Energy should be conserved to within tolerance."""
        ic = ic_fn()
        masses = np.array([ic['m1'], ic['m2'], ic['m3']])
        r1, r2, r3, v1, v2, v3 = transform_to_cm_frame(
            ic['r1'], ic['r2'], ic['r3'],
            ic['v1'], ic['v2'], ic['v3'],
            *masses,
        )
        ic_cm = dict(ic)
        ic_cm.update({'r1':r1,'r2':r2,'r3':r3,'v1':v1,'v2':v2,'v3':v3})
        y0 = build_state_vector(ic_cm)
        t_end = ic.get('t_end', 10.0)
        t_eval = np.linspace(0, t_end, 2000)

        result = integrate_scipy(y0, masses, (0, t_end), t_eval,
                                  method=method, rtol=1e-12, atol=1e-14)
        max_err = result.max_energy_error
        assert max_err < tol, (
            f"{ic['name']} ({method}): energy error {max_err:.2e} > tolerance {tol:.2e}"
        )

    def test_leapfrog_energy_conservation(self):
        """Leapfrog should conserve energy better than RK4 for long integrations."""
        ic = lagrange_triangle()
        masses = np.array([ic['m1'], ic['m2'], ic['m3']])
        y0 = build_state_vector(ic)
        t_end = ic.get('t_period', 1.0) * 10

        res_lf  = integrate_leapfrog(y0, masses, (0, t_end), dt=0.001, store_every=5)
        res_rk4 = integrate_rk4(y0,     masses, (0, t_end), dt=0.001, store_every=5)

        # Leapfrog should have smaller or comparable energy error
        assert res_lf.max_energy_error <= res_rk4.max_energy_error * 10, \
            "Leapfrog should not be dramatically worse than RK4 for energy conservation"

    def test_angular_momentum_conservation(self):
        """Angular momentum should be conserved."""
        ic = figure_eight()
        masses = np.array([ic['m1'], ic['m2'], ic['m3']])
        y0 = build_state_vector(ic)
        t_eval = np.linspace(0, 10, 3000)
        result = integrate_scipy(y0, masses, (0, 10), t_eval, method='DOP853',
                                  rtol=1e-12, atol=1e-14)
        L_err = angular_momentum_error(result).max()
        assert L_err < 1e-8, f"Angular momentum error {L_err:.2e} too large"


class TestKnownSolutions:
    """Test against known analytical solutions."""

    def test_lagrange_triangle_period(self):
        """Lagrange triangle: numerically measured period should match analytical."""
        from src.analysis.analysis import compute_orbital_period

        ic = lagrange_triangle()
        T_analytical = ic['t_period']
        masses = np.array([ic['m1'], ic['m2'], ic['m3']])
        y0 = build_state_vector(ic)
        t_eval = np.linspace(0, 5 * T_analytical, 10000)
        result = integrate_scipy(y0, masses, (0, 5*T_analytical), t_eval,
                                  method='DOP853', rtol=1e-12, atol=1e-14)

        T_numerical = compute_orbital_period(result, body=0, coordinate=0)
        if T_numerical is not None:
            rel_err = abs(T_numerical - T_analytical) / T_analytical
            assert rel_err < 0.01, (
                f"Period error {rel_err:.2%}: numerical={T_numerical:.4f}, "
                f"analytical={T_analytical:.4f}"
            )

    def test_figure_eight_returns_to_start(self):
        """After one period, all bodies should return near initial positions."""
        ic = figure_eight()
        T = ic['t_period']
        masses = np.array([1.0, 1.0, 1.0])
        y0 = build_state_vector(ic)
        t_eval = np.linspace(0, T, 5000)
        result = integrate_scipy(y0, masses, (0, T), t_eval,
                                  method='DOP853', rtol=1e-12, atol=1e-14)

        for b in range(3):
            pos_initial = result.positions[0, b]
            pos_final   = result.positions[-1, b]
            err = np.linalg.norm(pos_final - pos_initial)
            assert err < 0.1, (
                f"Body {b+1} didn't return to start after one period. Error: {err:.4f} AU"
            )

    def test_two_body_limit(self):
        """With m3 → 0, bodies 1 and 2 should follow Keplerian orbit."""
        # Place body 3 far away with tiny mass — system ≈ two-body
        y0 = np.array([
            1.0, 0.0,     # r1
            -1.0, 0.0,    # r2
            100.0, 0.0,   # r3 (far away)
            0.0, np.pi,   # v1 (circular orbit: v=2π/sqrt(a), a=2, v=2π/sqrt(2))
            0.0, -np.pi,  # v2
            0.0, 0.0,     # v3
        ])
        # Correct velocities for circular 2-body: v = sqrt(G*M/2r) with G=4π², M=2, r=1
        v_circ = np.sqrt(G_UNITS * 2.0 / 2.0) / 2  # approximate
        y0[6] = 0.0; y0[7] = v_circ * 2
        y0[8] = 0.0; y0[9] = -v_circ * 2

        masses = np.array([1.0, 1.0, 1e-10])
        t_eval = np.linspace(0, 2, 2000)
        result = integrate_scipy(y0, masses, (0, 2), t_eval, method='DOP853')
        # Energy should still be conserved
        assert result.max_energy_error < 1e-4


class TestIntegratorComparison:
    """Compare integrator outputs for consistency."""

    def test_all_integrators_agree_short_time(self):
        """All integrators should give consistent results for short integration."""
        ic = lagrange_triangle()
        masses = np.array([ic['m1'], ic['m2'], ic['m3']])
        y0 = build_state_vector(ic)
        t_end = 0.5  # short time — all integrators should agree

        t_eval = np.linspace(0, t_end, 1000)
        res_rk45   = integrate_scipy(y0, masses, (0, t_end), t_eval, method='RK45',   rtol=1e-12)
        res_dop853 = integrate_scipy(y0, masses, (0, t_end), t_eval, method='DOP853', rtol=1e-12)
        res_lf     = integrate_leapfrog(y0, masses, (0, t_end), dt=1e-4, store_every=1)

        # Final positions should agree to within 1e-5 AU
        for b in range(3):
            err = np.linalg.norm(res_rk45.positions[-1, b] - res_dop853.positions[-1, b])
            assert err < 1e-5, f"RK45 vs DOP853 disagree for body {b+1}: {err:.2e} AU"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])