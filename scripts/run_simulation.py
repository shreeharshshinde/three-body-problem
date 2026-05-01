"""
run_simulation.py
=================
Command-line entry point for three-body simulations.

Usage
-----
    python scripts/run_simulation.py --config figure_eight
    python scripts/run_simulation.py --config pythagorean --integrator DOP853
    python scripts/run_simulation.py --config lagrange_triangle --animate
    python scripts/run_simulation.py --config solar_system --duration 100 --dt 0.001

Available configurations
------------------------
    figure_eight          Chenciner & Montgomery (2000) choreography
    lagrange_triangle     Lagrange equilateral triangle
    pythagorean           Burrau (1913) chaotic problem
    solar_system          Sun + Jupiter + Saturn
    figure_eight_perturbed  Figure-8 with chaos
    sun_earth_moon        Scaled Solar System
"""

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.core.initial_conditions import CATALOGUE, build_state_vector, list_configurations
from src.core.integrator import integrate_scipy, integrate_leapfrog, integrate_rk4
from src.core.equations import transform_to_cm_frame
from src.visualization.visualize import plot_trajectories, plot_energy, plot_separation, animate_trajectories
from src.analysis.analysis import print_summary, estimate_lyapunov_exponent


def parse_args():
    parser = argparse.ArgumentParser(
        description='Three-Body Problem Simulator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--config', type=str, default='figure_eight',
                        help='Initial condition configuration name')
    parser.add_argument('--integrator', type=str, default='DOP853',
                        choices=['RK45', 'DOP853', 'Leapfrog', 'RK4'],
                        help='Numerical integrator')
    parser.add_argument('--duration', type=float, default=None,
                        help='Simulation duration [yr] (overrides config default)')
    parser.add_argument('--dt', type=float, default=0.001,
                        help='Fixed timestep [yr] (for Leapfrog/RK4)')
    parser.add_argument('--n-points', type=int, default=5000,
                        help='Number of output time points (for scipy integrators)')
    parser.add_argument('--rtol', type=float, default=1e-10,
                        help='Relative tolerance (scipy integrators)')
    parser.add_argument('--animate', action='store_true',
                        help='Save animation to data/outputs/')
    parser.add_argument('--lyapunov', action='store_true',
                        help='Compute Lyapunov exponent')
    parser.add_argument('--list', action='store_true',
                        help='List all available configurations and exit')
    parser.add_argument('--save-plots', action='store_true',
                        help='Save plots to docs/figures/')
    parser.add_argument('--no-plot', action='store_true',
                        help='Skip all plotting')
    return parser.parse_args()


def main():
    args = parse_args()

    if args.list:
        list_configurations()
        return

    if args.config not in CATALOGUE:
        print(f"Error: unknown config '{args.config}'")
        print(f"Available: {list(CATALOGUE.keys())}")
        sys.exit(1)

    # Load initial conditions
    ic = CATALOGUE[args.config]()
    print(f"\n{'='*60}")
    print(f"  {ic['name']}")
    print(f"{'='*60}")
    print(f"  {ic['description']}")
    print(f"  Reference: {ic.get('reference', 'N/A')}")
    print()

    # Transform to center-of-mass frame
    r1, r2, r3, v1, v2, v3 = transform_to_cm_frame(
        ic['r1'], ic['r2'], ic['r3'],
        ic['v1'], ic['v2'], ic['v3'],
        ic['m1'], ic['m2'], ic['m3'],
    )
    ic_cm = dict(ic)
    ic_cm.update({'r1':r1,'r2':r2,'r3':r3,'v1':v1,'v2':v2,'v3':v3})

    y0 = build_state_vector(ic_cm)
    masses = np.array([ic['m1'], ic['m2'], ic['m3']])
    t_end = args.duration or ic.get('t_end', 20.0)
    t_span = (0.0, t_end)

    print(f"  Integrator : {args.integrator}")
    print(f"  Duration   : {t_end} yr")
    print(f"  Masses     : {masses} M_sun")
    print()

    # Run simulation
    print("Integrating...")
    if args.integrator in ('RK45', 'DOP853'):
        t_eval = np.linspace(0, t_end, args.n_points)
        result = integrate_scipy(y0, masses, t_span, t_eval,
                                  method=args.integrator, rtol=args.rtol)
    elif args.integrator == 'Leapfrog':
        result = integrate_leapfrog(y0, masses, t_span, dt=args.dt)
    else:
        result = integrate_rk4(y0, masses, t_span, dt=args.dt)

    print("Done.")
    print_summary(result)

    # Lyapunov exponent
    if args.lyapunov:
        print("Computing Lyapunov exponent (this takes ~30 seconds)...")
        times_l, lyap = estimate_lyapunov_exponent(ic_cm, t_max=min(t_end, 50))
        print(f"  λ_L ≈ {lyap[-1]:.4f} yr⁻¹")
        if lyap[-1] > 0.01:
            print("  → System is CHAOTIC (λ_L > 0)")
        else:
            print("  → System is REGULAR (λ_L ≈ 0)")
        print()

    if args.no_plot:
        return

    # Plotting
    os.makedirs('docs/figures', exist_ok=True)
    os.makedirs('data/outputs', exist_ok=True)

    name = args.config

    fig_traj = plot_trajectories(result, title=ic['name'])
    fig_energy = plot_energy(result)
    fig_sep = plot_separation(result)

    if args.save_plots:
        for tag, fig in [('trajectories', fig_traj), ('energy', fig_energy), ('separation', fig_sep)]:
            path = f'docs/figures/{name}_{tag}.png'
            fig.savefig(path, dpi=150, bbox_inches='tight', facecolor='#0D1117')
            print(f"Saved: {path}")

    if args.animate:
        anim_path = f'data/outputs/{name}_animation.gif'
        print(f"Rendering animation → {anim_path}")
        animate_trajectories(result, output_path=anim_path, fps=20)

    plt.show()


if __name__ == '__main__':
    main()