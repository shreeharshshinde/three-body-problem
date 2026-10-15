#!/usr/bin/env python3
"""
generate_reference_plots.py
============================
Phase 2.1 / 2.4 — Generate reference plots for all 6 configurations
and the Lyapunov convergence plot for the Pythagorean problem.

Outputs (all to docs/figures/):
  trajectory plots:
    figure_eight_trajectory.png
    lagrange_triangle_trajectory.png
    pythagorean_trajectory.png
    solar_system_trajectory.png
    sun_earth_moon_trajectory.png
    figure_eight_perturbed_trajectory.png
  Lyapunov plot:
    pythagorean_lyapunov.png

Usage:
    python scripts/generate_reference_plots.py

Run from the project root.
"""

import os
import sys

# Make sure the package is importable when run from the project root.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib
matplotlib.use('Agg')   # non-interactive backend — safe for scripts
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from src.core.initial_conditions import (
    figure_eight,
    lagrange_triangle,
    pythagorean,
    hierarchical_solar_system,
    sun_earth_moon,
    figure_eight_perturbed,
    build_state_vector,
)
from src.core.integrator import integrate_scipy, integrate_rk4
from src.analysis.analysis import estimate_lyapunov_exponent

# --- Output directory --------------------------------------------------------
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'docs', 'figures')
os.makedirs(OUT_DIR, exist_ok=True)

# --- Dark-theme style (mirrors visualize.py) ---------------------------------
BG      = '#0D1117'
PANEL   = '#161B22'
EDGE    = '#30363D'
FG      = '#C9D1D9'
MUTED   = '#8B949E'
COLORS  = ['#E8632A', '#2A7AE8', '#2AE87A']   # coral, electric blue, mint
LABELS  = ['Body 1', 'Body 2', 'Body 3']

plt.rcParams.update({
    'figure.facecolor': BG,
    'axes.facecolor':   BG,
    'axes.edgecolor':   EDGE,
    'axes.labelcolor':  FG,
    'xtick.color':      MUTED,
    'ytick.color':      MUTED,
    'text.color':       FG,
    'grid.color':       '#21262D',
    'grid.alpha':       0.5,
    'legend.framealpha': 0.25,
    'legend.edgecolor': EDGE,
    'font.family':      'monospace',
    'axes.titlesize':   12,
    'axes.labelsize':   10,
})

DPI = 150


# --- Helpers -----------------------------------------------------------------

def _faded_trail(ax, pos, colors=COLORS):
    """Draw fading trajectory trails (alpha increases toward end)."""
    N = pos.shape[0]
    alphas = np.linspace(0.05, 1.0, max(N - 1, 1))
    for b, color in enumerate(colors):
        x, y = pos[:, b, 0], pos[:, b, 1]
        pts  = np.array([x, y]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        lc   = LineCollection(segs, color=color, alpha=None, linewidths=1.0)
        lc.set_array(alphas)
        lc.set_clim(alphas[0], alphas[-1])
        ax.add_collection(lc)


def _markers(ax, pos, colors=COLORS):
    """Draw start (circle) and end (star) markers for each body."""
    for b, color in enumerate(colors):
        ax.scatter(pos[0, b, 0], pos[0, b, 1],
                   color=color, s=50, zorder=6,
                   marker='o', edgecolors='white', linewidths=0.5)
        ax.scatter(pos[-1, b, 0], pos[-1, b, 1],
                   color=color, s=70, zorder=6, marker='*')


def _info_box(ax, result, extra_lines=None):
    """Add a small text box with simulation metadata."""
    lines = [
        f"t = {result.duration:.1f} yr",
        f"N = {result.n_steps:,}",
        f"Max dE/E0 = {result.max_energy_error:.1e}",
        f"Integrator: {result.integrator}",
    ]
    if extra_lines:
        lines += extra_lines
    ax.text(0.02, 0.98, '\n'.join(lines),
            transform=ax.transAxes,
            fontsize=7.5, va='top', family='monospace', color=MUTED,
            bbox=dict(boxstyle='round,pad=0.4', facecolor=PANEL, alpha=0.8,
                      edgecolor=EDGE))


def save(fig, name):
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=DPI, bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print(f"  OK  {os.path.relpath(path)}")


def _wrap(text, width):
    """Naive word-wrap for a fixed character width."""
    words  = text.split()
    lines  = []
    current = ''
    for w in words:
        if current and len(current) + 1 + len(w) > width:
            lines.append(current)
            current = w
        else:
            current = (current + ' ' + w).strip()
    if current:
        lines.append(current)
    return '\n'.join(lines)


# --- Per-config plot function ------------------------------------------------

def plot_trajectory(ic_fn, filename,
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    t_override=None, n_eval=4000, figsize=(11, 6)):
    """
    Integrate one initial condition and save a trajectory plot.
    """
    ic = ic_fn() if callable(ic_fn) else ic_fn
    masses = np.array([ic['m1'], ic['m2'], ic['m3']])
    y0     = build_state_vector(ic)
    t_end  = t_override if t_override is not None else ic.get('t_end', 20.0)
    t_eval = np.linspace(0, t_end, n_eval)

    result = integrate_scipy(y0, masses, (0, t_end), t_eval,
                             method=method, rtol=rtol, atol=atol)

    pos = result.positions

    # Layout: trajectory (left) + info panel (right)
    fig, (ax, ax_info) = plt.subplots(
        1, 2, figsize=figsize,
        gridspec_kw={'width_ratios': [4, 1]},
    )
    fig.patch.set_facecolor(BG)

    _faded_trail(ax, pos)
    _markers(ax, pos)

    ax.autoscale()
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_xlabel('x  [AU]')
    ax.set_ylabel('y  [AU]')
    ax.set_title(ic.get('name', filename.replace('_trajectory.png', '')),
                 fontsize=13, pad=10)
    ax.grid(True, alpha=0.3)

    # Proxy artists for legend
    for b, (color, label) in enumerate(zip(COLORS, LABELS)):
        ax.plot([], [], color=color, lw=2, label=label)
    ax.legend(loc='upper right', fontsize=8)

    _info_box(ax, result)

    # Right panel: description
    ax_info.set_facecolor(BG)
    ax_info.axis('off')
    desc = ic.get('description', '')
    ref  = ic.get('reference', '')
    wrap_desc = _wrap(desc, 28)
    wrap_ref  = _wrap(ref, 28)
    text = f"{wrap_desc}\n\n-- {wrap_ref}"
    ax_info.text(0.05, 0.98, text,
                 transform=ax_info.transAxes,
                 fontsize=7.5, va='top', color=MUTED,
                 linespacing=1.6)

    fig.tight_layout(pad=1.5)
    save(fig, filename)


# --- Lyapunov convergence ----------------------------------------------------

def plot_lyapunov_pythagorean(filename='pythagorean_lyapunov.png', figsize=(10, 5)):
    """
    Compute and plot Lyapunov exponent convergence for the Pythagorean problem.
    Uses the Benettin et al. (1980) renormalization algorithm.
    """
    print("  Computing Lyapunov exponent for Pythagorean orbit (this may take ~30 s)...")
    ic = pythagorean()
    times, lce = estimate_lyapunov_exponent(
        ic, t_max=60.0, dt=0.005, epsilon=1e-6, n_renorm=200
    )

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor(BG)

    ax.plot(times, lce, color=COLORS[0], lw=1.5, label='lambda_L(t)  [running mean]')
    ax.axhline(lce[-1], color=COLORS[2], ls='--', lw=1.2,
               label=f'Converged:  lambda_L = {lce[-1]:.4f} yr-1')
    ax.axhline(0, color=MUTED, ls=':', lw=0.8)

    ax.set_xlabel('Time [yr]')
    ax.set_ylabel('Lyapunov Exponent  lambda_L  [yr-1]')
    ax.set_title('Maximal Lyapunov Exponent -- Pythagorean Problem (Burrau 1913)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.35)

    ax.text(0.98, 0.08,
            'lambda_L > 0  ->  chaotic system\n'
            'Same quantity controls\n'
            'photon ring widths in Kerr BH',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=8.5, color=MUTED,
            bbox=dict(boxstyle='round,pad=0.45', facecolor=PANEL,
                      alpha=0.85, edgecolor=EDGE))

    fig.tight_layout(pad=1.5)
    save(fig, filename)


# --- Main --------------------------------------------------------------------

def main():
    print("\n--- Phase 2.1 / 2.4 -- Generating reference plots ------------------\n")

    configs = [
        # (ic_function,             output_filename,                  kwargs)
        (figure_eight,             'figure_eight_trajectory.png',
         dict(t_override=20.0,   n_eval=6000)),

        (lagrange_triangle,        'lagrange_triangle_trajectory.png',
         dict(n_eval=4000)),

        (pythagorean,              'pythagorean_trajectory.png',
         dict(t_override=70.0,   n_eval=8000, rtol=1e-11, atol=1e-13)),

        (hierarchical_solar_system,'solar_system_trajectory.png',
         dict(t_override=100.0,  n_eval=6000)),

        (sun_earth_moon,           'sun_earth_moon_trajectory.png',
         dict(t_override=5.0,    n_eval=6000)),

        (figure_eight_perturbed,   'figure_eight_perturbed_trajectory.png',
         dict(t_override=30.0,   n_eval=6000)),
    ]

    for ic_fn, fname, kw in configs:
        print(f"  Plotting: {fname}")
        plot_trajectory(ic_fn, fname, **kw)

    print("\n  Plotting: pythagorean_lyapunov.png")
    plot_lyapunov_pythagorean()

    print(f"\n--- Done. {len(configs) + 1} figures written to docs/figures/ -----\n")


if __name__ == '__main__':
    main()
