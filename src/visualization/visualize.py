"""
visualize.py
============
Visualization tools for three-body simulation results.

Functions
---------
plot_trajectories       : Static trajectory plot with energy conservation
plot_energy             : Energy and angular momentum conservation over time
plot_phase_space        : Phase space portrait (position vs velocity)
plot_separation         : Pairwise body separations over time
animate_trajectories    : Animated orbit movie (saves MP4 or GIF)
plot_lyapunov           : Lyapunov exponent convergence plot
plot_poincare_section   : Poincaré surface of section (chaos indicator)

All functions return matplotlib Figure objects for further customization.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
from typing import Optional, Tuple
from ..core.integrator import SimulationResult


# ─── Color palette ────────────────────────────────────────────────────────────
BODY_COLORS = ['#E8632A', '#2A7AE8', '#2AE87A']   # coral, blue, green
BODY_LABELS = ['Body 1', 'Body 2', 'Body 3']

plt.rcParams.update({
    'figure.facecolor': '#0D1117',
    'axes.facecolor':   '#0D1117',
    'axes.edgecolor':   '#30363D',
    'axes.labelcolor':  '#C9D1D9',
    'xtick.color':      '#8B949E',
    'ytick.color':      '#8B949E',
    'text.color':       '#C9D1D9',
    'grid.color':       '#21262D',
    'grid.alpha':       0.6,
    'legend.framealpha': 0.3,
    'legend.edgecolor': '#30363D',
    'font.family':      'monospace',
})


def plot_trajectories(
    result: SimulationResult,
    title: Optional[str] = None,
    fade: bool = True,
    show_markers: bool = True,
    figsize: Tuple[int, int] = (10, 8),
) -> plt.Figure:
    """
    Plot three-body trajectories with optional fading trail.

    Parameters
    ----------
    result : SimulationResult
    title : str, optional
        Plot title. Defaults to integrator name.
    fade : bool
        Use alpha-faded trail to indicate time direction.
    show_markers : bool
        Show start/end position markers.
    figsize : tuple

    Returns
    -------
    fig : matplotlib.Figure
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize,
                              gridspec_kw={'width_ratios': [3, 1]})
    ax = axes[0]
    ax_info = axes[1]

    pos = result.positions
    N = result.n_steps

    if fade:
        # Fading trail: alpha increases toward end of trajectory
        alphas = np.linspace(0.05, 0.9, N - 1)
        for b, color in enumerate(BODY_COLORS):
            x = pos[:, b, 0]
            y = pos[:, b, 1]
            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, color=color, alpha=None)
            lc.set_array(alphas)
            lc.set_clim(0.05, 0.9)
            ax.add_collection(lc)
    else:
        for b, (color, label) in enumerate(zip(BODY_COLORS, BODY_LABELS)):
            ax.plot(pos[:, b, 0], pos[:, b, 1], color=color,
                    lw=0.8, alpha=0.7, label=label)

    if show_markers:
        for b, color in enumerate(BODY_COLORS):
            # Start: circle
            ax.scatter(pos[0, b, 0], pos[0, b, 1], color=color,
                       s=60, zorder=5, marker='o', edgecolors='white', linewidths=0.5)
            # End: star
            ax.scatter(pos[-1, b, 0], pos[-1, b, 1], color=color,
                       s=80, zorder=5, marker='*')

    ax.set_xlabel('x  [AU]')
    ax.set_ylabel('y  [AU]')
    ax.set_title(title or f'Three-Body Trajectories ({result.integrator})')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(BODY_LABELS, loc='upper right', fontsize=9)

    # Info panel
    ax_info.axis('off')
    masses = result.masses
    info = (
        f"Integrator: {result.integrator}\n"
        f"Duration: {result.duration:.2f} yr\n"
        f"Steps: {result.n_steps}\n\n"
        f"Masses [M☉]:\n"
        f"  m₁ = {masses[0]:.4f}\n"
        f"  m₂ = {masses[1]:.4f}\n"
        f"  m₃ = {masses[2]:.4f}\n\n"
        f"Max ΔE/E₀:\n"
        f"  {result.max_energy_error:.2e}\n"
    )
    ax_info.text(0.05, 0.95, info, transform=ax_info.transAxes,
                 fontsize=9, verticalalignment='top', family='monospace',
                 color='#C9D1D9', linespacing=1.8)

    fig.tight_layout()
    return fig


def plot_energy(result: SimulationResult, figsize=(12, 5)) -> plt.Figure:
    """
    Plot energy and angular momentum conservation over time.

    A well-integrated simulation should show flat lines.
    Any slope or oscillation indicates numerical error.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    from ..analysis.analysis import angular_momentum_error

    # Energy conservation
    ax1.semilogy(result.t, result.energy_error + 1e-20, color='#E8632A', lw=1.2)
    ax1.set_xlabel('Time [yr]')
    ax1.set_ylabel('|ΔE / E₀|')
    ax1.set_title('Energy Conservation Error')
    ax1.axhline(1e-6, color='#2AE87A', ls='--', lw=0.8, label='1e-6 threshold')
    ax1.legend(fontsize=8)
    ax1.grid(True)

    # Angular momentum conservation
    L_err = angular_momentum_error(result)
    ax2.semilogy(result.t, L_err + 1e-20, color='#2A7AE8', lw=1.2)
    ax2.set_xlabel('Time [yr]')
    ax2.set_ylabel('|ΔL / L₀|')
    ax2.set_title('Angular Momentum Conservation Error')
    ax2.grid(True)

    fig.suptitle(f'Conservation Laws — {result.integrator} integrator')
    fig.tight_layout()
    return fig


def plot_separation(result: SimulationResult, figsize=(10, 5)) -> plt.Figure:
    """
    Plot pairwise separation between bodies over time.

    Spikes in separation correspond to close encounters;
    diverging separations indicate escape/ejection.
    """
    from ..analysis.analysis import separation_matrix

    seps = separation_matrix(result)
    fig, ax = plt.subplots(figsize=figsize)

    pair_labels = ['|r₁ - r₂|', '|r₁ - r₃|', '|r₂ - r₃|']
    pair_colors = ['#E8632A', '#2A7AE8', '#2AE87A']

    for i, (label, color) in enumerate(zip(pair_labels, pair_colors)):
        ax.plot(result.t, seps[:, i], label=label, color=color, lw=1.0)

    ax.set_xlabel('Time [yr]')
    ax.set_ylabel('Separation [AU]')
    ax.set_title('Pairwise Body Separations')
    ax.legend(fontsize=9)
    ax.grid(True)
    fig.tight_layout()
    return fig


def plot_phase_space(result: SimulationResult, body: int = 0, figsize=(8, 6)) -> plt.Figure:
    """
    Phase space portrait: position vs. velocity for one body.

    A periodic orbit → closed curve.
    A chaotic orbit → space-filling curve (strange attractor).
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    pos = result.positions[:, body, :]
    vel = result.velocities[:, body, :]
    color = BODY_COLORS[body]

    axes[0].plot(pos[:, 0], vel[:, 0], color=color, lw=0.5, alpha=0.7)
    axes[0].set_xlabel('x [AU]')
    axes[0].set_ylabel('vx [AU/yr]')
    axes[0].set_title(f'Phase space (x, vx) — Body {body+1}')
    axes[0].grid(True)

    axes[1].plot(pos[:, 1], vel[:, 1], color=color, lw=0.5, alpha=0.7)
    axes[1].set_xlabel('y [AU]')
    axes[1].set_ylabel('vy [AU/yr]')
    axes[1].set_title(f'Phase space (y, vy) — Body {body+1}')
    axes[1].grid(True)

    fig.suptitle(f'Phase Space Portrait — {BODY_LABELS[body]}')
    fig.tight_layout()
    return fig


def plot_lyapunov(times: np.ndarray, lyapunov: np.ndarray, figsize=(9, 5)) -> plt.Figure:
    """
    Plot Lyapunov exponent convergence.

    A positive converging value → chaotic system.
    A value converging to 0 → regular/periodic orbit.

    This is the same Lyapunov exponent concept used in black hole
    photon ring theory (Gralla & Lupsasca 2020).
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(times, lyapunov, color='#E8632A', lw=1.5, label='λ_L(t)')
    ax.axhline(lyapunov[-1], color='#2AE87A', ls='--', lw=1,
               label=f'Converged: λ_L = {lyapunov[-1]:.4f} yr⁻¹')
    ax.axhline(0, color='#8B949E', ls=':', lw=0.8)

    ax.set_xlabel('Time [yr]')
    ax.set_ylabel('Lyapunov Exponent λ_L [yr⁻¹]')
    ax.set_title('Maximal Lyapunov Exponent — Convergence')
    ax.legend()
    ax.grid(True)

    # Annotation box
    ax.text(0.98, 0.05,
            'λ_L > 0 → chaotic\nλ_L = 0 → regular',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=9, color='#8B949E',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#21262D', alpha=0.7))
    fig.tight_layout()
    return fig


def animate_trajectories(
    result: SimulationResult,
    output_path: Optional[str] = None,
    fps: int = 30,
    trail_length: int = 200,
    figsize: Tuple[int, int] = (8, 8),
    interval: int = 20,
) -> animation.FuncAnimation:
    """
    Create an animated movie of the three-body orbits.

    Parameters
    ----------
    result : SimulationResult
    output_path : str, optional
        Save to file if given (e.g., 'orbit.mp4' or 'orbit.gif')
    fps : int
        Frames per second for saved animation
    trail_length : int
        Number of past positions to show as trail
    figsize : tuple
    interval : int
        Delay between frames in milliseconds

    Returns
    -------
    anim : matplotlib.animation.FuncAnimation
    """
    pos = result.positions
    N = result.n_steps

    # Compute axis limits
    margin = 0.2
    all_x = pos[:, :, 0].flatten()
    all_y = pos[:, :, 1].flatten()
    xlim = (all_x.min() - margin * abs(all_x.min()),
            all_x.max() + margin * abs(all_x.max()))
    ylim = (all_y.min() - margin * abs(all_y.min()),
            all_y.max() + margin * abs(all_y.max()))

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect('equal')
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.grid(True, alpha=0.2)

    # Initialize artists
    trails = [ax.plot([], [], '-', color=c, alpha=0.5, lw=1.0)[0]
              for c in BODY_COLORS]
    points = [ax.plot([], [], 'o', color=c, ms=8, zorder=5,
                      markeredgecolor='white', markeredgewidth=0.5)[0]
              for c in BODY_COLORS]
    time_text = ax.text(0.02, 0.97, '', transform=ax.transAxes,
                        fontsize=10, va='top', family='monospace')
    energy_text = ax.text(0.02, 0.90, '', transform=ax.transAxes,
                          fontsize=8, va='top', color='#8B949E', family='monospace')

    for b, label in enumerate(BODY_LABELS):
        ax.plot([], [], 'o', color=BODY_COLORS[b], label=label)
    ax.legend(loc='upper right', fontsize=9)

    def init():
        for t in trails: t.set_data([], [])
        for p in points: p.set_data([], [])
        time_text.set_text('')
        return trails + points + [time_text]

    def update(frame):
        start = max(0, frame - trail_length)
        for b in range(3):
            trails[b].set_data(pos[start:frame+1, b, 0], pos[start:frame+1, b, 1])
            points[b].set_data([pos[frame, b, 0]], [pos[frame, b, 1]])
        time_text.set_text(f't = {result.t[frame]:.2f} yr')
        energy_text.set_text(f'ΔE/E₀ = {result.energy_error[frame]:.1e}')
        return trails + points + [time_text, energy_text]

    skip = max(1, N // (fps * 30))   # limit to ~30 seconds of animation
    frames = range(0, N, skip)

    anim = animation.FuncAnimation(
        fig, update, frames=frames,
        init_func=init, blit=True, interval=interval,
    )

    if output_path:
        if output_path.endswith('.gif'):
            anim.save(output_path, writer='pillow', fps=fps)
        else:
            anim.save(output_path, writer='ffmpeg', fps=fps,
                      extra_args=['-vcodec', 'libx264'])
        print(f"Animation saved to {output_path}")

    return anim


def plot_all(result: SimulationResult, name: str = "simulation") -> None:
    """
    Generate and save all standard plots for a simulation result.

    Saves to docs/figures/{name}_*.png
    """
    import os
    os.makedirs('docs/figures', exist_ok=True)

    figs = {
        'trajectories': plot_trajectories(result),
        'energy':       plot_energy(result),
        'separation':   plot_separation(result),
        'phase_space':  plot_phase_space(result),
    }

    for tag, fig in figs.items():
        path = f'docs/figures/{name}_{tag}.png'
        fig.savefig(path, dpi=150, bbox_inches='tight', facecolor='#0D1117')
        print(f"Saved: {path}")
        plt.close(fig)