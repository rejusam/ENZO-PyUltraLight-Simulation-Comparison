"""
Visualize ULDM Soliton Formation Evolution

Creates publication-quality figures showing density evolution
and soliton core formation over time.

Date: 2025-12-22
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_pyul_snapshot, spherical_average


# Plot parameters
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.format': 'pdf',
    'font.family': 'serif',
})


def create_evolution_figure(snapshot_dir, snapshots_to_plot=[0, 30, 60, 100],
                            box_size_kpc=20.0, output_path=None):
    """
    Create multi-panel figure showing evolution.

    Parameters
    ----------
    snapshot_dir : str
        Directory with snapshots
    snapshots_to_plot : list
        Snapshot numbers to include
    box_size_kpc : float
        Box size in kpc
    output_path : str, optional
        Output file path
    """
    n_snaps = len(snapshots_to_plot)

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, n_snaps, hspace=0.3, wspace=0.25)

    for idx, snap_num in enumerate(snapshots_to_plot):
        snap_path = os.path.join(snapshot_dir, f'snap_{snap_num:04d}.h5')

        if not os.path.exists(snap_path):
            print(f"Warning: {snap_path} not found")
            continue

        density, phase, metadata = load_pyul_snapshot(snap_path)
        N = density.shape[0]
        dx_kpc = box_size_kpc / N

        time = snap_num * 0.001  # dt = 0.001

        # Row 1: Density slice (log scale)
        ax = fig.add_subplot(gs[0, idx])
        mid = N // 2
        im = ax.imshow(np.log10(density[mid, :, :] + 1e-3),
                      cmap='inferno', origin='lower',
                      extent=[0, box_size_kpc, 0, box_size_kpc])
        ax.set_title(f't = {time:.3f}\\n(Step {snap_num})', fontsize=12)
        if idx == 0:
            ax.set_ylabel('Density Slice\\n(log scale, kpc)', fontsize=12)
        ax.set_xlabel('kpc', fontsize=10)
        plt.colorbar(im, ax=ax, label=r'$\log_{10} \rho$')

        # Row 2: Phase slice
        ax = fig.add_subplot(gs[1, idx])
        im = ax.imshow(phase[mid, :, :],
                      cmap='twilight', origin='lower',
                      extent=[0, box_size_kpc, 0, box_size_kpc],
                      vmin=-np.pi, vmax=np.pi)
        if idx == 0:
            ax.set_ylabel('Phase Field\\n(radians, kpc)', fontsize=12)
        ax.set_xlabel('kpc', fontsize=10)
        plt.colorbar(im, ax=ax, label=r'$\arg(\psi)$')

        # Row 3: Radial profile
        ax = fig.add_subplot(gs[2, idx])
        r_bins, rho_profile = spherical_average(density)
        r_kpc = r_bins * dx_kpc

        ax.loglog(r_kpc[1:], rho_profile[1:], 'b-', linewidth=2)
        ax.set_xlabel('Radius (kpc)', fontsize=10)
        if idx == 0:
            ax.set_ylabel(r'$\rho(r)$ (code units)', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0.1, box_size_kpc/2)

        # Add statistics
        stats_text = f'Max: {np.max(density):.0f}\\nMean: {np.mean(density):.2f}'
        ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7),
               fontsize=9)

    plt.suptitle('PyUltraLight: Soliton Formation in ULDM Halo (m = $10^{-22}$ eV)',
                 fontsize=16, y=0.995)

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {output_path}")
    else:
        plt.show()

    plt.close()


def create_radial_evolution_figure(snapshot_dir, box_size_kpc=20.0, output_path=None):
    """
    Create figure showing radial profile evolution over time.

    Parameters
    ----------
    snapshot_dir : str
        Directory with snapshots
    box_size_kpc : float
        Box size
    output_path : str, optional
        Output path
    """
    snapshots = sorted(glob.glob(os.path.join(snapshot_dir, 'snap_*.h5')))

    # Select subset for clarity
    snap_indices = [0, 10, 30, 50, 70, 100]
    colors = plt.cm.plasma(np.linspace(0, 1, len(snap_indices)))

    fig, ax = plt.subplots(figsize=(10, 8))

    for i, snap_idx in enumerate(snap_indices):
        if snap_idx >= len(snapshots):
            continue

        snap_path = snapshots[snap_idx]
        density, _, _ = load_pyul_snapshot(snap_path)

        N = density.shape[0]
        dx_kpc = box_size_kpc / N

        r_bins, rho_profile = spherical_average(density)
        r_kpc = r_bins * dx_kpc

        time = snap_idx * 0.001
        label = f't = {time:.3f}'

        ax.loglog(r_kpc[1:], rho_profile[1:], '-', color=colors[i],
                 linewidth=2, label=label, alpha=0.8)

    ax.set_xlabel('Radius (kpc)', fontsize=14)
    ax.set_ylabel(r'Density $\rho(r)$ (code units)', fontsize=14)
    ax.set_title('Radial Density Profile Evolution\\nSoliton Core Formation',
                fontsize=16)
    ax.legend(loc='lower left', fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(0.1, box_size_kpc/2)

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved radial evolution to {output_path}")
    else:
        plt.show()

    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Visualize ULDM soliton evolution')
    parser.add_argument('snapshot_dir', help='Directory containing snap_*.h5 files')
    parser.add_argument('--box-size', type=float, default=20.0, help='Box size in kpc')
    parser.add_argument('--output-multi', default='../results/figures/evolution_multi.pdf',
                       help='Output path for multi-panel figure')
    parser.add_argument('--output-radial', default='../results/figures/evolution_radial.pdf',
                       help='Output path for radial evolution figure')

    args = parser.parse_args()

    print("Creating multi-panel evolution figure...")
    create_evolution_figure(args.snapshot_dir, box_size_kpc=args.box_size,
                           output_path=args.output_multi)

    print("\\nCreating radial evolution figure...")
    create_radial_evolution_figure(args.snapshot_dir, box_size_kpc=args.box_size,
                                   output_path=args.output_radial)

    print("\\n✓ Visualization complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
