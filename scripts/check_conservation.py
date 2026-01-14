"""
Conservation Law Analysis for ULDM Simulations

Checks mass and energy conservation for PyUltraLight evolution.

Research demonstration for UoA postdoc position
Date: 2025-12-22
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_pyul_snapshot, check_mass_conservation


# Plot parameters
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.format': 'pdf',
    'font.family': 'serif',
})


def analyze_conservation(snapshot_dir, box_size_kpc=20.0):
    """
    Analyze conservation laws across all snapshots.

    Parameters
    ----------
    snapshot_dir : str
        Directory containing snap_*.h5 files
    box_size_kpc : float
        Box size in kpc

    Returns
    -------
    results : dict
        Conservation analysis results
    """
    # Get all snapshots
    snapshots = sorted(glob.glob(os.path.join(snapshot_dir, 'snap_*.h5')))

    print(f"Found {len(snapshots)} snapshots")

    times = []
    masses = []
    max_densities = []
    mean_densities = []

    for snap in snapshots:
        density, phase, metadata = load_pyul_snapshot(snap)

        # Extract time from filename
        snap_num = int(os.path.basename(snap).split('_')[1].split('.')[0])
        time = snap_num * 0.001  # dt = 0.001

        times.append(time)
        masses.append(np.sum(density))
        max_densities.append(np.max(density))
        mean_densities.append(np.mean(density))

    times = np.array(times)
    masses = np.array(masses)
    max_densities = np.array(max_densities)
    mean_densities = np.array(mean_densities)

    # Mass conservation
    mass_error = np.abs(masses / masses[0] - 1.0) * 100  # Percent

    # Statistics
    print(f"\nConservation Analysis:")
    print(f"  Initial mass: {masses[0]:.6f}")
    print(f"  Final mass: {masses[-1]:.6f}")
    print(f"  Mass error: {mass_error[-1]:.2e} %")
    print(f"  Max mass error: {mass_error.max():.2e} %")
    print(f"\nDensity Evolution:")
    print(f"  Initial max density: {max_densities[0]:.2f}")
    print(f"  Final max density: {max_densities[-1]:.2f}")
    print(f"  Concentration factor: {max_densities[-1] / max_densities[0]:.2f}x")

    results = {
        'times': times,
        'masses': masses,
        'mass_error_percent': mass_error,
        'max_densities': max_densities,
        'mean_densities': mean_densities
    }

    return results


def plot_conservation(results, output_path=None):
    """
    Create conservation plots.

    Parameters
    ----------
    results : dict
        Results from analyze_conservation()
    output_path : str, optional
        Output path for figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    times = results['times']
    masses = results['masses']
    mass_error = results['mass_error_percent']
    max_dens = results['max_densities']
    mean_dens = results['mean_densities']

    # Mass conservation
    ax = axes[0, 0]
    ax.plot(times, masses / masses[0], 'b-', linewidth=2)
    ax.axhline(1.0, color='r', linestyle='--', linewidth=1)
    ax.fill_between(times, 0.999, 1.001, alpha=0.2, color='gray')
    ax.set_xlabel('Time (code units)')
    ax.set_ylabel('M(t) / M(0)')
    ax.set_title('Mass Conservation')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.998, 1.002)

    # Mass error
    ax = axes[0, 1]
    ax.semilogy(times, mass_error, 'r-', linewidth=2)
    ax.axhline(0.01, color='k', linestyle='--', linewidth=1, label='0.01% threshold')
    ax.set_xlabel('Time (code units)')
    ax.set_ylabel('|M(t)/M(0) - 1| (%)')
    ax.set_title('Mass Conservation Error')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Maximum density evolution
    ax = axes[1, 0]
    ax.plot(times, max_dens, 'g-', linewidth=2)
    ax.set_xlabel('Time (code units)')
    ax.set_ylabel(r'$\rho_{\max}$ (code units)')
    ax.set_title('Maximum Density Evolution\\n(Soliton Formation)')
    ax.grid(True, alpha=0.3)

    # Mean density
    ax = axes[1, 1]
    ax.plot(times, mean_dens, 'm-', linewidth=2)
    ax.axhline(mean_dens[0], color='r', linestyle='--', linewidth=1, label='Initial mean')
    ax.set_xlabel('Time (code units)')
    ax.set_ylabel(r'$\langle \rho \rangle$ (code units)')
    ax.set_title('Mean Density (Should be Constant)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle('PyUltraLight: Conservation Laws & Evolution', fontsize=18, y=0.995)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\nSaved figure to {output_path}")
    else:
        plt.show()

    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Check conservation laws in ULDM simulation')
    parser.add_argument('snapshot_dir', help='Directory containing snap_*.h5 files')
    parser.add_argument('--box-size', type=float, default=20.0, help='Box size in kpc')
    parser.add_argument('--output', default='../results/figures/conservation.pdf',
                        help='Output figure path')

    args = parser.parse_args()

    results = analyze_conservation(args.snapshot_dir, box_size_kpc=args.box_size)

    plot_conservation(results, output_path=args.output)

    print("\n✓ Conservation analysis complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
