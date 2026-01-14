"""
Soliton Core Profile Fitting for ULDM Simulations

Fits theoretical soliton density profile to PyUltraLight simulation data
and compares with theoretical predictions.

Theoretical profile: ρ(r) = ρ_c / [1 + (r/r_c)²]⁴

Date: 2025-12-22
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
import json
from scipy.optimize import curve_fit
import sys
import os

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(__file__))
from utils import load_pyul_snapshot, spherical_average, soliton_profile, nfw_profile


# Set publication-quality plot parameters
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.format': 'pdf',
    'font.family': 'serif',
    'text.usetex': False,
    'lines.linewidth': 2,
    'axes.grid': True,
    'grid.alpha': 0.3
})


def fit_soliton(r, rho, r_max_fit=None):
    """
    Fit soliton profile to radial density data.

    Parameters
    ----------
    r : ndarray
        Radial bins (physical units)
    rho : ndarray
        Density profile
    r_max_fit : float, optional
        Maximum radius for fitting. Default is half the box.

    Returns
    -------
    popt : tuple
        (rho_c, r_c) best-fit parameters
    pcov : ndarray
        Covariance matrix
    fit_profile : ndarray
        Best-fit profile
    """
    if r_max_fit is None:
        r_max_fit = r.max() / 4  # Fit inner quarter

    # Select fitting region
    mask = (r > 0) & (r < r_max_fit) & (rho > 0)

    r_fit = r[mask]
    rho_fit = rho[mask]

    # Initial guess
    rho_c_guess = rho_fit[0] if len(rho_fit) > 0 else 1.0
    r_c_guess = r_fit[np.argmax(rho_fit < rho_c_guess / 2)] if len(r_fit) > 1 else 0.01

    try:
        popt, pcov = curve_fit(
            soliton_profile,
            r_fit,
            rho_fit,
            p0=[rho_c_guess, r_c_guess],
            bounds=([0, 0], [np.inf, r_max_fit]),
            maxfev=10000
        )

        # Compute fit over full range
        fit_profile = soliton_profile(r, *popt)

        return popt, pcov, fit_profile

    except Exception as e:
        print(f"Fitting failed: {e}")
        return None, None, None


def compute_fit_quality(rho_data, rho_fit, mask):
    """
    Compute goodness-of-fit metrics.

    Parameters
    ----------
    rho_data : ndarray
        Observed density
    rho_fit : ndarray
        Fitted density
    mask : ndarray
        Boolean mask for fitting region

    Returns
    -------
    chi2_red : float
        Reduced chi-squared
    rmse : float
        Root mean squared error
    """
    residuals = rho_data[mask] - rho_fit[mask]
    n_points = np.sum(mask)
    n_params = 2  # rho_c, r_c

    # Chi-squared (assuming Poisson errors: σ² ~ ρ)
    sigma = np.sqrt(np.abs(rho_data[mask]))
    sigma[sigma == 0] = 1  # Avoid division by zero

    chi2 = np.sum((residuals / sigma)**2)
    chi2_red = chi2 / (n_points - n_params) if n_points > n_params else np.inf

    # RMSE
    rmse = np.sqrt(np.mean(residuals**2))

    return chi2_red, rmse


def analyze_snapshot(snap_path, box_size_kpc=20.0, output_dir=None):
    """
    Analyze a single PyUltraLight snapshot.

    Parameters
    ----------
    snap_path : str
        Path to snapshot file
    box_size_kpc : float
        Box size in kpc
    output_dir : str, optional
        Output directory for results

    Returns
    -------
    results : dict
        Analysis results
    """
    print(f"\nAnalyzing {os.path.basename(snap_path)}...")

    # Load data
    density, phase, metadata = load_pyul_snapshot(snap_path)
    N = density.shape[0]
    dx_kpc = box_size_kpc / N

    print(f"  Grid: {N}³")
    print(f"  Resolution: {dx_kpc:.3f} kpc/cell")
    print(f"  Mean density: {np.mean(density):.3f}")
    print(f"  Max density: {np.max(density):.3f}")

    # Compute radial profile
    r_bins, rho_profile = spherical_average(density)
    r_kpc = r_bins * dx_kpc

    # Fit soliton profile
    print("  Fitting soliton profile...")
    popt, pcov, fit_profile = fit_soliton(r_kpc, rho_profile, r_max_fit=box_size_kpc/4)

    if popt is None:
        print("  WARNING: Fit failed!")
        return None

    rho_c, r_c = popt
    rho_c_err, r_c_err = np.sqrt(np.diag(pcov))

    print(f"  ρ_c = {rho_c:.2f} ± {rho_c_err:.2f} (code units)")
    print(f"  r_c = {r_c:.4f} ± {r_c_err:.4f} kpc")

    # Compute fit quality
    mask = (r_kpc > 0) & (r_kpc < box_size_kpc/4)
    chi2_red, rmse = compute_fit_quality(rho_profile, fit_profile, mask)

    print(f"  χ²_red = {chi2_red:.3f}")
    print(f"  RMSE = {rmse:.3f}")

    # Store results
    results = {
        'snapshot': os.path.basename(snap_path),
        'N': int(N),
        'box_size_kpc': float(box_size_kpc),
        'dx_kpc': float(dx_kpc),
        'mean_density': float(np.mean(density)),
        'max_density': float(np.max(density)),
        'rho_c': float(rho_c),
        'rho_c_err': float(rho_c_err),
        'r_c_kpc': float(r_c),
        'r_c_err_kpc': float(r_c_err),
        'chi2_red': float(chi2_red),
        'rmse': float(rmse),
        'r_profile_kpc': r_kpc.tolist(),
        'rho_profile': rho_profile.tolist(),
        'fit_profile': fit_profile.tolist()
    }

    # Save results
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        json_path = os.path.join(output_dir, f"fit_{os.path.basename(snap_path).replace('.h5', '.json')}")
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"  Saved results to {json_path}")

    return results


def plot_fit_comparison(results, output_path=None):
    """
    Create publication-quality soliton fit plot.

    Parameters
    ----------
    results : dict
        Fit results from analyze_snapshot()
    output_path : str, optional
        Output path for figure
    """
    r = np.array(results['r_profile_kpc'])
    rho_data = np.array(results['rho_profile'])
    rho_fit = np.array(results['fit_profile'])

    rho_c = results['rho_c']
    r_c = results['r_c_kpc']
    chi2 = results['chi2_red']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True,
                                     gridspec_kw={'height_ratios': [3, 1]})

    # Main plot
    ax1.loglog(r, rho_data, 'o', markersize=4, alpha=0.6, label='Simulation data')
    ax1.loglog(r, rho_fit, 'r-', linewidth=2, label='Soliton fit')

    # Add NFW for comparison (arbitrary normalization)
    rho_s_nfw = rho_c / 10
    r_s_nfw = r_c * 5
    rho_nfw = nfw_profile(r, rho_s_nfw, r_s_nfw)
    ax1.loglog(r, rho_nfw, 'k--', linewidth=1.5, alpha=0.7, label='NFW (for comparison)')

    ax1.set_ylabel(r'Density $\rho$ (code units)', fontsize=16)
    ax1.legend(loc='best', fontsize=12)
    ax1.set_title(f'Soliton Core Profile: {results["snapshot"]}', fontsize=16)

    # Add text box with fit parameters
    textstr = f'$\\rho_c = {rho_c:.2f} ± {results["rho_c_err"]:.2f}$\n'
    textstr += f'$r_c = {r_c:.3f} ± {results["r_c_err_kpc"]:.3f}$ kpc\n'
    textstr += f'$\\chi^2_{{\\mathrm{{red}}}} = {chi2:.2f}$'

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.05, 0.05, textstr, transform=ax1.transAxes, fontsize=12,
             verticalalignment='bottom', bbox=props)

    # Residuals
    mask = (r > 0) & (rho_data > 0) & (rho_fit > 0)
    residuals = (rho_data[mask] - rho_fit[mask]) / rho_fit[mask] * 100  # Percent

    ax2.semilogx(r[mask], residuals, 'o', markersize=3, alpha=0.5)
    ax2.axhline(0, color='r', linestyle='--', linewidth=1)
    ax2.fill_between(r[mask], -10, 10, alpha=0.2, color='gray')

    ax2.set_xlabel('Radius $r$ (kpc)', fontsize=16)
    ax2.set_ylabel('Residuals (%)', fontsize=14)
    ax2.set_ylim(-50, 50)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {output_path}")
    else:
        plt.show()

    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Fit soliton profile to ULDM simulation data')
    parser.add_argument('snapshot', help='Path to PyUltraLight snapshot (e.g., snap_0100.h5)')
    parser.add_argument('--box-size', type=float, default=20.0, help='Box size in kpc (default: 20)')
    parser.add_argument('--output-dir', default='../results/data', help='Output directory for results')
    parser.add_argument('--figure-dir', default='../results/figures', help='Output directory for figures')
    parser.add_argument('--no-plot', action='store_true', help='Skip plotting')

    args = parser.parse_args()

    # Analyze snapshot
    results = analyze_snapshot(args.snapshot, box_size_kpc=args.box_size, output_dir=args.output_dir)

    if results is None:
        print("Analysis failed!")
        return 1

    # Plot
    if not args.no_plot:
        os.makedirs(args.figure_dir, exist_ok=True)
        snap_name = os.path.basename(args.snapshot).replace('.h5', '')
        fig_path = os.path.join(args.figure_dir, f'soliton_fit_{snap_name}.pdf')
        plot_fit_comparison(results, output_path=fig_path)

    print("\n✓ Analysis complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
