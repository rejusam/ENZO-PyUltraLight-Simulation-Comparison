"""
Utility Functions for ULDM Simulation Analysis

This module provides shared functionality for analyzing ultralight dark matter
simulations from Enzo and PyUltraLight codes.

Date: 2025-12-22
"""

import numpy as np
import h5py
import os
from typing import Tuple, Dict, Optional


# Physical Constants (CGS)
G_CGS = 6.674e-8          # Gravitational constant
HBAR_CGS = 1.05457e-27    # Reduced Planck constant
EV_TO_ERG = 1.60218e-12   # Electron-volt to erg conversion
C_CGS = 2.99792458e10     # Speed of light
KPC_TO_CM = 3.08567758e21 # Kiloparsec to cm
MSUN_TO_G = 1.989e33      # Solar mass to grams


def load_pyul_snapshot(snap_path: str) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Load a PyUltraLight simulation snapshot.

    Parameters
    ----------
    snap_path : str
        Path to HDF5 snapshot file

    Returns
    -------
    density : ndarray
        3D density field |ψ|²
    phase : ndarray
        3D phase field arg(ψ)
    metadata : dict
        Snapshot metadata (time, step, etc.)
    """
    with h5py.File(snap_path, 'r') as f:
        density = f['density'][()]
        phase = f['phase'][()]

        # Extract metadata if available
        metadata = {
            'file': os.path.basename(snap_path),
            'shape': density.shape
        }

        # Try to extract attributes
        for key in f.attrs.keys():
            metadata[key] = f.attrs[key]

    return density, phase, metadata


def load_enzo_snapshot(dd_path: str) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Load an Enzo FDM simulation snapshot.

    Parameters
    ----------
    dd_path : str
        Path to Enzo DD directory (e.g., "DD0001")

    Returns
    -------
    density : ndarray
        3D FDM density field
    wavefunction : complex ndarray
        3D complex wavefunction ψ = Re(ψ) + i*Im(ψ)
    metadata : dict
        Snapshot metadata
    """
    # Find the .cpu0000 file
    cpu_file = os.path.join(dd_path, f"{os.path.basename(dd_path)}.cpu0000")

    with h5py.File(cpu_file, 'r') as f:
        # Get grid key
        grid_key = list(f.keys())[0]
        ds = f[grid_key]

        # Load FDM fields
        if "FDMDensity" in ds:
            density = ds["FDMDensity"][()]
        else:
            # Compute from wavefunction
            re_psi = ds["Re_Psi"][()]
            im_psi = ds["Im_Psi"][()]
            density = re_psi**2 + im_psi**2

        # Get wavefunction
        re_psi = ds["Re_Psi"][()]
        im_psi = ds["Im_Psi"][()]
        wavefunction = re_psi + 1j * im_psi

        # Metadata
        metadata = {
            'file': cpu_file,
            'shape': density.shape
        }

        # Try to get time from attributes
        for key in ds.attrs.keys():
            metadata[key] = ds.attrs[key]

    return density, wavefunction, metadata


def spherical_average(data: np.ndarray, center: Optional[Tuple] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute spherically-averaged radial profile.

    Parameters
    ----------
    data : ndarray
        3D data field to average
    center : tuple, optional
        (i, j, k) center coordinates. Default is box center.

    Returns
    -------
    r_bins : ndarray
        Radial bin centers (in grid units)
    profile : ndarray
        Spherically-averaged profile
    """
    if center is None:
        center = tuple(s // 2 for s in data.shape)

    # Create radial distance grid
    indices = np.indices(data.shape)
    r = np.sqrt(sum((idx - c)**2 for idx, c in zip(indices, center)))
    r_int = r.astype(int)

    # Bin by radius
    profile = np.bincount(r_int.ravel(), weights=data.ravel())
    counts = np.bincount(r_int.ravel())
    counts[counts == 0] = 1  # Avoid division by zero

    r_bins = np.arange(len(profile))
    profile = profile / counts

    return r_bins, profile


def radial_profile_with_errors(data: np.ndarray, center: Optional[Tuple] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute spherically-averaged radial profile with standard deviation.

    Parameters
    ----------
    data : ndarray
        3D data field
    center : tuple, optional
        Center coordinates

    Returns
    -------
    r_bins : ndarray
        Radial bins
    profile : ndarray
        Mean profile
    std : ndarray
        Standard deviation in each bin
    """
    if center is None:
        center = tuple(s // 2 for s in data.shape)

    indices = np.indices(data.shape)
    r = np.sqrt(sum((idx - c)**2 for idx, c in zip(indices, center)))
    r_int = r.astype(int)

    # Mean
    profile = np.bincount(r_int.ravel(), weights=data.ravel())
    counts = np.bincount(r_int.ravel())
    counts[counts == 0] = 1
    profile = profile / counts

    # Variance
    variance = np.bincount(r_int.ravel(), weights=(data.ravel() - profile[r_int.ravel()])**2)
    std = np.sqrt(variance / counts)

    r_bins = np.arange(len(profile))

    return r_bins, profile, std


def compute_mass_profile(density: np.ndarray, dx: float, center: Optional[Tuple] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute enclosed mass profile M(r).

    Parameters
    ----------
    density : ndarray
        3D density field
    dx : float
        Grid spacing (physical units)
    center : tuple, optional
        Center coordinates

    Returns
    -------
    r : ndarray
        Radius array (physical units)
    M_enc : ndarray
        Enclosed mass M(<r)
    """
    r_bins, rho_profile = spherical_average(density, center)

    # Convert to physical radius
    r = r_bins * dx

    # Compute shell volumes: V_shell = 4π r² dr
    dr = np.diff(r)
    dr = np.append(dr, dr[-1])  # Extend for last bin

    shell_volumes = 4 * np.pi * r**2 * dr
    shell_masses = rho_profile * shell_volumes

    # Cumulative sum
    M_enc = np.cumsum(shell_masses)

    return r, M_enc


def get_snapshot_time(snap_path: str, code: str) -> float:
    """
    Extract simulation time from snapshot.

    Parameters
    ----------
    snap_path : str
        Path to snapshot
    code : str
        'pyul' or 'enzo'

    Returns
    -------
    time : float
        Simulation time (code units)
    """
    if code.lower() == 'pyul':
        # Extract from filename: snap_0050.h5 -> step 50
        basename = os.path.basename(snap_path)
        step = int(basename.split('_')[1].split('.')[0])
        dt = 0.001  # Hardcoded from PyUltraLight setup
        return step * dt

    elif code.lower() == 'enzo':
        # Try to read from hierarchy file
        dd_dir = os.path.dirname(snap_path)
        hierarchy_file = os.path.join(dd_dir, f"{os.path.basename(dd_dir)}.hierarchy")

        # For now, return placeholder
        # In practice, would parse hierarchy file
        return 0.0

    else:
        raise ValueError(f"Unknown code: {code}")


def soliton_profile(r: np.ndarray, rho_c: float, r_c: float) -> np.ndarray:
    """
    Theoretical soliton density profile.

    ρ(r) = ρ_c / [1 + (r/r_c)²]⁴

    Parameters
    ----------
    r : ndarray
        Radius array
    rho_c : float
        Central density
    r_c : float
        Core radius

    Returns
    -------
    rho : ndarray
        Density profile
    """
    return rho_c / (1 + (r / r_c)**2)**4


def nfw_profile(r: np.ndarray, rho_s: float, r_s: float) -> np.ndarray:
    """
    NFW (Navarro-Frenk-White) density profile.

    ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]

    Parameters
    ----------
    r : ndarray
        Radius array
    rho_s : float
        Scale density
    r_s : float
        Scale radius

    Returns
    -------
    rho : ndarray
        Density profile
    """
    x = r / r_s
    x[x == 0] = 1e-10  # Avoid division by zero
    return rho_s / (x * (1 + x)**2)


def compute_power_spectrum(density: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute 3D power spectrum P(k) from density field.

    Parameters
    ----------
    density : ndarray
        3D density field

    Returns
    -------
    k_bins : ndarray
        Wavenumber bins
    P_k : ndarray
        Power spectrum P(k)
    """
    N = density.shape[0]
    assert density.shape == (N, N, N), "Only cubic grids supported"

    # FFT
    delta_k = np.fft.fftn(density - np.mean(density))
    power_3d = np.abs(delta_k)**2

    # Create k-space grid
    kx = 2 * np.pi * np.fft.fftfreq(N)
    ky = 2 * np.pi * np.fft.fftfreq(N)
    kz = 2 * np.pi * np.fft.fftfreq(N)

    kx_grid, ky_grid, kz_grid = np.meshgrid(kx, ky, kz, indexing='ij')
    k_mag = np.sqrt(kx_grid**2 + ky_grid**2 + kz_grid**2)

    # Spherical binning
    k_int = k_mag.astype(int)
    k_max = k_int.max()

    k_bins = np.arange(k_max + 1)
    P_k = np.bincount(k_int.ravel(), weights=power_3d.ravel())
    counts = np.bincount(k_int.ravel())
    counts[counts == 0] = 1

    P_k = P_k / counts / N**6  # Normalization

    return k_bins, P_k


def check_mass_conservation(snapshots: list, code: str = 'pyul') -> Tuple[np.ndarray, np.ndarray]:
    """
    Check mass conservation across multiple snapshots.

    Parameters
    ----------
    snapshots : list
        List of snapshot paths
    code : str
        'pyul' or 'enzo'

    Returns
    -------
    times : ndarray
        Snapshot times
    mass_ratio : ndarray
        M(t) / M(0) - 1 (should be ~0)
    """
    masses = []
    times = []

    for snap in snapshots:
        if code == 'pyul':
            density, _, _ = load_pyul_snapshot(snap)
        else:
            density, _, _ = load_enzo_snapshot(snap)

        mass = np.sum(density)
        masses.append(mass)
        times.append(get_snapshot_time(snap, code))

    masses = np.array(masses)
    times = np.array(times)

    mass_ratio = masses / masses[0] - 1.0

    return times, mass_ratio


def compute_velocity_field(wavefunction: np.ndarray, dx: float, hbar_over_m: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute velocity field from wavefunction.

    v = (ℏ/m) * ∇φ / ρ

    where φ = arg(ψ)

    Parameters
    ----------
    wavefunction : complex ndarray
        Complex wavefunction ψ
    dx : float
        Grid spacing
    hbar_over_m : float
        ℏ/m coefficient

    Returns
    -------
    vx, vy, vz : ndarray
        Velocity components
    """
    rho = np.abs(wavefunction)**2
    phase = np.angle(wavefunction)

    # Gradient of phase
    grad_phi_x = np.gradient(phase, dx, axis=0)
    grad_phi_y = np.gradient(phase, dx, axis=1)
    grad_phi_z = np.gradient(phase, dx, axis=2)

    # Avoid division by zero
    rho_safe = rho.copy()
    rho_safe[rho_safe < 1e-10] = 1e-10

    vx = hbar_over_m * grad_phi_x / rho_safe
    vy = hbar_over_m * grad_phi_y / rho_safe
    vz = hbar_over_m * grad_phi_z / rho_safe

    return vx, vy, vz


if __name__ == "__main__":
    print("ULDM Simulation Utilities Module")
    print("=" * 50)
    print(f"Physical constants loaded:")
    print(f"  G = {G_CGS:.3e} cm³/g/s²")
    print(f"  ℏ = {HBAR_CGS:.3e} erg·s")
    print(f"  c = {C_CGS:.3e} cm/s")
    print("\nAvailable functions:")
    print("  - load_pyul_snapshot()")
    print("  - load_enzo_snapshot()")
    print("  - spherical_average()")
    print("  - soliton_profile()")
    print("  - nfw_profile()")
    print("  - compute_power_spectrum()")
    print("  - check_mass_conservation()")
    print("  - compute_velocity_field()")
