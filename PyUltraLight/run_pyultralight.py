import numpy as np
import h5py
import time
import os

# Physical Constants (CGS)
G_CGS = 6.674e-8
HBAR_CGS = 1.05457e-27
EV_TO_ERG = 1.60218e-12

def run_simulation():
    print("Initializing PyUltraLight v2...")
    
    # Load ICs
    with h5py.File('ic_uldm_dwarf.h5', 'r') as f:
        re_psi = f['psi_re'][()]
        im_psi = f['psi_im'][()]
        N = f.attrs['N']
        L_kpc = f.attrs['L']
        m_ev = f.attrs['m_axion_eV']
        L_unit = f.attrs['units_L']
        T_unit = f.attrs['units_T']
        D_unit = f.attrs['units_D']

    psi = re_psi + 1j * im_psi
    
    # Simulation Parameters
    dt = 0.001 # Code units
    n_steps = 100
    dump_interval = 10
    
    # Grid setup (Code Units: L=1, T=1, G=1? No, sticking to physics-scaled code units)
    # Actually, simpler to run in dimensionless units where G=1, hbar_code = ...
    # But let's stick to the generated units structure.
    # Code constants:
    # G_code = G_cgs * D_unit * T_unit**2 
    # But D_unit = M / L^3 and T = sqrt(L^3/GM). So G_code = G * (M/L^3) * (L^3/GM) = 1.
    G_code = 1.0
    
    # hbar_code = hbar_cgs / (M_unit * L_unit^2 / T_unit) ? 
    # Actually, Schrodinger eq: i dpsi/dt = -1/2 alpha laplacian psi + V psi
    # alpha = hbar / m. 
    # In code units: alpha_code = (hbar/m)_cgs / (L_unit^2 / T_unit)
    
    m_cgs = m_ev * EV_TO_ERG / (2.998e10**2) # E = mc^2 -> m = E/c^2 !! Correction here
    # Wait, EV_TO_ERG is energy. mass in grams = (eV * EV_TO_ERG) / c^2
    c_cgs = 2.99792458e10
    m_gram = (m_ev * EV_TO_ERG) / c_cgs**2
    
    alpha_cgs = HBAR_CGS / m_gram
    alpha_code = alpha_cgs / (L_unit**2 / T_unit)
    
    print(f"Alpha (code): {alpha_code:.2e}")
    
    # K-space
    dx = 1.0 / N
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
    k2 = kx**2 + ky**2 + kz**2
    k2[0,0,0] = 1e-10 # Avoid div zero for potential
    
    # Kinetic Operator
    # Operator: exp(-i * (alpha/2) * k^2 * dt)
    # Actually SE: i dpsi/dt = -alpha/2 del^2 psi.  FFT: i dpsi_k/dt = alpha/2 k^2 psi_k
    # dpsi_k/dt = -i alpha/2 k^2 psi_k -> psi(t) = exp(-i alpha/2 k^2 t) psi(0)
    # But wait, standard dimensionless form often uses alpha=1. Here alpha is huge?
    # Let's check T_unit. T_unit ~ 4e16 s ~ 1 Gyr.
    # Oscillation timescale ~ 1/m. very fast.
    # The Schrodinger-Poisson system is often rescaled so alpha ~ 1.
    # Current units might be stiff. 
    # Let's proceed, but note that dt must be < dx^2 / alpha.
    
    # Stability check
    dt_stab = dx**2 / alpha_code
    print(f"Stability dt < {dt_stab:.2e}")
    if dt > dt_stab:
        dt = dt_stab * 0.5
        print(f"Adjusted dt to {dt:.2e}")
        
    kinetic_op_half = np.exp(-1j * (alpha_code/2) * k2 * (dt/2))
    kinetic_op_full = np.exp(-1j * (alpha_code/2) * k2 * dt)
    
    os.makedirs('output', exist_ok=True)
    
    start_time = time.time()
    
    for step in range(n_steps + 1):
        # 1. Potential Step (Kick)
        rho = np.abs(psi)**2
        rho_mean = np.mean(rho)
        delta = rho - rho_mean
        
        # Poisson Solve: -k^2 Phi = 4 pi G delta
        # Phi_k = 4 pi G delta_k / k^2 (minus signs cancel? Del^2 = -k^2)
        # Del^2 Phi = 4 pi G rho.  -k^2 Phi_k = 4 pi G rho_k. Phi_k = -4 pi G rho_k / k^2.
        delta_k = np.fft.fftn(delta)
        phi_k = -4 * np.pi * G_code * delta_k / k2
        phi_k[0,0,0] = 0
        phi = np.real(np.fft.ifftn(phi_k))
        
        # Apply Potential: psi = exp(-i V/alpha dt) psi ??
        # Eq: i dpsi/dt = V psi -> psi(t) = exp(-i V dt) psi (if alpha=1 in V term?)
        # Full eq: i dpsi/dt = -alpha/2 del^2 psi + (1/alpha) Phi psi ???
        # Standard: i hbar dpsi/dt = H psi.  / hbar -> i dpsi/dt = (-hbar/2m del^2 + m/hbar Phi) psi
        # coeff 1 = hbar/2m = alpha/2.
        # coeff 2 = m/hbar = 1/alpha.
        # So Potential op is exp(-i * (1/alpha) * Phi * dt)
        
        # Kick (Full step if using Drift-Kick-Drift or similar? No, split step usually D(dt/2) K(dt) D(dt/2))
        # First step structure
        if step == 0:
             psi = np.fft.ifftn(kinetic_op_half * np.fft.fftn(psi))
        
        # Kick
        psi *= np.exp(-1j * (1.0/alpha_code) * phi * dt)
        
        # Drift
        psi = np.fft.ifftn(kinetic_op_full * np.fft.fftn(psi)) # This effectively does D(dt/2) for next step too
        
        # IO
        if step % dump_interval == 0:
            print(f"Step {step}/{n_steps} | Wall: {time.time()-start_time:.2f}s | Max Dens: {np.max(rho):.2f}")
            with h5py.File(f'output/snap_{step:04d}.h5', 'w') as f_out:
                f_out.create_dataset('density', data=rho)
                f_out.create_dataset('phase', data=np.angle(psi))

if __name__ == "__main__":
    run_simulation()
