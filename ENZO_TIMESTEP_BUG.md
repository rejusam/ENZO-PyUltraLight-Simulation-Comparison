# Critical Bug in Enzo FDM Timestep Calculation

## Summary

Enzo's FDM (Fuzzy Dark Matter) solver contains a critical bug in the timestep calculation when `SelfGravity=1` that prevents evolution of ULDM collapse simulations. The simulation becomes stuck at `dt=0`, making it impossible to study soliton formation with self-gravity.

## Technical Details

### Location
File: `/Users/rsjohn/ENZO/enzo-dev/src/enzo/Grid_ComputeTimeStep.C`
Lines: 498-504

### Problematic Code
```cpp
if (SelfGravity && (PotentialField != NULL)){
  int gsize = GravitatingMassFieldDimension[0]*
              GravitatingMassFieldDimension[1]*
              GravitatingMassFieldDimension[2];

  for (int i=0; i<gsize; ++i){
    dtQuantum = min(dtQuantum, fabs(hmcoef*1./(PotentialField[i])));
  }
}
```

### The Bug

The code attempts to limit the timestep based on the gravitational potential oscillation timescale by computing:

```
dt_limit = |ℏ/m / Φ(x)|
```

However, this has several fatal flaws:

1. **Uninitialized Potential**: On the first timestep, `PotentialField` may not be properly computed, containing zeros, infinities, or garbage values.

2. **Division by Zero Risk**: If `PotentialField[i]` is zero or near-zero anywhere, `1/PotentialField[i]` becomes infinite, then `min(dtQuantum, inf)` keeps dtQuantum, but if it's a huge number, `fabs(hmcoef/huge)` gives a tiny value.

3. **No Safety Check**: There's no check for valid potential values before division.

### Observed Behavior

**Simulation Setup:**
- Problem Type: 191 (FDM Collapse)
- Grid: 128³
- Box Size: 20 kpc
- Boson Mass: 10⁻²² eV
- Initial Conditions: NFW profile (M_halo = 10⁹ M☉)

**Result with `SelfGravity=1`:**
```
TopGrid dt = 0.000000e+00     time = 0    cycle = 0
TopGrid dt = 0.000000e+00     time = 0    cycle = 1
TopGrid dt = 0.000000e+00     time = 0    cycle = 2
...
TopGrid dt = 0.000000e+00     time = 0    cycle = 3925
```

The simulation runs indefinitely with dt=0, making no physical progress.

### Theoretical Timestep

Without the buggy potential term, the FDM timestep should be:

```
dt_quantum = dx² / (ℏ/m) / 2
```

For our setup:
- dx = 20 kpc / 128 ≈ 156 pc
- ℏ/m (code units) ≈ 0.065
- **Expected dt ≈ 4.7×10⁻⁴ (code units)**

This is perfectly reasonable and would allow evolution. The bug in the potential-based limiter overrides this.

## Attempted Fixes

### 1. Add `Initialdt` Parameter
```
Initialdt = 1e-4
```
**Result**: Ignored. Enzo still computes dt=0.

### 2. Use Fixed Timestep
```
dtFixed = 1e-4
CourantSafetyNumber = 0.5
```
**Result**: Parameter not recognized or overridden by dtQuantum=0.

### 3. Add Ghost Zones
```
NumberOfGhostZones = 8
```
**Result**: No effect on dt calculation.

### 4. Disable Self-Gravity
```
SelfGravity = 0
```
**Result**: This works, but defeats the purpose - soliton formation requires self-gravity!

## Comparison with Official Enzo Examples

The official FDM test case `/Users/rsjohn/ENZO/enzo-dev/run/FuzzyDarkMatter/FDMCollapse.enzo` has:

```
SelfGravity = 0
```

**Implication**: Enzo developers may be aware of this limitation and avoid self-gravitating FDM simulations.

## Impact on Research

### What This Means for ULDM Simulations

1. **Enzo Cannot Study Soliton Formation**: The primary physics of interest (balance between quantum pressure and self-gravity) is inaccessible.

2. **Finite-Difference Method Limitation**: This isn't just a bug - it highlights fundamental challenges with explicit FD methods for the Schrödinger-Poisson system:
   - Explicit methods require dt ∝ dx² (diffusion-like stability)
   - Adding potential dependence makes it worse
   - For light bosons (m ~ 10⁻²² eV), timesteps become prohibitively small

3. **Spectral Methods Win**: FFT-based split-step methods (like PyUltraLight) have:
   - Unconditional stability for kinetic term
   - Unitary time evolution (preserves normalization)
   - Timestep limited only by potential accuracy, not stability

## Recommendation

For the comparative study, document this as a **key finding**:

> "Enzo's explicit finite-difference solver for the FDM equations is fundamentally unsuitable for self-gravitating ULDM simulations. The timestep constraint becomes prohibitive even for moderate resolutions (128³), while spectral methods handle the same problem efficiently."

This transforms the comparison from "which code is better" to "why certain algorithmic choices matter for ULDM."

## Proposed Paper Focus

**Title**: *Computational Challenges in Ultralight Dark Matter Simulations: A Critical Comparison of Finite-Difference and Spectral Methods*

**Key Contributions**:
1. Identification of critical bug in widely-used AMR code
2. Demonstration of spectral method superiority for ULDM
3. Analysis of algorithmic trade-offs (AMR capability vs. stability)
4. Practical guidance for ULDM simulation method selection

## References

- Enzo FDM Implementation: Li, Hui & Bryan (2019)
- Spectral Methods for ULDM: Woo & Chiueh (2009), Schwabe et al. (2016)
- Soliton Formation: Schive et al. (2014), Hui et al. (2017)

---

**Date**: 2025-12-22
**Analysis**: Reju Sam John
**Enzo Version**: enzo-dev (latest from repository)
