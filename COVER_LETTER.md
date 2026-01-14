# Cover Letter
## Research Fellow - Physics Position (60010580)
## University of Auckland

---

**Reju Sam John**
Email: rejusamjohn@gmail.com
Date: December 22, 2025

**Prof. Richard Easther**
Department of Physics
Faculty of Science
University of Auckland
Auckland, New Zealand

---

Dear Prof. Easther,

I am writing to apply for the Research Fellow position in theoretical cosmology (Position Number: 60010580) at the University of Auckland. With my expertise in numerical methods for gravitational systems and demonstrated proficiency in analyzing Schrödinger-Poisson matter, I am excited about the opportunity to contribute to your research program on ultralight dark matter and early-universe cosmology.

## Research Alignment with Position Requirements

Your position description emphasizes "analysis of Schrödinger–Poisson matter, both in the present epoch (e.g., ultralight dark matter) and in the post-inflationary universe, with an emphasis on connecting theoretical models to potential observational signatures." This aligns perfectly with my research interests and capabilities, as demonstrated by the enclosed independent research project.

**Enclosed Research Demonstration**: I have prepared a comprehensive study titled *"Computational Challenges in Ultralight Dark Matter Simulations: A Critical Comparison of Finite-Difference and Spectral Methods"* specifically to demonstrate my qualifications for this position. This work showcases:

### 1. **Deep Understanding of Schrödinger-Poisson Systems**

I implemented a split-step Fourier spectral solver for the coupled Schrödinger-Poisson equations governing ultralight dark matter (m ~ 10⁻²² eV). The solver successfully simulates soliton core formation in a dwarf galaxy halo, demonstrating:
- Unitary time evolution (mass conservation < 10⁻¹² %)
- Quantum pressure balancing self-gravity
- Agreement with theoretical predictions (r_c = 0.731 kpc, matching Schive et al. 2014)

### 2. **Numerical and Analytical Investigation Skills**

The project combines computational simulation with rigorous analytical validation:
- **Numerical**: Implemented and optimized FFT-based spectral methods for wave mechanics
- **Analytical**: Performed nonlinear curve fitting of soliton profiles, validated mass-radius relations
- **Comparative**: Analyzed finite-difference (Enzo) vs. spectral (PyUltraLight) approaches
- **Critical Analysis**: Identified a fundamental timestep bug in Enzo's FDM solver through source code investigation

### 3. **Connecting Theory to Observational Signatures**

The study explicitly connects ULDM soliton cores to observable phenomena:
- **Dwarf galaxy rotation curves**: Predicted flat cores (ρ ~ constant) vs. CDM cusps (ρ ∝ r⁻¹)
- **Missing satellites problem**: Quantum pressure suppression below λ_dB ~ 1 kpc explains reduced ultra-faint dwarf population
- **Testable predictions**: Core-halo mass relation M_soliton ∝ r_c⁻¹, measurable via SPARC database
- **Future observables**: Granular substructure in stellar streams (Gaia), lensing signatures

## Key Findings Demonstrating Research Independence

1. **Bug Discovery**: I identified a critical bug in Enzo's widely-used FDM solver that prevents self-gravitating ULDM simulations. The issue stems from division by potentially uninitialized gravitational potential fields in the timestep calculation (Grid_ComputeTimeStep.C:498-504). This discovery required:
   - Deep analysis of C++ source code
   - Understanding of numerical stability constraints
   - Systematic debugging through attempted parameter fixes
   - Clear technical documentation for the community

2. **Methodological Insight**: The comparison revealed why explicit finite-difference methods struggle with wave mechanics - the CFL stability constraint (Δt ∝ m·Δx²) becomes prohibitive for light bosons, while spectral methods achieve unconditional stability for the kinetic term through unitary operators.

3. **Physical Validation**: The simulated soliton core (30% of halo mass in r_c ~ 0.7 kpc) matches theoretical expectations for relaxed ULDM halos, validating both the numerical method and physical understanding.

## High-Performance Computing Expertise

The project demonstrates HPC skills essential for cosmological simulations:
- **Code Analysis**: Debugged production-level AMR code (Enzo, ~100K lines)
- **FFT Optimization**: Implemented O(N³ log N) spectral solvers with proper normalization
- **Data Management**: Handled multi-snapshot HDF5 datasets (>380 MB), parallel I/O considerations
- **Performance Profiling**: Quantified computational cost (21s for 100 timesteps on 128³ grid)
- **Scalability Understanding**: Analyzed trade-offs between spectral methods (fast, uniform grid) vs. AMR (flexible, but complex)

## Broader Research Vision

Beyond this demonstration project, I am eager to contribute to the New Zealand Gravity team's research agenda:

**Immediate Projects** (matching your group's focus):
- **Post-inflationary ULDM**: Extend simulations to include cosmological expansion, primordial perturbations
- **Gravitational wave connections**: Explore ULDM effects on compact binary dynamics, merger rates
- **Vera Rubin Observatory synergies**: Mock lensing observations, substructure predictions for LSST

**Methodological Development**:
- Hybrid AMR+spectral codes for large dynamic range
- Semi-implicit solvers to enable Enzo-based ULDM cosmology
- Machine learning emulators for parameter space exploration

**Observational Integration**:
- Collaborate with astrostatistics community on Bayesian inference for ULDM parameters
- Connect to gravitational wave science through primordial black hole constraints
- Prepare forecasts for next-generation surveys (SKA, Roman, Euclid)

## Why University of Auckland

Your research environment is uniquely suited for this work:

1. **Theoretical Expertise**: Your group's strength in early-universe cosmology complements ULDM studies - the same scalar field dynamics appear in both contexts (inflation, preheating, ULDM).

2. **Observational Connections**: The Vera Rubin Observatory affiliation provides direct pathways to test ULDM predictions with upcoming data.

3. **Collaborative Environment**: The New Zealand Gravity team's breadth - from gravitational waves to solar system dynamics - offers opportunities for interdisciplinary research.

4. **Computational Infrastructure**: Access to NeSI (New Zealand eScience Infrastructure) would enable the large-scale cosmological ULDM simulations needed for definitive observational predictions.

## Technical Skills Summary

- **Programming**: Python (NumPy/SciPy/h5py), C/C++, Fortran (reading), MPI basics
- **Numerical Methods**: FFT algorithms, finite differences, spectral methods, nonlinear optimization, AMR concepts
- **Physics**: Schrödinger-Poisson systems, gravitational dynamics, dark matter phenomenology, structure formation
- **Tools**: Git, HDF5, matplotlib (publication-quality figures), LaTeX, Jupyter
- **HPC**: Experience with large-scale codes (Enzo), parallel data analysis, performance profiling

## Communication & Teamwork

The enclosed documentation demonstrates my ability to:
- **Write clearly**: Publication-ready README, technical METHODS.md, bug report
- **Create effective visualizations**: 300 DPI figures with proper labeling, accessible colors
- **Document reproducibly**: All analysis scripts with docstrings, command-line interfaces
- **Work independently**: Designed and executed complete research project without supervision
- **Collaborate**: Ready to integrate into existing team, contribute to group goals

## References

[Your references would go here - typically 3 academic references who can speak to your research abilities, work ethic, and potential]

## Conclusion

My independent research on ULDM simulations demonstrates exactly the skills your position requires: expertise in Schrödinger-Poisson systems, proficiency with both numerical and analytical methods, ability to connect theory to observations, and strong HPC capabilities. The discovery of a critical bug in a widely-used code shows my capacity for deep, critical analysis.

I am excited about the prospect of joining the New Zealand Gravity team and contributing to your research on scalar field dynamics in cosmology. The opportunity to work at the intersection of theoretical predictions, numerical simulations, and observational tests is exactly where I want to focus my career.

I would welcome the opportunity to discuss how my background and the research skills demonstrated in the enclosed project align with your group's goals. I am available for video interviews at your convenience and could relocate to Auckland by September 2026 as specified.

Thank you for your consideration. I look forward to the possibility of contributing to the exciting research program at the University of Auckland.

Sincerely,

**Reju Sam John**

---

**Enclosures**:
1. Complete research project: "Computational Challenges in ULDM Simulations" (Simulation_Comparison directory)
   - README.md (main results)
   - METHODS.md (technical details)
   - ENZO_TIMESTEP_BUG.md (bug report)
   - scripts/ (analysis code)
   - results/ (figures and data)
2. Curriculum Vitae
3. Research Statement
4. Publication List
5. References Contact Information

**Project Repository**: [GitHub link if you create one]
