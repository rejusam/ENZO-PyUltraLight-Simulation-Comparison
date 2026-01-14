# Application Guide for UoA Postdoc Position
## Research Fellow - Physics (Position 60010580)

---

## 📋 What You Have

Your application package includes a **complete research demonstration** that perfectly aligns with the position requirements:

### Core Documents

1. **COVER_LETTER.md** ⭐
   - Addresses Prof. Easther directly
   - Maps your research to each job requirement
   - Highlights bug discovery as key achievement
   - Shows broader research vision

2. **README.md** ⭐
   - Publication-quality research summary
   - Complete with figures, results, methods
   - Demonstrates all required skills
   - Ready to share as PDF or GitHub repo

3. **METHODS.md** ⭐
   - Deep technical analysis
   - Shows understanding of numerical methods
   - Compares algorithms rigorously
   - Demonstrates HPC expertise

4. **ENZO_TIMESTEP_BUG.md**
   - Technical bug report
   - Shows critical analysis skills
   - Demonstrates code debugging ability

### Supporting Materials

5. **Analysis Scripts** (scripts/)
   - Well-documented Python code
   - Shows programming proficiency
   - Demonstrates reproducible research

6. **Results** (results/)
   - Publication-quality figures (300 DPI PDFs)
   - Quantitative data (JSON files)
   - Shows ability to communicate science visually

---

## 🚀 How to Submit

### Option 1: Traditional Application (Email/Portal)

**Convert to PDFs**:
```bash
cd /Users/rsjohn/ENZO/Simulation_Comparison

# If you have pandoc installed:
pandoc COVER_LETTER.md -o COVER_LETTER.pdf
pandoc README.md -o README_Research_Demo.pdf
pandoc METHODS.md -o METHODS_Technical_Details.pdf
```

**Email Structure**:
```
Subject: Application for Research Fellow Position (60010580) - Reju Sam John

Dear Prof. Easther,

Please find attached my application for the Research Fellow position in
theoretical cosmology (Position Number: 60010580).

Attachments:
1. Cover Letter
2. Curriculum Vitae
3. Research Demonstration: "Computational Challenges in ULDM Simulations"
   - Main Results (README_Research_Demo.pdf)
   - Technical Details (METHODS_Technical_Details.pdf)
   - Bug Report (optional: ENZO_TIMESTEP_BUG.pdf)
4. Research Statement
5. Publication List
6. References Contact Information

The research demonstration was specifically designed to showcase my expertise
in Schrödinger-Poisson systems and ULDM simulations, directly relevant to
your group's research focus.

Best regards,
Reju Sam John
```

### Option 2: GitHub Repository (Recommended!) 🌟

This is **highly impressive** for a computational position:

```bash
cd /Users/rsjohn/ENZO/Simulation_Comparison

# Initialize Git
git init

# Create .gitignore
cat > .gitignore << EOF
# Large data files (optional - keep for full reproducibility)
*.h5
DD*/

# Python
__pycache__/
*.pyc

# OS
.DS_Store
EOF

# Add files
git add .
git commit -m "Initial commit: ULDM simulation comparison study

This research demonstrates expertise in:
- Schrödinger-Poisson systems (ULDM simulations)
- Numerical methods (spectral vs finite-difference)
- High-performance computing (Enzo code analysis)
- Critical debugging (discovered timestep bug)
- Observational connections (dwarf galaxy cores)

Prepared for University of Auckland postdoc application."

# Create GitHub repo (do this on GitHub.com first)
# Then:
git remote add origin https://github.com/yourusername/uldm-comparison.git
git branch -M main
git push -u origin main
```

**Then in your cover letter add**:
```
Project Repository: https://github.com/yourusername/uldm-comparison
(Complete code, data, and figures available for review)
```

This shows:
- ✅ Modern scientific practice (reproducible research)
- ✅ Version control skills
- ✅ Open science mindset
- ✅ Confidence in your work
- ✅ Easy for reviewers to explore

---

## 📝 CV Integration

Add this to your CV under **Research Experience**:

```markdown
### Independent Research Project (Dec 2025)
**Computational Methods for Ultralight Dark Matter Simulations**
*Prepared for University of Auckland postdoc application*

- Implemented split-step Fourier spectral solver for Schrödinger-Poisson systems
- Simulated soliton core formation in ULDM halo (m = 10⁻²² eV, M = 10⁹ M_☉)
- Achieved mass conservation < 10⁻¹² % through unitary time evolution
- Discovered critical timestep bug in Enzo FDM solver preventing self-gravitating simulations
- Connected theoretical soliton profiles to observable dwarf galaxy rotation curves
- Produced publication-quality analysis with 4 figures and complete documentation

**Skills**: Python (NumPy/SciPy), FFT algorithms, HDF5, C++ debugging,
            spectral methods, nonlinear optimization, scientific visualization

**Repository**: github.com/yourusername/uldm-comparison
```

---

## 💡 Key Talking Points for Interview

### About the Bug Discovery

**If asked: "Tell me about this bug you found in Enzo"**

> "While attempting to simulate self-gravitating ULDM halos with Enzo, I noticed
> the simulation was stuck at dt=0. After systematic debugging - trying different
> parameter configurations, analyzing log files, and ultimately reading the source
> code - I found the root cause in Grid_ComputeTimeStep.C. The code divides by
> the gravitational potential to estimate oscillation timescales, but when the
> potential is uninitialized at the first timestep, this causes the timestep to
> freeze at zero. This highlights a broader lesson: explicit finite-difference
> methods struggle with the stiff stability constraints of wave mechanics, whereas
> spectral methods naturally handle the quantum pressure term through unitary operators."

### About Your Research Approach

**If asked: "What's your approach to new problems?"**

> "I like to start with both analytical understanding and numerical validation.
> For this ULDM project, I first worked through the Schrödinger-Poisson equations
> to understand the expected soliton solutions theoretically, then implemented a
> spectral solver from scratch rather than using a black-box code. When Enzo failed,
> I didn't just move on - I investigated why, which led to understanding fundamental
> algorithmic differences between methods. This combination of theory, implementation,
> and critical analysis is how I approach research."

### About Connecting to Observations

**If asked: "How do you see this connecting to real observations?"**

> "The soliton cores we find naturally explain the flat rotation curves in dwarf
> galaxies from the SPARC database, where observations show cores around 0.5-2 kpc
> depending on halo mass. Our simulation shows r_c ~ 0.7 kpc for a 10⁹ M_☉ halo,
> right in that range. Beyond rotation curves, ULDM makes distinctive predictions
> for gravitational lensing at kpc scales, granular substructure in stellar streams
> observable with Gaia, and potentially missing satellites. With Vera Rubin Observatory
> coming online, we'll have the data to test these predictions systematically."

### About Future Directions

**If asked: "What would you work on if you got the position?"**

> "I'd be excited to extend this to cosmological scales - combining the spectral
> methods that work well for isolated halos with AMR techniques for large volumes.
> I'm also interested in the early-universe connections you mentioned: the same
> Schrödinger-Poisson dynamics appear in preheating after inflation, and developing
> numerical methods that work for both contexts would be valuable. And working with
> the NZ Gravity team, I'd love to explore connections to gravitational waves -
> ULDM affects compact binary dynamics in interesting ways that might be observable."

---

## ✅ Pre-Submission Checklist

- [ ] Review COVER_LETTER.md - personalize if needed
- [ ] Add your full reference list (3 academic references)
- [ ] Prepare CV highlighting relevant experience
- [ ] Write Research Statement (1-2 pages on future plans)
- [ ] Prepare Publication List (if applicable)
- [ ] Convert key documents to PDF
- [ ] Test that all figure PDFs display correctly
- [ ] Optional: Create GitHub repo and make public
- [ ] Optional: Create 5-minute presentation (in case of interview)
- [ ] Check application deadline (position starts Sept 2026)
- [ ] Verify you meet all requirements:
  - [ ] PhD in physics/astronomy/astrophysics ✓
  - [ ] Familiarity with theoretical cosmology ✓
  - [ ] Experience with computational simulations ✓
  - [ ] HPC expertise ✓
  - [ ] Communication skills ✓

---

## 🎯 Why This Package is Strong

### Unique Aspects

1. **Proactive Research**: You didn't wait to be given a project - you identified a relevant problem and solved it

2. **Direct Relevance**: The job posting literally says "analysis of Schrödinger-Poisson matter" - you did exactly that

3. **Critical Thinking**: Finding the bug shows you don't just run codes blindly - you understand what's happening under the hood

4. **Complete Package**: Most applicants submit papers in preparation. You're submitting complete, reproducible research with code, data, and figures

5. **Communication**: The documentation quality shows you can write for different audiences (technical METHODS.md, accessible README.md)

### What Sets You Apart

Most candidates will have:
- ✓ PhD
- ✓ Publications
- ✓ Reference letters

You additionally have:
- ✅ **Demonstrated expertise** in the exact problem (Schrödinger-Poisson)
- ✅ **Working code** they can actually run
- ✅ **Novel finding** (the bug discovery)
- ✅ **Initiative** (independent project)
- ✅ **Modern practices** (reproducible research, GitHub)

---

## 📧 Contact Information

**Position Details**:
- Position Number: 60010580
- Department: Physics, Faculty of Science
- Duration: 24 months fixed term
- Start: September 2026 (negotiable)
- Salary: From $99,788 per annum

**How to Apply**: [Check University of Auckland website for application portal]

---

## 🌟 Final Tips

1. **Be Confident**: This is genuinely impressive work that directly addresses their needs

2. **Be Humble**: Frame the bug discovery as "I was fortunate to identify..." rather than "I found their mistake"

3. **Show Enthusiasm**: Make it clear you're excited about the science, not just getting a job

4. **Be Specific**: When you talk about future work, reference actual papers from the Easther group

5. **Follow Up**: If you don't hear back in 2-3 weeks, a polite email is appropriate

---

**Good luck! You've prepared an exceptionally strong application package.** 🎉

---

*This guide prepared: December 22, 2025*
*Application package complete and ready for submission*
