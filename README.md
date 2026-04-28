
# Tailored 3D Speckle Optimization (MATLAB)

MATLAB code for a two-step optimization pipeline:
1) **Speckle field optimization**: optimize a complex field whose propagated intensities match prescribed intensity statistics (PDFs) across multiple axial planes.
2) **SLM optimization**: realize a target complex field using a **phase-only SLM** with carrier grating and **first-order diffraction filtering**.

This repository accompanies:  
Codes for "Tailored Speckle Illumination Microscopy with Enhanced Sectioning and Image Quality"

Authors: SeungYun Han, KyeoReh Lee, Young Seo Kim, Chuan Li, Nicholas Bender, Kabish Wisal, Taeyun Ku, Jerome Mertz, Hui Cao

Corresponding author: Hui Cao (hui.cao@yale.edu)

For code details: SeungYun Han (seungyun.han@yale.edu)

---

## Repository structure

- `001_speckle_field_optimization/`  
  Multi-plane field optimization (NA-limited pupil constraint + axial propagation + intensity PDF transformation).
- `002_slm_optimization/`  
  Phase-only SLM pattern optimization with first-order filtering and Nesterov accelerated gradient.
- `example/`  
  Example of an optimized 3D speckle field.

Top-level:
- `main.m`: one-click entry script that runs the minimal example pipeline.

---

## Requirements

- MATLAB: **R20XXa/b** (tested on R2023b)
- Toolboxes:
  - 'Image Processing Toolbox'
  - Optional (GPU): Parallel Computing Toolbox

---

## Quick start

1) Clone the repository and open MATLAB in the repo root:
```matlab
cd('path/to/repo')
```

(Adjust the save path to your actual behavior.)

### 2) Clarify inputs: does `main.m` create/load `E0`?
You mentioned earlier `run_slm_optimization` needs `E0` in the workspace (scripts). In README you should state one of:

- **Option A:** `main.m` loads/creates `E0` automatically.  
- **Option B:** User must load `E0` before running `run_slm_optimization`.

If B, add a snippet:

```markdown
If you run `002_slm_optimization/run_slm_optimization.m` directly, you must define `E0` first (target complex field at the sample plane).
```

## Reproducibility notes
- Random seed is set in `main.m` via `rng(0)`.
- GPU is auto-detected; set `useGPU = 0` to force CPU.
- Scripts assume you run from the repo root (`cd` to the repo root before `main`).

## License
MIT License. See `LICENSE`.

## Citation
If you use this code, please cite:

SeungYun Han et al., *Tailored Speckle Illumination Microscopy with Enhanced Sectioning and Image Quality*, arXiv:2604.20112 (2026). [arXiv link](https://arxiv.org/abs/2604.20112)
