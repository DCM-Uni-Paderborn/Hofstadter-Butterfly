# Hofstadter Butterfly

Manuscript archive for the twisted-graphene project and a separate conceptual perspective on Bellissard's gap labeling.

The working manuscript is also available in [Overleaf](https://www.overleaf.com/project/6aaa40e72ddf45865b2e3bce). Updates between the two services are manual; this repository is not automatically synchronized with Overleaf.

The PRL-format Letter and its supplement are by Vladislav Efremkin, Thomas D. Kühne, and Emil Prodan. The separate perspective is by Thomas D. Kühne, Vladislav Efremkin, and Emil Prodan. These are draft manuscripts. The source DFT image data have not been altered.

## Contents

- [Manuscript PDF](output/pdf/twisted_graphene_letter.pdf) and [Supplemental Material PDF](output/pdf/twisted_graphene_supplement.pdf).
- `main.tex`: Letter in REVTeX 4.2.
- `supplement.tex`: input audit, complete DFT survey, model normalization and reproducibility, and necessary controls.
- `bellissard_perspective.tex`: independent conceptual article, *From Bellissard's gap labeling to atomistic quasiperiodic bilayers*.
- [Perspective PDF](output/pdf/bellissard_perspective.pdf), with finite-model spectral moments, refined butterfly spectra, integrated state counts, and compact resolution controls.
- [Butterfly and state-count figure](figures/bellissard_butterfly_ids.pdf): two refined angle scans and their DOS/counting-function cross sections, including an uncoupled-layer reference.
- `bellissard_references.bib`: additional references for the perspective.
- `references.bib`: DOI-checked primary literature.
- `figures/`: PDF figures used by the two documents.
- `data/original_panels/`: the 13 original DFT image panels.
- `data/structures_d2p15/`: 51 supplied coordinate/input configurations.
- `data/dft_audit/`: September 2026 CP2K inventory, checksummed spin-resolved spectral extracts, monolayer DOS, and the SI coverage table.
- `data/tb_reproduction/`: independently calculated dimensionless eigenvalues.
- `data/bellissard_butterfly/`: 241-angle model spectra, DOS, exact finite state counts, and numerical checks for the perspective.
- `data/bellissard_robustness/`: three-radius, three-width controls, uncoupled-layer spectra, and per-case numerical comparisons.
- `analysis/`: model reproduction, figure generation, and consistency checks.

Private correspondence, account information, mail exports, and internal working or validation notes are not part of the manuscript package. Cluster submission scripts are excluded from the distributable archive; they are not required to compile the paper.

## Build

Requirements: TeX Live with REVTeX 4.2 and latexmk. From this directory:

```sh
latexmk -pdf -outdir=build main.tex
latexmk -pdf -outdir=build supplement.tex
latexmk -pdf -outdir=build bellissard_perspective.tex
```

The perspective is a separate root document in the same Overleaf project. Select `bellissard_perspective.tex` for its compilation without replacing `main.tex`. It currently uses the standard LaTeX article class and can be transferred to the proceedings template when supplied.

The figures are already included. To reproduce the model calculations, use Python with NumPy, SciPy, Matplotlib, and Pillow:

```sh
python analysis/reproduce_tb.py
python analysis/reproduce_tb.py --distances 0.97 --xi 0.1
python analysis/reproduce_tb.py --distances 0.97 --xi 0.3
python analysis/check_model.py
python analysis/make_figures.py
python analysis/bellissard_model_moments.py
python analysis/bellissard_butterfly.py --calculate-only
python analysis/bellissard_robustness.py
python analysis/bellissard_butterfly.py --plot-only
```

Existing model-result files are preserved rather than silently recalculated. Delete or move a specific result only if deliberately requesting a fresh run. Default calculations are five sets of 41 dense 3020-dimensional diagonalizations and may take appreciable time on a different machine.

The separate butterfly script refines two of those distance series to 241 angles, retaining the 41 archived spectra in each series and calculating 200 additional angles per distance. It uses two worker processes with three numerical-library threads each. Checkpoints are stored under `tmp/`, outside the archive. Both DOS maps use the same Gaussian width of 0.02 in nearest-neighbor hopping units. The normalized finite state count is calculated without broadening. The figure additionally requires the uncoupled reference from the robustness script. All required data are included in this archive. For a fresh calculation without archived control data, run `python analysis/bellissard_butterfly.py --calculate-only`, then `python analysis/bellissard_robustness.py`, then `python analysis/bellissard_butterfly.py --plot-only`.

The robustness script tests R/a = 20, 25, 30 (979, 1510, 2167 sites per layer), angles 0, 15, 30 degrees, both distances 0.99 and 0.97, and Gaussian widths 0.01, 0.02, 0.04. It retains the six matching R/a = 25 spectra, calculates twelve additional bilayer spectra and three monolayer spectra, and archives normalized DOS arrays on an energy grid of spacing 0.0025. The uncoupled bilayer has two copies of the monolayer spectrum, hence the same per-state DOS. The exact integral of the absolute difference of state counts is independently checked against the one-dimensional Wasserstein distance. DOS normalization, Gaussian second moments, and the interlayer trace identity are also checked. These selected-angle controls assess finite-size sensitivity, not a bulk extrapolation or topological invariant.

Verified development environment: Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, and Pillow 12.3.0; TeX Live 2025. The survey output headers identify CP2K 2023.1 except for the supplied d = 3.30 angstrom, theta = 49.8 degree record (2026.1 development version). The ancillary 3.35 angstrom series has separate versions and settings and is not part of the ten-distance image survey.

## Data provenance and scientific scope

The July 23, 2026 package supplies the original image panels, 51 detailed input configurations at d = 2.15 angstrom, and the model notebooks. The September 16, 2026 package additionally supplies CP2K inputs, outputs, PDOS files, plotting and structure-generation notebooks, and seven monolayer DOS files. The original bilayer images remain unchanged. The state-resolved SZV output for d = 2.15 angstrom is still absent from the additional package. Some other PDOS files are incomplete, including truncated files associated with converged outputs. The CSV inventory and SI table distinguish these cases from audited columns.

The bilayer plotting notebooks read **alpha-spin p_z weights on lower-layer atoms 1--1250**, not the full layer PDOS or the sum of both spins. Their Hartree-to-eV factor is 27.211384. They select eigenvalues in [-7, 3] eV relative to each file's header Fermi energy and sum `exp(-(E-e)^2/eta^2)/(2*pi*eta)` with eta = 0.02 eV on 6000 points. This corresponds to Gaussian standard deviation 0.01414 eV and FWHM 0.03330 eV. Each angular column is divided by its own maximum. No occupation factor is applied. The kernel's non-unit integral cancels under this maximum normalization. The highest printed unoccupied state can be below +3 eV, especially for DZVP, so dark regions above that limit are not physical gap evidence.

### Reproduce the DFT postprocessing

The source directory and image named `2-85` actually have a planar separation of **2.84 angstrom** in all 51 coordinate files. The SI label and coverage table use the coordinate-derived value; source filenames are retained for traceability. No geometry or calculation was changed.

No new electronic-structure or model calculation is needed for these commands:

```sh
python analysis/plot_dft_references.py
python analysis/replot_pdos.py --series SZV/IL_25_z2-20 --output tmp/szv_d2p20_replot
```

The first command reproduces the SI monolayer figure from the archived binned total DOS, adding both spins and using a common Gaussian standard deviation of 0.10 eV. Only the header Fermi energy is subtracted, without the empirical size-dependent shifts in the original notebook. Each curve is normalized to its maximum within [-7, 3] eV. The source size labels are retained because the monolayer input geometries are not supplied.

The second command reconstructs the supplied bilayer display convention from the spectral extracts, with explicit masks for incomplete columns and energies above the last printed state. It writes a diagnostic PDF and numerical arrays without replacing any original manuscript figure. Gaussian-tail truncation at eight kernel widths is checked against direct summation at 31 energies per column. The `coverage` array uses 0 for unavailable columns, 1 for the printed-state region, and 2 above the last printed state. Colors remain separately maximum-normalized for each angle, not an absolute DOS.

`data/dft_audit/spectra/*.npz` contains arrays keyed by configuration index and spin. Their four columns are the printed eigenvalue in hartree, occupation, lower-layer p_z weight, and summed lower-layer orbital weight. Separate scalar keys store the Fermi energy in hartree and angle in degrees. The CSV provides row counts, output-orbital counts, convergence and termination flags, energy limits, and source SHA-256 hashes. Records with incomplete PDOS are retained and flagged, not silently treated as full spectra. The monolayer archive retains all six numerical columns of the original `.dos` files and both header Fermi energies.

To regenerate the inventory and extracts from the complete author-supplied data folder, run:

```sh
python analysis/audit_dft_data.py --source /path/to/data_for_Thomas
```

This reads raw data only and copies no private correspondence, account-specific logs, cluster scripts, or access links into the public archive. The full raw archive is retained by the authors. A monolayer reference on its own Fermi-energy axis is not a vacuum-aligned, matched noninteracting-layer subtraction.

The tight-binding code implements the supplied all-pairs exponential-hopping model, divided by its nearest-neighbor intralayer hopping. It is not a fitted graphene Hamiltonian. The small raw scale exp(-1/0.03) is handled by an exact energy and DOS change of variables. No electronvolt conversion, topological index, or bulk convergence is inferred from it.

The perspective uses the same finite-model eigenvalue archive for a distinct spectral-moment analysis, tabulated in `data/bellissard_model_moments.csv`, and refines two distance series for a separate butterfly/state-count illustration. These new maps use finer angular sampling and a common broadening, but share the underlying Hamiltonian and part of the spectral data with the Letter. The shared model and data provenance should be disclosed when submitting the related manuscripts. The perspective does not include the DFT survey or establish a bulk gap or topological invariant.

The checks in `data/model_checks.json` validate exact energy rescaling, angle-reflection symmetry, the DOS Jacobian, sorted finite eigenvalues, and zero matrix trace. They do not validate the DFT data or establish a topological gap.
