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

Verified development environment: Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, and Pillow 12.3.0; TeX Live 2025. The original CP2K version requires confirmation from run headers. The supplied job script names CP2K/2023.1.

## Data provenance and scientific scope

The principal data package was circulated by Vladislav Efremkin on July 23, 2026. It supersedes earlier subsets for this draft. The 51 detailed configurations are for d = 2.15 angstrom; input settings at the other distances and for DZVP are not independently established by the available files. The DOS maps cannot be reconstructed quantitatively without the original output and plotting workflow.

The tight-binding code implements the supplied all-pairs exponential-hopping model, divided by its nearest-neighbor intralayer hopping. It is not a fitted graphene Hamiltonian. The small raw scale exp(-1/0.03) is handled by an exact energy and DOS change of variables. No electronvolt conversion, topological index, or bulk convergence is inferred from it.

The perspective uses the same finite-model eigenvalue archive for a distinct spectral-moment analysis, tabulated in `data/bellissard_model_moments.csv`, and refines two distance series for a separate butterfly/state-count illustration. These new maps use finer angular sampling and a common broadening, but share the underlying Hamiltonian and part of the spectral data with the Letter. The shared model and data provenance should be disclosed when submitting the related manuscripts. The perspective does not include the DFT survey or establish a bulk gap or topological invariant.

The checks in `data/model_checks.json` validate exact energy rescaling, angle-reflection symmetry, the DOS Jacobian, sorted finite eigenvalues, and zero matrix trace. They do not validate the DFT data or establish a topological gap.
