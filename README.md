# Hofstadter Butterfly

Manuscript archive for the twisted-graphene project and a separate conceptual perspective on Bellissard's gap labeling.

The working manuscript is also available in [Overleaf](https://www.overleaf.com/project/6aaa40e72ddf45865b2e3bce). Updates between the two services are manual; this repository is not automatically synchronized with Overleaf.

The PRL-format Letter and its supplement are by Vladislav Efremkin, Thomas D. Kühne, and Emil Prodan. The separate perspective is by Thomas D. Kühne, Vladislav Efremkin, and Emil Prodan. These are draft manuscripts. The source DFT image data have not been altered.

## Contents

- [Manuscript PDF](output/pdf/twisted_graphene_letter.pdf) and [Supplemental Material PDF](output/pdf/twisted_graphene_supplement.pdf).
- `main.tex`: Letter in REVTeX 4.2.
- `supplement.tex`: input audit, complete DFT survey, model normalization and reproducibility, and necessary controls.
- `bellissard_perspective.tex`: independent conceptual article, *From Bellissard's gap labeling to atomistic quasiperiodic bilayers*.
- [Perspective PDF](output/pdf/bellissard_perspective.pdf), with a finite-model spectral-moment illustration and no first-principles results.
- `bellissard_references.bib`: additional references for the perspective.
- `references.bib`: DOI-checked primary literature.
- `figures/`: PDF figures used by the two documents.
- `data/original_panels/`: the 13 original DFT image panels.
- `data/structures_d2p15/`: 51 supplied coordinate/input configurations.
- `data/tb_reproduction/`: independently calculated dimensionless eigenvalues.
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
```

Existing model-result files are preserved rather than silently recalculated. Delete or move a specific result only if deliberately requesting a fresh run. Default calculations are five sets of 41 dense 3020-dimensional diagonalizations and may take appreciable time on a different machine.

Verified development environment: Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, and Pillow 12.3.0; TeX Live 2025. The original CP2K version requires confirmation from run headers. The supplied job script names CP2K/2023.1.

## Data provenance and scientific scope

The principal data package was circulated by Vladislav Efremkin on July 23, 2026. It supersedes earlier subsets for this draft. The 51 detailed configurations are for d = 2.15 angstrom; input settings at the other distances and for DZVP are not independently established by the available files. The DOS maps cannot be reconstructed quantitatively without the original output and plotting workflow.

The tight-binding code implements the supplied all-pairs exponential-hopping model, divided by its nearest-neighbor intralayer hopping. It is not a fitted graphene Hamiltonian. The small raw scale exp(-1/0.03) is handled by an exact energy and DOS change of variables. No electronvolt conversion, topological index, or bulk convergence is inferred from it.

The perspective uses the same finite-model eigenvalue archive for a distinct spectral-moment analysis, tabulated in `data/bellissard_model_moments.csv`. Its figure does not reproduce the Letter's DOS maps. The shared model and data provenance should be disclosed when submitting the related manuscripts. The perspective does not include the DFT survey or establish a bulk gap or topological invariant.

The checks in `data/model_checks.json` validate exact energy rescaling, angle-reflection symmetry, the DOS Jacobian, sorted finite eigenvalues, and zero matrix trace. They do not validate the DFT data or establish a topological gap.
