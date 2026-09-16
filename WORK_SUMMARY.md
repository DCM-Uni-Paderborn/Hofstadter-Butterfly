# Work completed

## Manuscript and Supplemental Material

- Prepared a four-page PRL-format manuscript and a nine-page Supplemental Material document from the supplied Efremkin/Prodan correspondence, figures, structures, CP2K inputs, and model notebooks.
- Integrated the discussion from the earlier analysis of Emil's energy scale.
- Assembled the full supplied DFT image survey, including the SZV/DZVP comparison. Original DFT image data were retained without recoloring or interpolation.
- Added geometry illustrations, an angle-dependent atom-count analysis, and independently regenerated tight-binding DOS figures.
- Checked 19 primary-literature references and distinguished the proposed contribution from prior acoustic, mechanical-bilayer, and compressed-graphene work.

## Methods and numerical checks

- Audited all 51 supplied configurations at an interlayer distance of 2.15 angstrom, including the changing upper-layer carbon and hydrogen counts, cell geometry, and explicitly specified CP2K settings.
- Checked the distinction between operators and their matrix representations, the generalized Kohn-Sham eigenproblem, and the Lowdin-type projection used in the CP2K 2023.1 PDOS implementation.
- Derived a consistent dimensionless energy and DOS normalization for the supplied exponential-hopping model. The small raw hopping scale is not, by itself, evidence of numerical failure; the rescaling does not establish a physical conversion to DFT energies in electronvolts.
- Independently diagonalized five model parameter sets, each at 41 nonnegative angles with 3020 eigenvalues per angle. The negative-angle spectra follow from an explicitly justified layer-exchange symmetry.
- Tested energy rescaling, angle-reflection symmetry, the DOS Jacobian, eigenvalue ordering and finiteness, and the zero-trace sum rule. Detailed numerical results are in `data/model_checks.json`.

## Scientific scope and remaining work

- Kept the observed spectral reconstruction separate from unproven claims of bulk topological gaps, Chern numbers, fractality, or quantized pumping.
- Documented the missing DFT raw output and postprocessing information, and the need for matched monolayer references, convergence and boundary checks, compressed-geometry stability checks, and an appropriate topological analysis.
- Collected the required follow-up decisions in `AUTHOR_REVIEW.md`. No new DFT calculations or topological-invariant calculations have been performed.

## Author details, validation, and availability

- Set the author order to Vladislav Efremkin, Thomas D. Kuehne, and Emil Prodan. Thomas D. Kuehne carries the corresponding-author star with `tkuehne@cp2k.org` in both documents.
- Shortened the AI-assisted-tools acknowledgment as requested.
- Compiled and visually checked the PDFs; the build has no undefined references or overfull boxes.
- Uploaded the manuscript, supplement, data, analysis code, and PDFs to the existing Overleaf project.
- Created `DCM-Uni-Paderborn/Hofstadter-Butterfly` as a version-controlled archive and subsequently made it public at the author's request. Private email exports, credentials, and cluster submission scripts are excluded.

Overleaf and GitHub are not automatically synchronized. No journal submission or email to the coauthors has been sent.
