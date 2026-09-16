# Validation record, 16 September 2026

- Both documents compile with TeX Live 2025 / REVTeX 4.2.
- Letter: four pages, three figures, 19 primary-literature references plus the Supplemental Material entry.
- Supplemental Material: nine pages, five figures, one settings table, 17 numbered equations.
- No undefined citations/references or overfull boxes in the final LaTeX logs.
- Every PDF page was rendered and visually reviewed; figure labels, equations, table boundaries, and bibliography were checked. Model contrast was improved using an explicitly documented common color scale; the original DFT images were not changed.
- The 51 geometry files were parsed to count lower-layer carbon, upper-layer carbon, and hydrogen.
- Five model parameter sets were diagonalized independently, each at 41 nonnegative angles with 3020 eigenvalues per angle. Negative angles are obtained by an exactly justified symmetry.
- Independent small-matrix tests verify the energy rescaling and angle-reflection spectra to below 5e-15 absolute error. The DOS Jacobian test is below 7e-16. All saved eigenvalues are finite and ordered, and the maximum absolute trace residual is below 1e-12.
- Primary bibliographic metadata were checked against the DOI records. The DFT projection formula was checked against the CP2K 2023.1 PDOS source, which forms S^(1/2) C.

These checks establish technical consistency of the draft and model reproduction, not DFT convergence, experimental feasibility, or a topological invariant. The missing scientific evidence is listed in AUTHOR_REVIEW.md. No new DFT calculation or journal submission was made.
