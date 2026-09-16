"""Audit the supplied CP2K data and archive publication-safe spectral extracts.

No electronic-structure or tight-binding calculation is performed. Inputs,
outputs, and PDOS files are read only. Account-specific logs and job scripts
are never copied to the manuscript archive.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HARTREE_EV = 27.211384  # Exact conversion used in the supplied plotting notebooks.


def match(text, pattern, default=""):
    found = re.search(pattern, text, re.MULTILINE)
    return found.group(1).strip() if found else default


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""


def read_pdos(path):
    with path.open() as stream:
        header = stream.readline()
        labels = stream.readline().split()
    orbital_labels = labels[labels.index("Occupation") + 1:]
    assert orbital_labels[:4] == ["s", "py", "pz", "px"], labels
    fermi = float(match(header, r"E\(Fermi\)\s*=\s*([-+\d.Ee]+)"))
    values = np.loadtxt(path, skiprows=2, ndmin=2)
    assert np.all(np.diff(values[:, 1]) >= 0), path
    assert np.all(values[:, 3:] >= -1e-7), path
    assert np.all(values[:, 3:].sum(axis=1) <= 1.00001), path
    # Preserve printed Hartree energies and occupations, p_z and all-orbital
    # lower-layer weights. No interpolation or occupation weighting.
    extracted = np.column_stack((values[:, 1:3], values[:, 5], values[:, 3:].sum(axis=1)))
    return fermi, extracted


def audit(source, destination):
    destination.mkdir(parents=True, exist_ok=True)
    spectra_dir = destination / "spectra"
    spectra_dir.mkdir(exist_ok=True)
    rows, bundles, settings = [], {}, {}
    candidates = sorted((source / "Results_from_cluster").glob("*/*/*angle*"))
    for directory in candidates:
        if not directory.is_dir():
            continue
        series = "/".join(directory.relative_to(source / "Results_from_cluster").parts[:2])
        index, angle = directory.name.split("_angle_")
        input_path = directory / "input_proj_dos.inp"
        output_path = directory / "result.out"
        raw_input = input_path.read_text() if input_path.exists() else ""
        active_input = "\n".join(line.split("!")[0].split("#")[0] for line in raw_input.splitlines())
        output = output_path.read_text(errors="replace") if output_path.exists() else ""
        row = dict(series=series, index=int(index), angle_deg=float(angle),
                   input_present=input_path.exists(), output_present=output_path.exists(),
                   cp2k_version=match(output, r"CP2K\| version string:.*CP2K version (\S+)"),
                   run_type=match(active_input, r"^\s*RUN_TYPE\s+(\S+)"),
                   basis=match(active_input, r"^\s*BASIS_SET\s+(\S+)"),
                   scf_converged=bool(re.search(r"SCF run converged in", output)),
                   program_ended="PROGRAM ENDED AT" in output,
                   scf_steps=match(output, r"SCF run converged in\s+(\d+)"),
                   multiplicity=match(output, r"DFT\| Multiplicity\s+(\d+)"),
                   relative_cutoff_ha=match(output, r"Relative density cutoff \[a.u.\]:\s+([\d.]+)"),
                   input_sha256=sha256(input_path), output_sha256=sha256(output_path))
        xyz_path = directory / "atoms.xyz"
        if xyz_path.exists():
            atom_lines = xyz_path.read_text().splitlines()
            atoms = [line.split() for line in atom_lines[2:] if line.strip()]
            assert len(atoms) == int(atom_lines[0]), xyz_path
            xyz = np.array([[float(v) for v in line[1:4]] for line in atoms])
            elements = np.array([line[0] for line in atoms])
            lower = np.isclose(xyz[:, 2], xyz[:, 2].min(), atol=1e-6)
            row.update(atoms=len(atoms), lower_c=int(np.sum(lower & (elements == "C"))),
                       upper_c=int(np.sum(~lower & (elements == "C"))),
                       hydrogen=int(np.sum(elements == "H")),
                       separation_angstrom=float(np.ptp(xyz[:, 2])),
                       xyz_sha256=sha256(xyz_path))
        for keyword in ("CUTOFF", "NGRIDS", "EPS_DEFAULT", "EPS_SCF", "MAX_SCF",
                        "SCF_GUESS", "ADDED_MOS", "ELECTRONIC_TEMPERATURE", "PERIODIC"):
            value = match(active_input, rf"^\s*{keyword}\s+(.+)$")
            settings.setdefault(keyword, set()).add(value)
            row["input_" + keyword.lower()] = value
        settings.setdefault("basis", set()).add(row["basis"])
        settings.setdefault("version", set()).add(row["cp2k_version"])
        settings.setdefault("run_type", set()).add(row["run_type"])
        mos = [int(v) for v in re.findall(r"Number of molecular orbitals:\s+(\d+)", output)]
        electrons = re.findall(r"Number of electrons:\s+(\d+)", output)
        bundle = bundles.setdefault(series, {})
        for spin_index, spin in enumerate(("alpha", "beta")):
            path = directory / f"graphene_pdos-{spin.upper()}_list1-1.pdos"
            row[f"{spin}_present"] = path.exists() and path.stat().st_size > 0
            if not row[f"{spin}_present"]:
                continue
            try:
                fermi, values = read_pdos(path)
            except (ValueError, AssertionError):
                row[f"{spin}_parse_error"] = True
                continue
            energies = HARTREE_EV * (values[:, 0] - fermi)
            row.update({f"{spin}_rows": len(values), f"{spin}_fermi_ha": fermi,
                        f"{spin}_emin_ev": float(energies.min()),
                        f"{spin}_emax_ev": float(energies.max()),
                        f"{spin}_sha256": sha256(path),
                        f"{spin}_output_mos": mos[spin_index] if len(mos) > spin_index else "",
                        f"{spin}_electrons": electrons[spin_index] if len(electrons) > spin_index else ""})
            row[f"{spin}_rows_match_output"] = len(mos) > spin_index and mos[spin_index] == len(values)
            key = f"i{int(index):03d}_{spin}"
            bundle[key] = values
            bundle[key + "_fermi_ha"] = np.array(fermi)
            bundle[key + "_angle_deg"] = np.array(float(angle))
        rows.append(row)

    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (destination / "calculations.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for series, bundle in bundles.items():
        np.savez_compressed(spectra_dir / (series.replace("/", "__") + ".npz"), **bundle)

    summary = {"conversion_ha_to_ev": HARTREE_EV,
               "spectral_columns": ["eigenvalue_ha", "occupation", "lower_pz_weight", "lower_all_orbital_weight"],
               "settings_present_in_inputs": {key: sorted(values) for key, values in settings.items()},
               "series": {}}
    for series in sorted(bundles):
        group = [row for row in rows if row["series"] == series]
        displayed = [row for row in group if 30 <= row["angle_deg"] <= 60]
        complete = [row for row in displayed if row["scf_converged"] and row["program_ended"]
                    and row.get("alpha_rows_match_output", False)
                    and row.get("beta_rows_match_output", False)]
        complete_alpha = [row for row in displayed if row["scf_converged"] and row["program_ended"]
                          and row.get("alpha_rows_match_output", False)]
        cutoffs = [row["alpha_emax_ev"] for row in complete_alpha]
        summary["series"][series] = dict(
            directories=len(group), displayed_interval_directories=len(displayed),
            converged_with_alpha_pdos=len(complete_alpha),
            converged_with_two_pdos=len(complete),
            incomplete_angles=sorted(row["angle_deg"] for row in displayed if row not in complete),
            available_angles_deg=sorted(row["angle_deg"] for row in group),
            missing_survey_angles_deg=[round(30 + 0.6 * index, 1) for index in range(51)
                                      if not any(abs(row["angle_deg"] - (30 + 0.6 * index)) < 1e-8 for row in displayed)],
            all_distances=sorted({round(row["separation_angstrom"], 6) for row in group if "separation_angstrom" in row}),
            printed_alpha_upper_ev=[min(cutoffs), max(cutoffs)] if cutoffs else [],
            multiplicities=sorted({row["multiplicity"] for row in complete}),
            relative_cutoff_ha=sorted({row["relative_cutoff_ha"] for row in complete}))
    # The monolayer files contain binned total DOS, not state-resolved PDOS.
    mono, mono_metadata = {}, {}
    for path in sorted((source / "Single_layer").glob("*/*.dos")):
        header = path.read_text().splitlines()[0]
        fermis = [float(v) for v in header.split("=")[-1].split()]
        assert len(fermis) == 2, path
        key = path.parent.name + "_" + re.search(r"ideal_(\d+)", path.name).group(1)
        values = np.loadtxt(path, skiprows=2)
        mono[key] = values
        mono[key + "_fermi_ha"] = np.array(fermis)
        mono_metadata[key] = dict(sha256=sha256(path), fermi_ha=fermis, rows=len(values),
                                  maximum_relative_energy_ev=float(HARTREE_EV * (values[-1, 0] - fermis[0])))
    np.savez_compressed(destination / "monolayer_spectra.npz", **mono)
    summary["monolayers"] = mono_metadata
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    table = [r"\begin{tabular}{llcc}", r"Basis & $d$ (\AA) & Audited $\alpha$ columns / 51 & Printed upper energy (eV) \\", r"\hline"]
    for basis, distances in (("SZV", (3.30, 3.00, 2.84, 2.70, 2.55, 2.30, 2.25, 2.20, 2.15, 2.00)),
                             ("DZVP", (2.25, 2.20, 2.15))):
        for distance in distances:
            matching = [value for key, value in summary["series"].items()
                        if key.startswith(basis) and value["all_distances"] == [distance]]
            assert len(matching) <= 1, (basis, distance)
            if matching:
                value = matching[0]
                count = str(value["converged_with_alpha_pdos"])
                limits = value["printed_alpha_upper_ev"]
                upper = f"{limits[0]:.2f}--{limits[1]:.2f}" if limits else "---"
            else:
                count, upper = "not supplied", "---"
            table.append(f"{basis} & {distance:.2f} & {count} & {upper} " + r"\\")
    table.append(r"\end{tabular}")
    (destination / "coverage_table.tex").write_text("\n".join(table) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="Author-supplied folder containing Results_from_cluster and Single_layer")
    parser.add_argument("--output", type=Path, default=ROOT / "data/dft_audit")
    arguments = parser.parse_args()
    audit(arguments.source, arguments.output)
