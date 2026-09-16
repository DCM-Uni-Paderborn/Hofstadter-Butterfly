"""Reconstruct audited alpha-spin, lower-layer p_z DOS columns from the archive.

The original figure images remain unchanged. This diagnostic replot additionally
masks energies above the last printed state and columns with incomplete records.
"""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
HARTREE_EV = 27.211384


def spectral_column(states, fermi, grid, eta=0.02):
    energy = HARTREE_EV * (states[:, 0] - fermi)
    selected = (energy >= -7) & (energy <= 3)
    energy, weights = energy[selected], states[selected, 2]
    density = np.zeros_like(grid)
    # Omit Gaussian tails smaller than exp(-64). This is numerical acceleration
    # of the supplied kernel, not a new broadening convention or histogram.
    for center, weight in zip(energy, weights):
        lower, upper = np.searchsorted(grid, (center - 8 * eta, center + 8 * eta))
        density[lower:upper] += weight * np.exp(-((grid[lower:upper] - center) / eta) ** 2)
    density /= 2 * np.pi * eta
    probes = np.linspace(0, len(grid) - 1, 31, dtype=int)
    direct = (np.exp(-((grid[probes, None] - energy[None, :]) / eta) ** 2)
              @ weights) / (2 * np.pi * eta)
    assert np.allclose(density[probes], direct, rtol=1e-12, atol=1e-12)
    return density / density.max()


def replot(folder, series, output):
    archive = np.load(folder / "spectra" / (series.replace("/", "__") + ".npz"), allow_pickle=False)
    with (folder / "calculations.csv").open() as stream:
        rows = [row for row in csv.DictReader(stream) if row["series"] == series]
    angles = np.arange(51) * 0.6 + 30
    grid = -7 + np.arange(6000) * (10 / 6000)
    intensity = np.full((len(grid), len(angles)), np.nan)
    coverage = np.zeros_like(intensity, dtype=np.uint8)
    # coverage: 0 unavailable column; 1 within printed states; 2 above last state.
    for row in rows:
        angle = float(row["angle_deg"])
        columns = np.flatnonzero(np.isclose(angles, angle, rtol=0, atol=1e-8))
        if len(columns) != 1 or not all(row.get(key) == "True" for key in
                                       ("scf_converged", "program_ended", "alpha_rows_match_output")):
            continue
        column = int(columns[0])
        key = f'i{int(row["index"]):03d}_alpha'
        states = archive[key]
        intensity[:, column] = spectral_column(states, float(archive[key + "_fermi_ha"]), grid)
        coverage[:, column] = np.where(grid <= float(row["alpha_emax_ev"]), 1, 2)
    visible = np.ma.masked_where(coverage != 1, intensity)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(5.0, 3.7), layout="constrained")
    ax.set_facecolor("#cccccc")
    image = ax.pcolormesh(angles, grid, visible, vmin=0, vmax=1, cmap="viridis", shading="auto", rasterized=True)
    for column in np.flatnonzero(np.all(coverage == 0, axis=0)):
        ax.axvspan(angles[column] - 0.3, angles[column] + 0.3, color="white", lw=0)
    ax.set(xlabel=r"$\theta$ (deg)", ylabel=r"$E-E_F$ (eV)", xlim=(29.7, 60.3), ylim=(-7, 3))
    ax.set_title(series.replace("/", " / ").replace("_", " "), fontsize=10)
    fig.colorbar(image, ax=ax, label=r"Relative lower-layer $p_z$ weight ($\alpha$)")
    fig.text(0.5, 0.01, "White: unavailable column. Gray: beyond printed states.", ha="center", fontsize=8)
    fig.get_layout_engine().set(rect=(0, 0.05, 1, 0.95))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)
    np.savez_compressed(output.with_suffix(".npz"), energy_ev=grid, angles_deg=angles,
                        intensity=intensity, coverage=coverage)
    print(json.dumps(dict(series=series, audited_columns=int(np.sum(np.any(coverage == 1, axis=0))),
                          gaussian_eta_ev=0.02, gaussian_fwhm_ev=2 * np.sqrt(np.log(2)) * .02)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data/dft_audit")
    parser.add_argument("--series", default="SZV/IL_25_z2-20")
    parser.add_argument("--output", type=Path, default=ROOT / "tmp/dft_replot")
    args = parser.parse_args()
    replot(args.data, args.series, args.output)
