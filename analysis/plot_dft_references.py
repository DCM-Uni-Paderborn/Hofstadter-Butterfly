"""Render supplied monolayer DOS using common broadening and header Fermi levels.

This is postprocessing of existing CP2K output, not a new model calculation.
The figure shows shapes with individual maximum normalization, not absolute DOS.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
HARTREE_EV = 27.211384


def monolayer_figure(folder, output):
    source = np.load(folder / "monolayer_spectra.npz", allow_pickle=False)
    sigma_ev = 0.10
    curves, checks = {}, {}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "axes.linewidth": 0.6})
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.65), sharey=True, layout="constrained")
    colors = {10: "#0072B2", 20: "#D55E00", 25: "#009E73", 30: "#CC79A7"}
    for ax, basis in zip(axes, ("szv", "dzvp")):
        for size in (10, 20, 25, 30):
            key = f"{basis}_{size}"
            if key not in source:
                continue
            values = source[key]
            fermis = source[key + "_fermi_ha"]
            assert np.allclose(values[:, 0], values[:, 3], atol=1e-10, rtol=0)
            assert fermis[0] == fermis[1]
            energy = HARTREE_EV * (values[:, 0] - fermis[0])
            step = float(np.mean(np.diff(energy)))
            assert np.allclose(np.diff(energy), step, atol=1e-8, rtol=0)
            # Both spin columns are densities on the same uniform CP2K grid.
            density = values[:, 1] + values[:, 4]
            smooth = gaussian_filter1d(density, sigma_ev / step, mode="constant", truncate=8)
            window = (energy >= -7) & (energy <= 3)
            maximum = float(smooth[window].max())
            plotted = smooth / maximum
            curves[key] = np.column_stack((energy, plotted))
            checks[key] = dict(grid_spacing_ev=step, gaussian_sigma_ev=sigma_ev,
                               fermi_ha=float(fermis[0]), manual_shift_ev=0,
                               maximum_in_display_window=float(plotted[window].max()))
            ax.plot(energy, plotted, lw=0.95, color=colors[size], label=f"{size}")
        ax.set_xlim(-7, 3)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel(r"$E-E_F$ (eV)")
        ax.set_title(f"({'a' if basis == 'szv' else 'b'}) {basis.upper()}", loc="left")
        ax.axvline(0, color="0.5", lw=0.5, ls="--", zorder=0)
        ax.legend(title="Size label", frameon=False, fontsize=8, title_fontsize=8,
                  ncol=2, loc="upper right")
    axes[0].set_ylabel("Normalized total DOS")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)
    np.savez_compressed(folder / "monolayer_curves.npz", **curves)
    (folder / "monolayer_plot_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data/dft_audit")
    parser.add_argument("--figure", type=Path, default=ROOT / "figures/monolayer_reference.pdf")
    args = parser.parse_args()
    monolayer_figure(args.data, args.figure)
