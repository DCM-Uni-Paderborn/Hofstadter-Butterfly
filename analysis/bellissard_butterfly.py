"""Refined finite-bilayer butterflies and normalized spectral state counts.

Reuses the archived 41-angle spectra without changing them. Additional angles
are independently diagonalized with the same all-pairs hopping Hamiltonian.
The output is a finite open-disk calculation, with energy in units of t_nn.
"""
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
import shutil
import time

# Bound the CPU use of each of the two independent distance calculations.
for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '3'

import numpy as np
from scipy.linalg import eigvalsh
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'bellissard_butterfly'
SCRATCH = ROOT / 'tmp' / 'butterfly_calculation'
DISTANCES = (0.99, 0.97)
ANGLE_COUNT = 241
SIGMA = 0.02


def calculate(distance):
    target = OUT / f'R25_d{distance:g}_xi0.03_n{ANGLE_COUNT}.npz'
    if target.exists():
        with np.load(target) as result:
            assert result['eigenvalues'].shape == (ANGLE_COUNT, 3020)
            assert np.isfinite(result['eigenvalues']).all()
            assert np.isclose(result['d'], distance)
        print(f'Retained {target.name}', flush=True)
        return str(target)

    with np.load(ROOT / 'data' / 'tb_reproduction' /
                 f'R25_d{distance:g}_xi0.03.npz') as archive:
        points = archive['sites'].copy()
        old_angles = archive['theta'].copy()
        old_spectra = archive['eigenvalues'].copy()
        assert np.isclose(archive['R'], 25.0)
        assert np.isclose(archive['a'], 1.0)
        assert np.isclose(archive['xi'], 0.03)
    angles = np.linspace(0, np.pi / 3, ANGLE_COUNT)
    n = len(points)
    intra = np.exp((1.0 - cdist(points, points)) / 0.03)
    np.fill_diagonal(intra, 0.0)
    intra_moment = np.sum(intra * intra) / n
    spectra = np.full((ANGLE_COUNT, 2 * n), np.nan)
    residuals = np.full(ANGLE_COUNT, np.nan)
    checkpoint = SCRATCH / f'd{distance:g}_n{ANGLE_COUNT}.npz'
    if checkpoint.exists():
        with np.load(checkpoint) as saved:
            assert np.array_equal(saved['theta'], angles)
            spectra[:] = saved['eigenvalues']
            residuals[:] = saved['moment_residual']
    start = time.monotonic()
    reused = 0
    for index, angle in enumerate(angles):
        if np.isfinite(spectra[index]).all() and np.isfinite(residuals[index]):
            continue
        rotation = np.array([[np.cos(angle), -np.sin(angle)],
                             [np.sin(angle), np.cos(angle)]])
        cross = np.exp((1.0 - np.sqrt(
            cdist(points, points @ rotation.T)**2 + distance**2)) / 0.03)
        match = np.flatnonzero(np.isclose(old_angles, angle, atol=1e-13, rtol=0))
        if match.size:
            values = old_spectra[int(match[0])]
            reused += 1
        else:
            matrix = np.block([[intra, cross], [cross.T, intra]])
            values = eigvalsh(matrix, overwrite_a=True,
                             check_finite=False, driver='evr')
        spectra[index] = values
        residuals[index] = abs(np.mean(values**2) - intra_moment
                               - np.sum(cross * cross) / n)
        assert residuals[index] < 1e-10
        assert abs(np.mean(values)) < 1e-11
        assert np.all(np.diff(values) >= 0)
        if index % 20 == 0 or index == ANGLE_COUNT - 1:
            np.savez_compressed(checkpoint, theta=angles, eigenvalues=spectra,
                                moment_residual=residuals)
            print(f'd/a={distance:.2f}: {index+1}/{ANGLE_COUNT}, '
                  f'{time.monotonic()-start:.0f} s', flush=True)
    assert np.isfinite(spectra).all()
    np.savez_compressed(target, theta=angles, eigenvalues=spectra,
                        R=25.0, a=1.0, d=distance, xi=0.03, sites=points,
                        moment_residual=residuals,
                        units='energy / exp(-a/xi); angle in radians')
    print(f'Saved {target.name}; {reused} archived angles reused; '
          f'maximum trace-identity error {residuals.max():.3g}', flush=True)
    return str(target)


def density(energies, grid, sigma):
    return np.mean(np.exp(-0.5 * ((grid[:, None] - energies[None, :])
                                / sigma)**2), axis=1) / (np.sqrt(2*np.pi)*sigma)


def make_figure(paths):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.titlesize': 9, 'axes.labelsize': 9,
                         'legend.fontsize': 8, 'xtick.labelsize': 8,
                         'ytick.labelsize': 8, 'pdf.fonttype': 42,
                         'ps.fonttype': 42})
    grid = np.linspace(-6.0, 6.0, 2401)
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.0), layout='constrained',
                             gridspec_kw={'height_ratios': [1.22, 1.0]})
    inset = axes[1, 1].inset_axes([0.56, 0.12, 0.40, 0.32])
    inset.set(xlim=(2.6, 3.5), ylim=(0.91, 1.005),
              xticks=[2.6, 3.0, 3.4], yticks=[0.92, 0.96, 1.0])
    inset.tick_params(labelsize=6, length=2, pad=1)
    colors = ('#2271B2', '#B23056')
    report = {'angle_count': ANGLE_COUNT, 'angle_step_degrees': 0.25,
              'energy_unit': 't_nn', 'gaussian_sigma': SIGMA,
              'energy_grid_step': float(grid[1]-grid[0]),
              'normalization': 'one state (2N=3020)',
              'cross_section_degrees': 30.0, 'series': []}
    for col, (path, color) in enumerate(zip(paths, colors)):
        with np.load(path) as result:
            angles = np.degrees(result['theta'])
            spectra = result['eigenvalues']
            distance = float(result['d'])
            residual = float(result['moment_residual'].max())
        rho = np.column_stack([density(values, grid, SIGMA) for values in spectra])
        integrals = np.trapezoid(rho, grid, axis=0)
        assert np.max(abs(integrals - 1.0)) < 1e-8
        ax = axes[0, col]
        mesh = ax.pcolormesh(angles, grid, rho, shading='nearest', cmap='magma',
                             norm=Normalize(0, 0.5), rasterized=True)
        ax.set(xlim=(0, 60), ylim=(-6, 6), xlabel=r'Twist angle $\theta$ (deg)',
               ylabel=r'Energy $\epsilon=E/t_{\mathrm{nn}}$',
               xticks=[0, 15, 30, 45, 60])
        ax.set_title(rf'({chr(97+col)}) $d/a={distance:.2f}$', loc='left')
        ax.axvline(30, color='white', linestyle='--', linewidth=0.75, alpha=0.9)
        pick = int(np.argmin(abs(angles - 30.0)))
        assert np.isclose(angles[pick], 30.0)
        values = spectra[pick]
        axes[1, 0].plot(grid, rho[:, pick], color=color, linewidth=1.0,
                        label=rf'$d/a={distance:.2f}$')
        # Exact, unbroadened empirical state count. Values below/above the
        # spectrum are 0/1, and degeneracies count with their multiplicity.
        positions = np.r_[grid[0], values, grid[-1]]
        counts = np.r_[0.0, np.arange(1, len(values)+1)/len(values), 1.0]
        axes[1, 1].step(positions, counts, where='post', color=color,
                        linewidth=1.2, label=rf'$d/a={distance:.2f}$')
        inset.step(positions, counts, where='post', color=color, linewidth=0.9)
        report['series'].append({'d_over_a': distance,
                                'spectral_min': float(spectra.min()),
                                'spectral_max': float(spectra.max()),
                                'maximum_moment_residual': residual,
                                'maximum_DOS_integral_error':
                                    float(np.max(abs(integrals-1))),
                                'maximum_DOS': float(rho.max()),
                                'color_clipped_fraction': float(np.mean(rho > 0.5))})
        np.savez_compressed(OUT / f'DOS_IDS_d{distance:g}.npz',
                            theta_degrees=angles, energy_grid=grid, dos=rho,
                            sigma=SIGMA, cross_section_eigenvalues=values,
                            count_energy=positions, normalized_count=counts)
        print(f'Rendered d/a={distance:.2f}; maximum DOS integral error '
              f'{np.max(abs(integrals-1)):.3g}', flush=True)
    fig.colorbar(mesh, ax=list(axes[0]), label=r'$t_{\mathrm{nn}}\rho_\eta(E)$',
                 shrink=0.93, pad=0.025, extend='max')
    axes[1, 0].set(xlim=(-4.2, 4.2), ylim=(0, 0.5), xlabel=r'Energy $\epsilon$',
                   ylabel=r'$t_{\mathrm{nn}}\rho_\eta(E)$')
    axes[1, 0].set_title(r'(c) DOS at $\theta=30^\circ$', loc='left')
    axes[1, 1].set(xlim=(-4.2, 4.2), ylim=(0, 1), xlabel=r'Energy $\epsilon$',
                   ylabel=r'Normalized state count $F_R(\epsilon)$')
    axes[1, 1].set_title(r'(d) Integrated state count at $\theta=30^\circ$', loc='left')
    for ax in axes[1]:
        ax.legend(frameon=False, loc='upper left')
        ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(ROOT / 'figures' / 'bellissard_butterfly_ids.pdf', dpi=400)
    (ROOT / 'output' / 'pdf').mkdir(parents=True, exist_ok=True)
    shutil.copyfile(ROOT / 'figures' / 'bellissard_butterfly_ids.pdf',
                    ROOT / 'output' / 'pdf' / 'bellissard_butterfly_ids.pdf')
    fig.savefig(ROOT / 'output' / 'bellissard_butterfly_ids.png', dpi=250)
    plt.close(fig)
    (OUT / 'numerical_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print('Saved bellissard_butterfly_ids.pdf and PNG preview', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot-only', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    (ROOT / 'output').mkdir(exist_ok=True)
    if args.plot_only:
        files = [str(OUT / f'R25_d{d:g}_xi0.03_n{ANGLE_COUNT}.npz')
                 for d in DISTANCES]
    else:
        with ProcessPoolExecutor(max_workers=2) as pool:
            files = list(pool.map(calculate, DISTANCES))
    make_figure(files)
