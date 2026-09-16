"""Finite-radius, Gaussian-width, and decoupled-layer controls for the perspective.

Uses exactly the archived honeycomb site construction and normalized hopping.
No bulk extrapolation or topological assignment is made. Existing result files
are validated and retained. Run this before bellissard_butterfly.py --plot-only.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import os

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '3'

import numpy as np
from scipy.linalg import eigvalsh
from scipy.spatial.distance import cdist
from scipy.stats import wasserstein_distance

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'bellissard_robustness'
RADII = (20.0, 25.0, 30.0)
DISTANCES = (0.99, 0.97)
ANGLES = (0.0, 15.0, 30.0)
WIDTHS = (0.01, 0.02, 0.04)
XI = 0.03


def sites(radius):
    extent = int(2 * radius) + 2
    v1 = np.array([np.sqrt(3), 0.0])
    v2 = np.array([np.sqrt(3) / 2, 1.5])
    basis = np.array([[0.0, 0.0], [0.0, 1.0]])
    return np.array([i*v1+j*v2+b for i in range(-extent, extent+1)
                     for j in range(-extent, extent+1) for b in basis
                     if np.linalg.norm(i*v1+j*v2+b) <= radius])


def calculate(radius):
    target = OUT / f'R{radius:g}_controls.npz'
    points = sites(radius)
    n = len(points)
    if target.exists():
        with np.load(target) as result:
            assert np.array_equal(result['sites'], points)
            assert np.array_equal(result['angles_degrees'], ANGLES)
            assert np.array_equal(result['distances'], DISTANCES)
            assert result['eigenvalues'].shape == (2, 3, 2*n)
            assert result['monolayer_eigenvalues'].shape == (n,)
            assert np.isfinite(result['eigenvalues']).all()
            assert float(result['R']) == radius
            assert float(result['a']) == 1.0
            assert float(result['xi']) == XI
        print(f'Retained {target.name}', flush=True)
        return str(target)
    intra = np.exp((1.0 - cdist(points, points)) / XI)
    np.fill_diagonal(intra, 0.0)
    mono = eigvalsh(intra, check_finite=False, driver='evr')
    intra_moment = np.sum(intra*intra) / n
    spectra = np.empty((2, 3, 2*n))
    inter_moments = np.empty((2, 3))
    for di, distance in enumerate(DISTANCES):
        for ti, degrees in enumerate(ANGLES):
            angle = np.radians(degrees)
            rotation = np.array([[np.cos(angle), -np.sin(angle)],
                                 [np.sin(angle), np.cos(angle)]])
            cross = np.exp((1.0 - np.sqrt(
                cdist(points, points @ rotation.T)**2 + distance**2)) / XI)
            inter_moments[di, ti] = np.sum(cross*cross) / n
            if radius == 25:
                with np.load(ROOT / 'data' / 'bellissard_butterfly' /
                             f'R25_d{distance:g}_xi0.03_n241.npz') as archive:
                    assert np.array_equal(points, archive['sites'])
                    index, = np.flatnonzero(np.isclose(
                        archive['theta'], angle, atol=1e-13, rtol=0))
                    values = archive['eigenvalues'][index].copy()
            else:
                matrix = np.block([[intra, cross], [cross.T, intra]])
                values = eigvalsh(matrix, overwrite_a=True,
                                 check_finite=False, driver='evr')
            assert np.all(np.diff(values) >= 0)
            assert abs(np.mean(values)) < 1e-11
            assert abs(np.mean(values**2) - intra_moment
                       - inter_moments[di, ti]) < 1e-10
            spectra[di, ti] = values
            print(f'R/a={radius:g}, N={n}, d/a={distance:g}, '
                  f'theta={degrees:g}: validated', flush=True)
    np.savez_compressed(target, sites=points, R=radius, a=1.0, xi=XI,
                        distances=DISTANCES, angles_degrees=ANGLES,
                        eigenvalues=spectra, monolayer_eigenvalues=mono,
                        intra_second_moment=intra_moment,
                        inter_second_moment=inter_moments,
                        energy_unit='t_nn', normalization='one state')
    return str(target)


def density(values, grid, sigma):
    # Chunking bounds temporary storage even at the largest radius.
    result = np.empty_like(grid)
    for start in range(0, len(grid), 400):
        energies = grid[start:start+400, None]
        result[start:start+400] = np.mean(
            np.exp(-0.5*((energies-values[None, :])/sigma)**2), axis=1
        ) / (np.sqrt(2*np.pi)*sigma)
    return result


def count_distance(left, right):
    """Exact integral |F_left-F_right|, independently checked against W1."""
    knots = np.unique(np.r_[left, right])
    counts_left = np.searchsorted(left, knots[:-1], side='right') / len(left)
    counts_right = np.searchsorted(right, knots[:-1], side='right') / len(right)
    exact = np.sum(np.diff(knots)*abs(counts_left-counts_right))
    assert abs(exact - wasserstein_distance(left, right)) < 1e-12
    return float(exact)


def analyze(paths):
    assert abs(count_distance(np.array([0.0, 1.0]),
                              np.array([0.3, 1.3])) - 0.3) < 1e-14
    grid = np.linspace(-7.0, 7.0, 5601)
    results = []
    checks = {'energy_unit': 't_nn', 'normalization': 'one state',
              'radii': list(RADII), 'angles_degrees': list(ANGLES),
              'distances': list(DISTANCES), 'gaussian_widths': list(WIDTHS),
              'energy_grid_step': float(grid[1]-grid[0]),
              'maximum_DOS_integral_error': 0.0,
              'maximum_DOS_second_moment_error': 0.0,
              'maximum_spectral_moment_residual': 0.0,
              'maximum_decoupled_DOS_identity_error': 0.0,
              'sizes': [], 'size_comparisons': [], 'coupling_comparisons': []}
    for path in paths:
        with np.load(path) as data:
            result = {key: data[key].copy() for key in data.files}
        spectra = result['eigenvalues']
        mono = result['monolayer_eigenvalues']
        n = len(mono)
        assert np.all(np.diff(mono) >= 0) and np.isfinite(mono).all()
        assert abs(mono.mean()) < 1e-11
        residual = float(np.max(abs(np.mean(spectra**2, axis=2)
                             - result['intra_second_moment']
                             - result['inter_second_moment'])))
        assert residual < 1e-10
        assert abs(np.mean(mono**2)-result['intra_second_moment']) < 1e-10
        checks['maximum_spectral_moment_residual'] = max(
            checks['maximum_spectral_moment_residual'], residual)
        checks['sizes'].append({'R_over_a': float(result['R']),
                                'sites_per_layer': n})
        rho = np.empty((2, 3, 3, len(grid)))
        mono_rho = np.empty((3, len(grid)))
        for wi, width in enumerate(WIDTHS):
            mono_rho[wi] = density(mono, grid, width)
            # The 2N-state block-diagonal spectrum duplicates each eigenvalue.
            doubled = density(np.repeat(mono, 2), grid, width)
            error = float(np.max(abs(mono_rho[wi]-doubled)))
            checks['maximum_decoupled_DOS_identity_error'] = max(
                checks['maximum_decoupled_DOS_identity_error'], error)
            for di in range(2):
                for ti in range(3):
                    values = spectra[di, ti]
                    rho[di, ti, wi] = density(values, grid, width)
                    area = np.trapezoid(rho[di, ti, wi], grid)
                    error = abs(float(area)-1.0)
                    checks['maximum_DOS_integral_error'] = max(
                        checks['maximum_DOS_integral_error'], error)
                    moment = np.trapezoid(grid**2*rho[di, ti, wi], grid)
                    error = abs(float(moment)-np.mean(values**2)-width**2)
                    checks['maximum_DOS_second_moment_error'] = max(
                        checks['maximum_DOS_second_moment_error'], error)
        assert np.max(abs(np.trapezoid(mono_rho, grid, axis=-1)-1)) < 1e-10
        result['dos'] = rho
        result['mono_dos'] = mono_rho
        results.append(result)
        np.savez_compressed(OUT / f'R{float(result["R"]):g}_DOS.npz',
                            energy_grid=grid, widths=WIDTHS, dos=rho,
                            monolayer_dos=mono_rho, angles_degrees=ANGLES,
                            distances=DISTANCES,
                            axis_order='distance, angle, width, energy')
    for key in ('maximum_DOS_integral_error', 'maximum_DOS_second_moment_error',
                'maximum_decoupled_DOS_identity_error'):
        assert checks[key] < 1e-10, (key, checks[key])
    for left, right in zip(results[:-1], results[1:]):
        for di, distance in enumerate(DISTANCES):
            for ti, angle in enumerate(ANGLES):
                row = {'R_left': float(left['R']), 'R_right': float(right['R']),
                       'd_over_a': distance, 'angle_degrees': angle,
                       'integrated_count_distance': count_distance(
                           left['eigenvalues'][di, ti], right['eigenvalues'][di, ti])}
                for wi, width in enumerate(WIDTHS):
                    row[f'DOS_L1_sigma_{width:g}'] = float(np.trapezoid(abs(
                        left['dos'][di, ti, wi]-right['dos'][di, ti, wi]), grid))
                checks['size_comparisons'].append(row)
    for result in results:
        for di, distance in enumerate(DISTANCES):
            for ti, angle in enumerate(ANGLES):
                row = {'R_over_a': float(result['R']), 'd_over_a': distance,
                       'angle_degrees': angle, 'coupled_vs_decoupled_count_distance':
                       count_distance(result['eigenvalues'][di, ti],
                                      result['monolayer_eigenvalues'])}
                for wi, width in enumerate(WIDTHS):
                    row[f'DOS_L1_sigma_{width:g}'] = float(np.trapezoid(abs(
                        result['dos'][di, ti, wi]-result['mono_dos'][wi]), grid))
                checks['coupling_comparisons'].append(row)
    checks['comparison_summary'] = []
    checks['grid_refinement_checks'] = []
    for radius in RADII[:-1]:
        rows = [row for row in checks['size_comparisons'] if row['R_left'] == radius]
        summary = {'R_left': radius, 'R_right': rows[0]['R_right'],
                   'count_distance_min': min(row['integrated_count_distance'] for row in rows),
                   'count_distance_max': max(row['integrated_count_distance'] for row in rows)}
        for width in WIDTHS:
            key = f'DOS_L1_sigma_{width:g}'
            summary[f'max_{key}'] = max(row[key] for row in rows)
            worst = max(rows, key=lambda row: row[key])
            di = DISTANCES.index(worst['d_over_a'])
            ti = ANGLES.index(worst['angle_degrees'])
            left = next(result for result in results if result['R'] == radius)
            right = next(result for result in results if result['R'] == worst['R_right'])
            refined_grid = np.linspace(-7.0, 7.0, 11201)
            refined = float(np.trapezoid(abs(
                density(left['eigenvalues'][di, ti], refined_grid, width)
                - density(right['eigenvalues'][di, ti], refined_grid, width)),
                refined_grid))
            error = abs(refined-worst[key])
            assert error < 5e-5
            checks['grid_refinement_checks'].append({
                'R_left': radius, 'R_right': worst['R_right'],
                'd_over_a': worst['d_over_a'],
                'angle_degrees': worst['angle_degrees'], 'sigma': width,
                'DOS_L1_grid_0.0025': worst[key],
                'DOS_L1_grid_0.00125': refined, 'absolute_difference': error})
        checks['comparison_summary'].append(summary)
    for name in ('size_comparisons', 'coupling_comparisons'):
        with (OUT / f'{name}.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(checks[name][0]),
                                    lineterminator='\n')
            writer.writeheader()
            writer.writerows(checks[name])
    (OUT / 'numerical_checks.json').write_text(json.dumps(checks, indent=2)+'\n')
    print(json.dumps({key: value for key, value in checks.items()
                      if key not in ('size_comparisons', 'coupling_comparisons')},
                     indent=2), flush=True)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=2) as pool:
        files = list(pool.map(calculate, RADII))
    analyze(files)
