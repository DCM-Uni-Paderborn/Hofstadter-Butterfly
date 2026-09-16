"""Finite bilayer spectral moments for the Bellissard perspective.

Uses the archived, dimensionless eigenvalues and independently evaluates the
second moment from the hopping matrix. No gap label or invariant is computed.
"""
from pathlib import Path
import csv
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
rows = []
curves = []
max_error = 0.0
for distance in (0.99, 0.98, 0.97):
    with np.load(root / 'data' / 'tb_reproduction' /
                 f'R25_d{distance:g}_xi0.03.npz') as result:
        points = result['sites']
        angles = result['theta']
        eigenvalues = result['eigenvalues']
        a, xi = float(result['a']), float(result['xi'])
        assert np.isclose(float(result['d']), distance)
        assert np.isclose(float(result['R']), 25.0)
    n = len(points)
    intra = np.exp((a - cdist(points, points)) / xi)
    np.fill_diagonal(intra, 0.0)
    moment_intra = np.sum(intra * intra) / n
    cross_moments = []
    for theta, spectrum in zip(angles, eigenvalues):
        rotation = np.array([[np.cos(theta), -np.sin(theta)],
                             [np.sin(theta), np.cos(theta)]])
        separations = np.sqrt(cdist(points, points @ rotation.T)**2
                              + distance**2)
        cross = np.exp((a - separations) / xi)
        moment_cross = np.sum(cross * cross) / n
        moment_spectrum = float(np.mean(spectrum**2))
        error = abs(moment_spectrum - moment_intra - moment_cross)
        max_error = max(max_error, error)
        assert error < 1e-10, (distance, theta, error)
        assert spectrum.size == 2 * n
        assert abs(float(np.mean(spectrum))) < 1e-12
        assert np.isfinite(spectrum).all()
        cross_moments.append(moment_cross)
        rows.append((distance, xi, np.degrees(theta), n,
                     moment_intra, moment_cross, moment_spectrum, error))
    curves.append((distance, np.degrees(angles), np.array(cross_moments)))

destination = root / 'data' / 'bellissard_model_moments.csv'
with destination.open('w', newline='') as stream:
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(('d_over_a', 'xi_over_a', 'theta_degrees', 'sites_per_layer',
                     'intralayer_second_moment', 'interlayer_second_moment',
                     'total_second_moment', 'identity_absolute_error'))
    writer.writerows(rows)

plt.rcParams.update({'font.size': 10, 'axes.labelsize': 10,
                     'legend.fontsize': 9, 'pdf.fonttype': 42,
                     'ps.fonttype': 42})
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.65), layout='constrained')
colors = ('#2271B2', '#E69F00', '#B23056')
distances = np.linspace(0.955, 1.02, 300)
axes[0].plot(distances, np.exp((1.0 - distances) / 0.03), color='#303030')
for (distance, angles, cross_moments), color in zip(curves, colors):
    axes[0].scatter([distance], [np.exp((1.0 - distance) / 0.03)],
                    color=color, s=32, zorder=3)
    axes[1].plot(angles, cross_moments, '.-', markersize=3, linewidth=1,
                 color=color, label=rf'$d/a={distance:.2f}$')
axes[0].axhline(1.0, color='0.7', linewidth=0.7, linestyle='--')
axes[0].set(xlabel=r'Layer separation $d/a$',
            ylabel=r'$t_\perp^{\mathrm{max}}/t_{\mathrm{nn}}$',
            xlim=(0.955, 1.02), ylim=(0, 4.6))
axes[1].set(xlabel=r'Twist angle $\theta$ (degrees)',
            ylabel=r'Interlayer contribution to $\mu_2$',
            xlim=(0, 60), ylim=(0, None))
axes[1].legend(frameon=False, loc='upper right')
for letter, axis in zip(('a', 'b'), axes):
    axis.text(0.0, 1.03, f'({letter})', transform=axis.transAxes,
              ha='left', va='bottom', fontweight='bold')
    axis.spines[['top', 'right']].set_visible(False)
    axis.tick_params(direction='out')
fig.savefig(root / 'figures' / 'bellissard_model_moments.pdf')
plt.close(fig)
print(f'{len(rows)} spectra checked, {n} sites/layer')
print(f'Intralayer second moment: {moment_intra:.12g}')
print(f'Maximum second-moment identity error: {max_error:.3g}')
print(f'Saved {destination.name} and bellissard_model_moments.pdf')
