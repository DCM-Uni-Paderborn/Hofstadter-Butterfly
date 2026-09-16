"""Independent reproduction of the supplied Efremkin/Prodan notebook.

The Hamiltonian is divided by the nearest-neighbor intralayer hopping.
This is an exact scalar rescaling, not a fit to DFT. Only eigenvalues are
computed. No topological invariant is inferred from these data.
"""
from pathlib import Path
import argparse
import time
import numpy as np
from scipy.linalg import eigvalsh
from scipy.spatial.distance import cdist

parser = argparse.ArgumentParser()
parser.add_argument('--radius', type=float, default=25.0)
parser.add_argument('--distances', type=float, nargs='+', default=[0.99, 0.98, 0.97])
parser.add_argument('--xi', type=float, default=0.03)
parser.add_argument('--angles', type=int, default=41)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
out = root / 'data' / 'tb_reproduction'
out.mkdir(parents=True, exist_ok=True)
a = 1.0
M = int(2 * args.radius / a) + 2
v1 = np.array([np.sqrt(3) * a, 0.0])
v2 = np.array([np.sqrt(3) * a / 2, 3 * a / 2])
basis = np.array([[0.0, 0.0], [0.0, a]])
P = np.array([i*v1+j*v2+b for i in range(-M, M+1)
              for j in range(-M, M+1) for b in basis
              if np.linalg.norm(i*v1+j*v2+b) <= args.radius])
theta = np.linspace(0, np.pi/3, args.angles)
intra = np.exp((a - cdist(P, P)) / args.xi)
np.fill_diagonal(intra, 0.0)
print(f'{len(P)} sites/layer; {2*len(P)} matrix dimension', flush=True)
for d in args.distances:
    target = out / f'R{args.radius:g}_d{d:g}_xi{args.xi:g}.npz'
    if target.exists():
        print(f'Existing result retained: {target.name}', flush=True)
        continue
    spectra = []
    start = time.monotonic()
    for i, th in enumerate(theta):
        rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        cross = np.exp((a - np.sqrt(cdist(P, P @ rot.T)**2+d*d)) / args.xi)
        H = np.block([[intra, cross], [cross.T, intra]])
        energies = eigvalsh(H, overwrite_a=True, check_finite=False, driver='evr')
        spectra.append(energies)
        if i % 5 == 0:
            print(f'd={d:g}: {i+1}/{len(theta)} angles, {time.monotonic()-start:.1f} s', flush=True)
    np.savez_compressed(target, theta=theta, eigenvalues=np.array(spectra),
                        R=args.radius, a=a, d=d, xi=args.xi, sites=P,
                        units='energy / exp(-a/xi); angle in radians')
    print(f'Saved {target.name}', flush=True)
