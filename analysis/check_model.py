"""Numerical consistency tests; not a topology or bulk-convergence test."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial.distance import cdist
from scipy.linalg import eigvalsh

root = Path(__file__).resolve().parents[1]
v1 = np.array([np.sqrt(3), 0.])
v2 = np.array([np.sqrt(3)/2, 1.5])
points = np.array([i*v1+j*v2+b for i in range(-8,9) for j in range(-8,9)
                   for b in [np.array([0.,0.]),np.array([0.,1.])]
                   if np.linalg.norm(i*v1+j*v2+b) <= 3.5])
xi, d, angle = .03, .97, .317

def matrix(theta, normalized):
    rot = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
    sites = np.vstack([np.c_[points, np.zeros(len(points))],
                       np.c_[points @ rot.T, np.full(len(points), d)]])
    distance = cdist(sites, sites)
    h = np.exp(((1 if normalized else 0)-distance)/xi)
    np.fill_diagonal(h,0.)
    return h

h = matrix(angle, True)
tnn = np.exp(-1/xi)
original_scaled = matrix(angle, False)/tnn
rescaling_error = float(np.max(np.abs(h-original_scaled)))
mirror_error = float(np.max(np.abs(eigvalsh(h)-eigvalsh(matrix(-angle,True)))))
assert np.allclose(h, h.T, atol=0, rtol=0)
assert rescaling_error < 1e-12
assert mirror_error < 1e-11
epsilon = np.linspace(-2,2,1001)
eta = .13
scaled_gaussian = np.exp(-.5*(epsilon/eta)**2)/(np.sqrt(2*np.pi)*eta)
original_gaussian = np.exp(-.5*((epsilon*tnn)/(eta*tnn))**2)/(np.sqrt(2*np.pi)*eta*tnn)
jacobian_error = float(np.max(np.abs(scaled_gaussian-tnn*original_gaussian)))
assert jacobian_error < 1e-12
arrays=[]
for path in sorted((root/'data/tb_reproduction').glob('*.npz')):
    data=np.load(path)
    ev=data['eigenvalues']
    assert ev.shape == (41,3020)
    assert np.isfinite(ev).all()
    assert (np.diff(ev,axis=1)>=0).all()
    trace_error=float(np.max(np.abs(ev.sum(axis=1))))
    assert trace_error < 1e-8
    arrays.append({'file':path.name,'shape':list(ev.shape),'maximum_trace_error':trace_error,
                   'minimum_eigenvalue':float(ev.min()),'maximum_eigenvalue':float(ev.max())})
report={'small_test_matrix_dimension':len(h),'energy_rescaling_error':rescaling_error,
        'angle_reflection_eigenvalue_error':mirror_error,'dos_jacobian_error':jacobian_error,
        'verified_arrays':arrays,
        'scope':'Consistency only; no DFT validation, topology, or finite-size convergence.'}
(root/'data/model_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
