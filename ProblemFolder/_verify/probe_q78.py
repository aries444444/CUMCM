import os, sys
sys.path.insert(0, os.getcwd())
import numpy as np
import q78_sensitivity as s
NZ, DZ, dt, K = 81, 1.0/80.0, 2.5, 1000
T = np.full(NZ, s.T0); C = np.full(NZ, s.C0)
R = s.R0
snap = []
for k in range(K):
    res = s.step(T, C, k*dt, (k+1)*dt, dt, NZ, DZ, R, s.rho4, s.cp4, s.k4, s.D4)
    T, C = res[0], res[1]
    if (k+1) % 50 == 0:
        snap.append(np.concatenate([T, C]))
np.savez('probe_q78.npz', T=T, C=C, snap=np.array(snap))
print('q78 done:', T[0], T[-1], C[0], C[-1])
