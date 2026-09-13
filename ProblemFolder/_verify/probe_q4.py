import os, sys
sys.path.insert(0, os.getcwd())
import numpy as np
import q4_solver as s
dt, K = 5.0, 1000
T = np.full(s.NZ, s.T0); C = np.full(s.NZ, s.C0)
snap = []
for k in range(K):
    res = s.step(T, C, k*dt, (k+1)*dt, dt)
    T, C = res[0], res[1]
    if (k+1) % 50 == 0:
        snap.append(np.concatenate([T, C]))
np.savez('probe_q4.npz', T=T, C=C, snap=np.array(snap))
print('q4 done:', T[0], T[-1], C[0], C[-1])
