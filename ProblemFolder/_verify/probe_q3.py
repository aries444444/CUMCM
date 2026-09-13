import os, sys
sys.path.insert(0, os.getcwd())
import numpy as np
import q3_solver as s
dt, K = 2.5, 4000
T = np.full(s.N, s.T0); C = np.full(s.N, s.C0)
snap = []
for k in range(K):
    res = s.step(T, C, k*dt, (k+1)*dt, dt)
    T, C = res[0], res[1]
    if (k+1) % 200 == 0:
        snap.append(np.concatenate([T, C]))
np.savez('probe_q3.npz', T=T, C=C, snap=np.array(snap))
print('q3 done:', T[0], T[-1], C[0], C[-1])
