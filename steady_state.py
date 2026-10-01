import numpy as np
import matplotlib.pyplot as plt


K, k, h, r = 1, 8, 1, 1
N_max = 10000
runs = 2000
t_run = 200

# Calculating the rates
def tplus(n):
    return k

def tminus(n):
    return n*r + K**h/(K**h + n**h)


logphi = np.zeros(N_max+1)

for n in range(1, N_max+1):
    logphi[n] = logphi[n-1] + np.log(tminus(n)) - np.log(tplus(n))

P = np.exp(-logphi)
P /= np.sum(P)

# Using the Gillespie Algorithm
def run_ssa():
    t, n = 0, 0
    while t < t_run:
        a = tplus(n) + tminus(n)
        t += -np.log(np.random.rand())/a
        
        if np.random.rand() < tplus(n)/a:
            n += 1
        else:
            if n > 0:
                n -= 1
    return n

samples = np.array([run_ssa() for _ in range(runs)])


n = np.arange(20)

plt.plot(n, P[:20], label="Theory")
plt.scatter(n, [np.mean(samples==i) for i in n], label="SSA")
plt.xlabel("Protein number")
plt.ylabel("P(n)")
plt.legend()
plt.show()
