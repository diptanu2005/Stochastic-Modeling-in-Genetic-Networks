import numpy as np


K, k, h, r = 1, 8, 1, 1
N_max = 10000

# Rates
def tplus(n):
    return k

def tminus(n):
    return n*r + K**h/(K**h + n**h)

def phi(n):
    prod = 1.0
    for z in range(1, n+1):
        prod *= tminus(z)/tplus(z)
    return prod

#steady-state P(n)
P = np.zeros(N_max+1)

for n in range(N_max+1):
    P[n] = 1/phi(n)

P /= np.sum(P)

mean = np.sum(np.arange(N_max+1)*P)
target = int(mean)

#MFPT (Tn)
def Tn(final):
    total = 0.0
    for y in range(final):
        inner = 0.0
        for z in range(y+1):
            inner += 1/(tplus(z)*phi(z))
        total += phi(y)*inner
    return total

mfpt = Tn(target)

print("Mean =", mean)
print("Target =", target)
print("Analytical MFPT =", mfpt)
