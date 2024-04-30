# Resolver la ecuación de Schrödinger unidimensional para un potencial cuadrado. Comprobar que
# se conserva la norma

import numpy as np

#Initial parameters
N = 1000
N_CICLOS = 100
VAR_LAMBDA = 1.0
TIEMPO =1000

#Generate s, k0, Vj, Phi_j,0 (phi_0_0=phi_N_0=0) y alpha

K0 = 2*np.pi*N_CICLOS/N
S = 1 /(4*K0) 
V_j = np.zeros(N)
Phi_j = np.zeros(TIEMPO,N)
norms = np.zeros(TIEMPO)
alpha_j = np.zeros(N-1)
beta_j = np.zeros(TIEMPO,N-1)
chi_j = np.zeros(TIEMPO,N-1)
gamma_j = np.zeros(N-1)
A_J_MINUS = 1 
A_J_PLUS = 1
X_INF = 2*N/5
X_SUP = 3*N/5

for j in range(N+1):
    if X_INF < j < X_SUP:
        V_j[j] = VAR_LAMBDA * K0**2

for j in range(1, N):
    Phi_j[0,j] = np.exp(1j*K0*j)*np.exp(-8*(4*j-N)**2/N**2)

A_J_0 = -2 + 2j/S - V_j
for j in range(N-2, -1, -1):
    gamma_inv = A_J_0[j+1] + A_J_PLUS * alpha_j[j+1]
    gamma_j[j+1] =  1/gamma_inv
    alpha_j[j] = -A_J_MINUS*gamma_j[j+1]

for n in range(TIEMPO):

#Calculate beta from equation 22.
    for j in range(N-2, -1, -1):
        beta_j[n,j] = gamma_j[j+1]*(4j*Phi_j[n,j+1]/S - A_J_PLUS*beta_j[n,j+1])

#Calulate chi from equation 20.
    for j in range(1, N):
        chi_j[n,j] = alpha_j[n,j-1]*chi_j[n,j-1] + beta_j[n,j-1]

#Calculate Φj,n+1 from equation 15.
    for j in range (1,N):
        Phi_j[n+1,j] = chi_j[n,j] - Phi_j[n,j]

#Compute the conservation of the norm.
for n in range (TIEMPO):
    norms[n] = np.sum(np.abs(Phi_j[n])**2)
