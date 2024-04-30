# Resolver la ecuación de Schrödinger unidimensional para un potencial cuadrado. Comprobar que
# se conserva la norma

import numpy as np

# 1. Dar los parámetros iniciales: N, nciclos y λ. Generar s˜,
# ˜k0, V˜j# , Φj,0 (incluyendo las condiciones de contorno
# Φ0,0 = ΦN,0 = 0) y α.

#Initial parameters
N = 1000
N_CICLOS = 100
VAR_LAMBDA = 1.0

#Generate s, k0, Vj, Phi_j,0 (phi_0_0=phi_N_0=0) y alpha

K0 = 2*np.pi*N_CICLOS/N
S = 1 /(4*K0) 
V_j = np.zeros(N)
Phi_j_0 = np.zeros(N)
alpha_j = np.zeros(N-1)
gamma_j = np.zeros(N-1)
A_J_MINUS = 1 
A_J_0 = -2 + 2j/S - V_j
A_J_PLUS = 1
X_INF = 2*N/5
X_SUP = 3*N/5

for j in range(N+1):
    if X_INF < j < X_SUP:
        V_j[j] = VAR_LAMBDA * K0**2

for j in range(1, N):
    Phi_j_0 = np.exp(1j*K0*j)*np.exp(-2*(4*j-N)**2/N**2)

for j in range(N-2, -1, -1):
    gamma_inv = A_J_0[j+1] + A_J_PLUS * alpha_j[j+1]
    gamma_j [j+1] =  np.linalg.inv(gamma_inv)
    alpha_j [j] = -A_J_MINUS*gamma_j [j+1]
    
# 2. Calcular β utilizando la recurrencia (22).
# 3. Calcular χ a partir de (20).
# 4. Calcular Φj,n+1 de (15).
# 5. n = n + 1, ir a al paso 2.