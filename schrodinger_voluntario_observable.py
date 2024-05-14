# Resolver la ecuación de Schrödinger unidimensional para un potencial cuadrado. Comprobar que
# se conserva la norma

import numpy as np

#Initial parameters
N = 1000
N_CICLOS = 100
VAR_LAMBDA = 0.3
TIME = 1000
m = 1
h = 1

#Generate s, k0, Vj, Phi_j,0 (phi_0_0=phi_N_0=0) y alpha

K0 = 2*np.pi*N_CICLOS/N
S = 1 /(4*K0) 
V_j = np.zeros(N)
Phi_j = np.zeros((TIME,N), dtype=complex)
norms = np.zeros(TIME)
alpha_j = np.zeros(N-1, dtype=complex)
beta_j = np.zeros((TIME,N-1), dtype=complex)
chi_j = np.zeros((TIME,N-1), dtype=complex)
gamma_j = np.zeros(N-1, dtype=complex)
A_J_MINUS = 1 
A_J_PLUS = 1
X_INF = 2*N/5
X_SUP = 3*N/5

for j in range(N+1):
    if X_INF < j < X_SUP:
        V_j[j] = VAR_LAMBDA * K0**2

for j in range(1, N):
    Phi_j[0,j] = 1/np.sqrt(np.sqrt(np.pi)*N/16)*np.exp(1j*K0*j)*np.exp(-8*(4*j-N)**2/N**2)

A_J_0 = -2 + 2j/S - V_j
for j in range(N-3, -1, -1):
    gamma_inv = A_J_0[j+1] + A_J_PLUS * alpha_j[j+1]
    gamma_j[j+1] =  1/gamma_inv
    alpha_j[j] = -A_J_MINUS*gamma_j[j+1]

for n in range(TIME-1):

#Calculate beta from equation 22.
    for j in range(N-3, -1, -1):
        beta_j[n,j] = gamma_j[j+1]*(4j*Phi_j[n,j+1]/S - A_J_PLUS*beta_j[n,j+1])

#Calulate chi from equation 20.
    for j in range(1, N-1):
        chi_j[n,j] = alpha_j[j-1]*chi_j[n,j-1] + beta_j[n,j-1]

#Calculate Φj,n+1 from equation 15.
    for j in range (1,N-1):
        Phi_j[n+1,j] = chi_j[n,j] - Phi_j[n,j]

#Compute the conservation of the norm.
for n in range (TIME):
    norms[n] = np.sum(np.abs(Phi_j[n])**2)

# Calculate the derivatives of the wave function
dPhi_j = np.gradient(Phi_j, axis=1)
ddPhi_j = np.gradient(dPhi_j, axis=1)

# Calculate the average position, average position squared, average momentum, average energy, average kinetic energy and average potential energy.
x_avg = np.sum(np.abs(Phi_j)**2 * np.arange(N), axis=1)
x2_avg = np.sum(np.abs(Phi_j)**2 * np.arange(N)**2, axis=1)
p_avg = np.sum(np.conj(Phi_j) * (-1j * h) * dPhi_j, axis=1)
E_avg = np.sum(np.conj(Phi_j) * (-h**2 / (2 * m)) * ddPhi_j + V_j * np.abs(Phi_j)**2, axis=1)
K_avg = np.sum(np.conj(Phi_j) * (-h**2 / (2 * m)) * ddPhi_j, axis=1)
V_avg = np.sum(V_j * np.abs(Phi_j)**2, axis=1)

# Save the results
np.save('resultados/voluntario/observable/x_avg.npy', x_avg)
np.save('resultados/voluntario/observable/x2_avg.npy', x2_avg)
np.save('resultados/voluntario/observable/p_avg.npy', p_avg)
np.save('resultados/voluntario/observable/E_avg.npy', E_avg)
np.save('resultados/voluntario/observable/K_avg.npy', K_avg)
np.save('resultados/voluntario/observable/V_avg.npy', V_avg)
np.save('resultados/voluntario/observable/wave_function.npy', Phi_j)
np.save('resultados/voluntario/observable/norms.npy', norms)
