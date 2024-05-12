# Quiero calcular el coeficiente de transmision.

import numpy as np

def simulate(n : int = 1000, lamb : float = 0.1):
    '''If called as main script, saves the wave function and the norms.
    If called as a module, function returns the transmission coefficient.'''
    #Initial parameters
    N = n
    N_CICLOS = 100
    VAR_LAMBDA = lamb
    TIME = 2*N
    n_d = N//4
    m_t = 0.0
    m_times = 10**3

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

    #Generate the potential
    for j in range(N+1):
        if X_INF < j < X_SUP:
            V_j[j] = VAR_LAMBDA * K0**2

    #Generate the initial wave function
    def generate_initial_wave_function():
        nonlocal Phi_j
        for j in range(1, N):
            Phi_j[0,j] = 1/np.sqrt(np.sqrt(np.pi)*N/16)*np.exp(1j*K0*j)*np.exp(-8*(4*j-N)**2/N**2)

    #Calculate the coefficients alpha and gamma
    A_J_0 = -2 + 2j/S - V_j
    for j in range(N-3, -1, -1):
        gamma_inv = A_J_0[j+1] + A_J_PLUS * alpha_j[j+1]
        gamma_j[j+1] =  1/gamma_inv
        alpha_j[j] = -A_J_MINUS*gamma_j[j+1]

    #Compute the system 10**3 times
    for i in range(m_times):
        generate_initial_wave_function()

        #Compute the wave function for each time step.
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

            if n % n_d == 0 and n != 0:
                print(f'Iteración {n} de {TIME}')

                # Calculate the probability of being in the right side between j=4*N/5 and j=N
                probability_right = np.sum(np.abs(Phi_j[n, 4*N//5:N])**2)

                #Generate a random number between 0 and 1
                random_number = np.random.rand()
                if random_number < probability_right:
                    m_t += 1
                    print(f'Hey, I am in the right side at time {n}!')
                    break

                Phi_j[n,4*N//5:N] = 0
                k_norm = np.sum(np.abs(Phi_j[n])**2)
                Phi_j[n] = Phi_j[n]/np.sqrt(k_norm)

                #Calculate the probability of being in the left side between j=0 and j=N/5
                probability_left = np.sum(np.abs(Phi_j[n, 0:N//5])**2)

                #Generate a random number between 0 and 1
                random_number = np.random.rand()
                if random_number < probability_left:
                    print(f'Hey, I am in the left side at time {n}!')
                    break

                Phi_j[n,0:N//5] = 0
                k_norm = np.sum(np.abs(Phi_j[n])**2)
                Phi_j[n] = Phi_j[n]/np.sqrt(k_norm)

    #Calculate the transmission coefficient
    T_coefficient = m_t/m_times
    print(f'El coeficiente de transmisión es: {T_coefficient}')

    #Compute the conservation of the norm.
    for n in range (TIME):
        norms[n] = np.sum(np.abs(Phi_j[n])**2)

    if __name__ == '__main__':
        np.save('resultados/voluntario/wave_function.npy', Phi_j)
        np.save('resultados/voluntario/norms.npy', norms)
        return
    
    return T_coefficient

if __name__ == '__main__':
    simulate()