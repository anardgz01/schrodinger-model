import schrodinger_voluntario as schrodinger
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from itertools import product

NS = [500, 1000, 2000]
LAMBDAS = [0.1, 0.3, 0.5, 1, 5, 10]
worker_cores = 12

data = np.zeros((len(LAMBDAS)+1, len(NS)+1))
data[1:, 0] = LAMBDAS
data[0, 1:] = NS

def run_simulations(n_lamb_tuple):
    n, lamb = n_lamb_tuple
    result = schrodinger.simulate(n, lamb)
    print(f"Simulation for n={n} and lambda={lamb} finished with T. coefficient {result}.")
    return (n, lamb, result)

with ProcessPoolExecutor(max_workers=worker_cores) as executor:
    results = list(executor.map(run_simulations, product(NS, LAMBDAS)))

for n, lamb, result in results:
    data[LAMBDAS.index(lamb)+1, NS.index(n)+1] = result

np.save('schrodinger_voluntario_results.npy', data)