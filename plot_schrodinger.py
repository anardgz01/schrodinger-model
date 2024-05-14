import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

def myplot():
    norms = np.load('resultados/norms.npy')
    wave_function = np.load('resultados/wave_function.npy')
    # norms = np.load('resultados/norms_t_100000.npy')
    # wave_function = np.load('resultados/wave_function_t_100000.npy')
    wave_function_norm = np.abs(wave_function)**2
    
    # Create an array for the time points
    time = np.arange(len(norms))
    position = np.arange(wave_function.shape[1])
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Plot the norms on the first subplot
    ax1.plot(time, norms)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Norm')
    ax1.set_xlim(np.min(time), np.max(time)+1)
    ax1.set_ylim(0.5,1.5)
    print(norms)

    # Initialize the line object for the wave function plot on the second subplot
    line, = ax2.plot(position, np.abs(wave_function_norm[0, :]))
    ax2.set_xlabel('Position')
    ax2.set_ylabel('Wave Function')
    ax2.set_xlim(np.min(position), np.max(position)+1)
    ax2.set_ylim(np.min(wave_function_norm),np.max(wave_function_norm)+0.005)

    # Initialize the text annotation for the current time
    time_text = ax2.text(0.02, 0.95, '', transform=ax2.transAxes)

    # Animation update function
    def update(i):
        line.set_ydata(np.abs(wave_function_norm[i, :]))  # update the y-data of the line object
        time_text.set_text(f'Time = {i}')  # update the text annotation
        return line, time_text,

    # Create the animation with double the default framerate
    ani = FuncAnimation(fig, update, frames=len(norms), blit=True, interval=1000 / (0.3 * len(norms)))

    # Save the animation as a video file
    # writer = PillowWriter(fps=60)
    # ani.save("wave_function_ani_t_100000.gif", writer=writer)

    # Show the plot
    plt.show()

def observableplot():
    x_avg = np.load('resultados/voluntario/observable/x_avg.npy')
    x2_avg = np.load('resultados/voluntario/observable/x2_avg.npy')
    p_avg = np.load('resultados/voluntario/observable/p_avg.npy')
    E_avg = np.load('resultados/voluntario/observable/E_avg.npy')
    K_avg = np.load('resultados/voluntario/observable/K_avg.npy')
    V_avg = np.load('resultados/voluntario/observable/V_avg.npy')

    observables = [x_avg, x2_avg, p_avg, E_avg, K_avg, V_avg]

    x_avg_medida = np.load('resultados/voluntario/observable/medida/x_avg.npy')
    x2_avg_medida = np.load('resultados/voluntario/observable/medida/x2_avg.npy')
    p_avg_medida = np.load('resultados/voluntario/observable/medida/p_avg.npy')
    E_avg_medida = np.load('resultados/voluntario/observable/medida/E_avg.npy')
    K_avg_medida = np.load('resultados/voluntario/observable/medida/K_avg.npy')
    V_avg_medida = np.load('resultados/voluntario/observable/medida/V_avg.npy')

    observables_medida = [x_avg_medida, x2_avg_medida, p_avg_medida, E_avg_medida, K_avg_medida, V_avg_medida]

    time = np.arange(len(x_avg))

    fig_x, axes_x = plt.subplots(2)
    fig_x2, axes_x2 = plt.subplots(2)
    fig_p, axes_p = plt.subplots(2)
    fig_E, axes_E = plt.subplots(2)
    fig_K, axes_K = plt.subplots(2)
    fig_V, axes_V = plt.subplots(2)

    axes = [axes_x, axes_x2, axes_p, axes_E, axes_K, axes_V]
    names = ['<x>', '<x^2>', '<p>', '<E>', '<K>', '<V>']

    for ax, observable, observable_medida, name in zip(axes, observables, observables_medida, names):
        observable = np.real(observable)
        observable_medida = np.real(observable_medida)

        ax[0].plot(time, observable)
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel(name)
        ax[0].set_xlim(np.min(time), np.max(time)+1)
        ax[0].set_ylim(np.min(observable), np.max(observable)+0.005)
        ax[0].set_title(f'Temporal evolution of {name} without measurement')

        ax[1].plot(time, observable_medida)
        ax[1].set_xlabel('Time')
        ax[1].set_ylabel(name)
        ax[1].set_xlim(np.min(time), np.max(time)+1)
        ax[1].set_ylim(np.min(observable_medida), np.max(observable_medida)+0.005)
        ax[1].set_title(f'Temporal evolution of {name} with measurement')

    plt.show()

    
if __name__ == '__main__':
    observableplot()