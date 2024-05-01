import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

def myplot():
    norms = np.load('resultados/norms.npy')
    wave_function = np.load('resultados/wave_function.npy')
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
    ax1.set_ylim(100,120)

    # Initialize the line object for the wave function plot on the second subplot
    line, = ax2.plot(position, np.abs(wave_function_norm[0, :]))
    ax2.set_xlabel('Position')
    ax2.set_ylabel('Wave Function')
    ax2.set_xlim(np.min(position), np.max(position)+1)
    ax2.set_ylim(np.min(wave_function_norm),np.max(wave_function_norm)+0.15)

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
    writer = PillowWriter(fps=60)
    ani.save("wave_function_ani_t_100000.gif", writer=writer)

    # Show the plot
    plt.show()
    
if __name__ == '__main__':
    myplot()