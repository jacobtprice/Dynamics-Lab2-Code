#%% imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%
L = 18 * 0.0254 #length of rod, in meters
#array of tested frequencies, correlated with the trial number
freqs = [3.7, 3.7, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1, 3.0, 2.9, 2.8, 2.7, 3.8,
        3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7]

points = [] #initialized array of points to be plotted
for i in range(1, 24):
    data = pd.read_csv(f"week2_trial{i}.csv", header = None) #importing
    time = np.array(data.iloc[:,0]) 
    disp = np.array(data.iloc[:,1])*0.0254
    ang = np.arcsin(disp/L)
    ang = ang - np.average(ang) #mean-centering displacement
    
    
    y = abs(np.fft.rfft(disp)) #fourier fast transform of displacement
    y = y*2 / len(disp) 
    m = np.round(max(y), 4) #finding max val in y array
    points.append(m)

#plotting
plt.scatter(x = freqs, y = points, marker = "o", s = 15, label = "Raw Data")
plt.axvline(x=3.7, color = "r", linestyle = "--", alpha = 0.6, 
            label = "Resonance (3.7)") #resonance frequency
plt.grid()
plt.xlabel("Frequency")
plt.ylabel(r"$\frac{{w_n^2}}{{w_n^2 - w^2}}$")
plt.legend()
plt.tight_layout()
plt.savefig("part2_plot.jpg")
plt.show()