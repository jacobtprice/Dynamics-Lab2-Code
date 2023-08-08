#%% imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.optimize import curve_fit

#%% Calculations for spring constant k
spring = pd.read_csv("displacement_force_lab2.csv", header = 1)
meter = spring.iloc[:,0]*0.0254
force = spring.iloc[:,2]
#creating line of best fit equation
meter_model = np.linspace(np.min(meter), np.max(meter), 1000)
force_model = np.polyfit(meter, force, 1)
k = round(force_model[0], 4)
intercept = force_model[1]
#plotting
plt.plot(meter_model, k*meter_model + intercept, label = f"k ={k}")
plt.scatter(meter, force, marker = "o", label = "Raw Data")
plt.xlabel("Displacement (meters)")
plt.ylabel("Force Reading (newtons)")
plt.grid()
plt.legend()
plt.savefig("force_disp.jpg")
plt.show()

#%%
g = 9.81 #gravity
L = 18 * 0.0254 #length of rod, in meters
washL = L - (0.3*0.0254) # distance from pin to washer location
L_B = 0.2928 #distance from pin to spring attachment
m1 = 0.02688 #mass of washers in systems 2, in kg
m2 = 0.01336 #mass of washers in system 3 and 4, in kg
m_rod = 0.15289 #mass of rod in kg
total_m = [m_rod, m1+m_rod, m2+m_rod, m2+m_rod] #total system mass for each system

rgo = [] 
I0 = []   #initialized center of mass and moment of intertia arrays
#calculating rgo and I0 for each system
for i in range(4):
    if i == 0:
        rgo.append(L/2)
        I0.append(1/3*m_rod*L**2)
    else:
        if i == 1:
            rgo.append((m_rod*(L/2)+m1*washL)/(m_rod+m1))
            I0.append(1/3*m_rod*L**2 + m1*(washL)**2)
        else:
            rgo.append((m_rod*(L/2)+m2*washL)/(m_rod+m2))
            I0.append(1/3*m_rod*L**2 + m2*(washL)**2)

#%%
system1 = pd.read_csv("system1_trial2.csv", header = None)
system2 = pd.read_csv("system2_trial1.csv", header = None)
system3 = pd.read_csv("system3_trial1.csv", header = None)
system4 = pd.read_csv("system4_trial2.csv", header = None) 
systems = [system1, system2, system3, system4] #array of each system

errors = []

for i, system in enumerate(systems):
    time = np.array(system.iloc[:, 0])
    disp = np.array(system.iloc[:, 1]*0.0254)
    if i == 2:
        disp = disp[271:550]
        time = time[271:550]
    angle = np.arcsin(disp/L) #converting displacement to angle
    angle = angle - np.average(angle) #mean centering the data, making the mean at 0
    peaks = sig.find_peaks(angle) #finding indices for all peaks
    peak_time = []
    peak_angle = []
    
    for j in peaks[0]: #finding time and angle at each peak
        pt = time[j]
        pa = angle[j]
        peak_time.append(pt)
        peak_angle.append(pa)
        
    expo = curve_fit(lambda t, a, d: a*np.exp(d*t), peak_time, peak_angle) 
    #creating exponential fit to peaks
    coeff = expo[0][1] #coefficient
    curve = expo[0][0] *np.exp(np.multiply(coeff, peak_time)) #creating curve
    
    c = -2*I0[i]*coeff #calculation of damping coefficent
    print(f"Damping Coefficient for System {i+1}: ", c)
       
    sum_peaks = 0
    n = 0
    for w in range(len(peak_time)-1):
        bt_peaks = peak_time[w+1]- peak_time[w] #finding time between peaks
        sum_peaks += bt_peaks #summing these together
        n +=1 #total number of peaks

    periods = sum_peaks / n #average period
    damped = 2*np.pi / periods #damped frequency calculation
    print(f"Damped Frequency for System {i+1}: ", damped)
    
    z = np.sqrt((coeff**2)/((4*(damped**2))+(coeff**2))) #finding z coefficient
    natural = -coeff / (2*z) #natural frequency
    print(f"Experimental Natural Frequency for System {i+1}: ", natural)
    
    #Analytical Solutions
    if i==2:
        ana_wd = np.sqrt((total_m[i]*g*rgo[i] + 2*k*L_B**2)/I0[i])
    else:
        ana_wd = np.sqrt((total_m[i]*g*rgo[i]) / I0[i])
    print(f"Analytically Determined Natural Frequency for System {i+1}: ", ana_wd)
    
    #Percent Error
    error = ((ana_wd - natural) / ana_wd)*100
    print(f"Percent Error of Natural Frequency for System {i+1}: {round(error, 5)}%")
    print("")
    
    #plotting
    plt.scatter(time, angle, marker = "o", s = 2, label = "raw data")
    plt.scatter(peak_time, peak_angle, s = 40, color = "red", label = "peak locations")
    plt.plot(peak_time, curve, color = "black", label = "Exponential Envelope")
    plt.plot(peak_time, -curve, color = "black")
    plt.grid()
    plt.xlabel("Time (seconds)")
    plt.ylabel("Angle (radians)")
    plt.legend(loc = "best")
    plt.savefig(f"system{i+1}_plot.jpg")
    plt.show()
    
    errors.append(error)

sys_num = [1, 2, 3, 4]
plt.bar(sys_num, errors)
plt.xlabel("System Number")
plt.ylabel("Percent Error (%)")
plt.savefig("errorplot.jpg")
plt.show()
           