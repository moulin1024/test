import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import windfarm as wf
# Author: Mou Lin 

# TODO: Implement a GUI for parameter inputs

Lx = 1.28*4 
Ly = 1.28
Nx = int(64*4)
Ny = int(64)
inflow = (4.88,0.07)
wt_location = [(0.16,0.64)]
yaw_angle = [0.0]

# Instantiate a wind farm with the given parameter
wf1 = wf.WindFarm(wt_location,(Lx,Ly),(Nx,Ny),inflow)
# Superpose the wake generate by each wind turbine
wf1.wake_superposition()
# Output the flow field (streamwise velocity and turbulence intensity)
wf1.flowfield_contour() 

# ========== Post processing ===========
x_coord = np.linspace(0,Lx,Nx)
y_coord = np.linspace(0,Ly,Ny)
xv, yv = np.meshgrid(x_coord, y_coord)

# Read the outputed flow field 
velocity = wf1.velocity_field.value
turb_intensity = wf1.turb_intensity_field.value

# Plot the contours
fig = plt.figure()
cs = plt.contourf(xv, yv, velocity,10)
plt.axis('scaled')
cbar = fig.colorbar(cs)
plt.show()

fig = plt.figure()
cs = plt.contourf(xv, yv,turb_intensity,10)
plt.axis('scaled')
cbar = fig.colorbar(cs)
plt.show()