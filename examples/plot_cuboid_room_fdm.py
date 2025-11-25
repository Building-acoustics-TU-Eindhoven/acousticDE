"""
The Finite Difference Method
============================

Simulate the energy decay in a cuboid room using the Finite difference diffusion
equation.

As per text below, the inputs need to be prepared. For that, follow the instructions in
`Finite Different Method Use Documentation — Inputs <https://building-acoustics-tu-eindhoven.github.io/acousticDE/Finite%20Difference%20Method%20Use.html#inputs>`_.
"""
# %%
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from acousticDE.FiniteDifferenceMethod.FDM import run_fdm_sim
import tempfile

# %%
# Use a temporary directory for the example
temp_dir = tempfile.TemporaryDirectory()
# You can replace the temporary directory with a specific path if desired
script_dir = temp_dir.name

#%%
###############################################################################
#General input variables
###############################################################################
input_data = {
    "room_dim": [3.0, 3.0, 3.0],
    "coord_source": [1.5, 1.5, 1.5], #source coordinates x,y,z
    "coord_rec": [2.0, 1.5, 1.5], #rec coordinates x,y,z
    "alpha_1": [0.10, 0.15, 0.20, 0.25, 0.25, 0.30], #Absorption coefficient for Surface1 - Floor
    "alpha_2": [0.07, 0.10, 0.13, 0.15, 0.15, 0.16], #Absorption coefficient for Surface2 - Ceiling
    "alpha_3": [0.08, 0.09, 0.11, 0.15, 0.14, 0.14], #Absorption coefficient for Surface3 - Wall Front
    "alpha_4": [0.08, 0.09, 0.11, 0.15, 0.14, 0.14], #Absorption coefficient for Surface4 - Wall Back
    "alpha_5": [0.08, 0.09, 0.11, 0.15, 0.14, 0.14], #Absorption coefficient for Surface5 - Wall Left
    "alpha_6": [0.08, 0.09, 0.11, 0.15, 0.14, 0.14], #Absorption coefficient for Surface6 - Wall Right
    "fc_low": 125, #lowest frequency
    "fc_high": 4000, #highest frequency
    "num_octave": 1, # 1 or 3 depending on how many octave you want
    "dx": 0.5,
    "dt": 1/8000, #time discretization
    "m_atm": 0, #air absorption coefficient [1/m]
    "th": 3, #int(input("Enter type Absortion conditions (option 1,2,3):")) # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
    "tcalc": "decay" #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
}

#%%
###############################################################################
#Creation of json
###############################################################################
fname_input_configuration = "cube_input_fdm.json"
with open(os.path.join(script_dir, fname_input_configuration), "w") as f:
    json.dump(input_data, f, indent=4)

print("Input file successfully created: cube_input_fdm.json")

#%%
###############################################################################
#Run simulation
###############################################################################
result = run_fdm_sim(os.path.join(script_dir, fname_input_configuration))

print("Reverberation time T30 band values:", result["t30_band"])
print("Early decay time EDT band values:", result["edt_band"])
print("Clarity C80 band values:", result["c80_band"])
print("Definition D50 band values:", result["d50_band"])
print("Centre time Ts band values:", result["ts_band"])

#%%
###############################################################################
#Plotting Energy decay curve
###############################################################################
times = result['t'][:len(result['t'])//2]
energy_decay_curve = np.array(result['w_rec_off_band'])
# %%

plt.plot(
    times, 
    10*np.log10(np.abs(energy_decay_curve.T/np.max(energy_decay_curve))), 
    label=[f'{int(fc)} Hz' for fc in result['center_freq']])
plt.grid(True)
plt.ylim([-65, 5])
plt.ylabel('Energy Decay Curve (dB)')
plt.xlabel('Time (s)')
plt.legend()

# %%
temp_dir.cleanup()
