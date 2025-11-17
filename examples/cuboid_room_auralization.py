"""
Creating an Auralization
========================

Simulate the energy decay in a cuboid room using the Finite Volume diffusion
equation.
"""
# %%
import os
import matplotlib.pyplot as plt
from IPython.display import Audio
import numpy as np
from acousticDE.Auralization.Auralization import run_auralization

script_dir = os.path.dirname(os.path.abspath(__file__))
anechoic_file_path = os.path.join(script_dir, 'anechoic_file.wav') # Full path to the file
results_fvm_path = os.path.join(script_dir, 'resultsFVM.pkl') # Full path to the file

# %%
results = run_auralization(anechoic_file_path,results_fvm_path)

# %%
plt.figure()
plt.plot(
    results['t_conv'],
    results['sh_conv_normalized']/np.max(np.abs(results['sh_conv_normalized'])))
plt.grid(True)
plt.ylabel('Normalized sound pressure (-)')
plt.xlabel('Time (s)')
plt.show()
# %%
Audio(
    results['sh_conv_normalized']/np.max(np.abs(results['sh_conv_normalized'])), 
    rate=results['fs'])
