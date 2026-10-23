"""
Creating an Auralization
========================

Simulate the energy decay in a cuboid room using the Finite Volume diffusion
equation.
"""
# %%
import tempfile
import matplotlib.pyplot as plt
from IPython.display import Audio
import numpy as np
import pooch
from acousticDE.Auralization.Auralization import run_auralization

tmp_dir = tempfile.TemporaryDirectory()
script_dir = tmp_dir.name
#%%
###############################################################################
#Import anechoic file
###############################################################################
anechoic_file_path = pooch.retrieve(
    url="https://github.com/Building-acoustics-TU-Eindhoven/acousticDE/raw/refs/heads/master/examples/anechoic_file.wav",
    known_hash=None,
    path=script_dir,
    fname="anechoic_file.wav"
)

#%%
###############################################################################
#Import FVM or FDM results
###############################################################################

results_fvm_path = pooch.retrieve(
    url="https://github.com/Building-acoustics-TU-Eindhoven/acousticDE/raw/eba9989df48768f8c2574b4fd50c5d4d75bf0213/examples/resultsFVM.pkl",
    known_hash=None,
    path=script_dir,
    fname="resultsFVM.pkl"
)

#%%
###############################################################################
#Run script
###############################################################################
results = run_auralization(anechoic_file_path,results_fvm_path)

#%%
###############################################################################
#Plotting convolved signal
###############################################################################
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
