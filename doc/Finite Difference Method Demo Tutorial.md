# Finite Difference Method Demo Tutorial

After installing the package and reading the documentation in [Finite Difference Method Use Documentation](https://building-acoustics-tu-eindhoven.github.io/acousticDE/Finite%20Difference%20Method%20Use.html), the software can be tested.

The inputs need to be prepared. For that, follow the instruction in [Finite Difference Method Use Documentation, Paragraph 'Inputs'](https://building-acoustics-tu-eindhoven.github.io/acousticDE/Finite%20Difference%20Method%20Use.html#inputs).

Once the files is created (json file), the main acoustics simulation can be run in cmd window with the following commands:

```
python
>>> import acousticDE
>>> from acousticDE.FiniteDifferenceMethod.FVM import run_fdm_sim
>>> results = run_fdm_sim('C:\....\input_fvm.json')
```

The software should provide results of Sound Pressure Levels at the receiver position, Reverberation time, Clarity and other energetic parameters in a pickle file called _resultsFDM.pkl_.


If using the following input data (Navarro et al., 2012):

```
input_data = {
    "room_dim": [8.0, 8.0, 8.0], #dimension of the room x,y,z
    "coord_source": [4.0, 4.0, 4.0], #source coordinates x,y,z
    "coord_rec": [2.0, 2.0, 2.0], #rec coordinates x,y,z
    "alpha_1": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface1 - Floor
    "alpha_2": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface2 - Ceiling
    "alpha_3": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface3 - Wall Front
    "alpha_4": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface4 - Wall Back
    "alpha_5": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface5 - Wall Left
    "alpha_6": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface6 - Wall Right
    "fc_low": 125, #lowest frequency
    "fc_high": 4000, #highest frequency
    "num_octave": 1, # 1 or 3 depending on how many octave you want
    "dx": 0.5,
    "dt": 1/8000, #time discretization
    "m_atm": 0, #air absorption coefficient [1/m]
    "th": 3, #int(input("Enter type Absortion conditions (option 1,2,3):")) # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
    "tcalc": "decay" #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
}
```

test if the software provides the following results:

- Reverberation time (RT): [1.22, 1.22, 1.22, 1.22, 1.22] s; 
- Early Decay Time (EDT): [1.22, 1.22, 1.22, 1.22, 1.22] s;
- Clarity ($C_{80}$): [1.66, 1.66, 1.66, 1.66, 1.66] dB;
- Definition ($D_{50}$): [43.13, 43.13, 43.13, 43.13, 43.13] %
- Centre Time ($T_{s}$): [88.58, 88.58, 88.58, 88.58, 88.58] ms

The result file is a pickle file called _resultsFDM.pkl_. All the results are included in this file.

## References
- J. M. Navarro, J. Escolano, J. J. Lopez, Implementation and evaluation of a diffusion equation model based on finite difference schemes for sound field prediction in rooms, Applied Acoustics 73 (6-7) (2012).

