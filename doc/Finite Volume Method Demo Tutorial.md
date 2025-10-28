# Finite Volume Method Demo Tutorial

After installing the package and reading the documentation in [Finite Volume Method Use Documentation](https://building-acoustics-tu-eindhoven.github.io/acousticDE/Finite%20Volume%20Method%20Use.html), the software can be tested with the following files:
- _cube.skp_
- _cube.geo_
- _cube.msh_

The following files can be found in the example folder.

The file _cube.skp_ is a SketchUp file with a room of volume 3x3x3 $m^3$. The files _cube.geo_ and _cube.msh_ are both important files created and saved according to the instructions in [Finite Volume Method Use Documentation, Paragraph 'Geometry & Mesh'](https://building-acoustics-tu-eindhoven.github.io/acousticDE/Finite%20Volume%20Method%20Use.html#geometry-mesh). 

The inputs need to be prepared. For that, follow the instruction in [Finite Volume Method Use Documentation, Paragraph 'General Inputs'](https://building-acoustics-tu-eindhoven.github.io/acousticDE/Finite%20Volume%20Method%20Use.html#general-inputs). This will create a csv file that you need to fill with the absorption coefficients.

Once all the files are created (msh file, json file and csv file), the main acoustics simulation can be run as below:

```
results = run_fvm_sim('C:\....\cube.msh', 'C:\....\cube_input_fvm.json', 'C:\....\absorption_coefficients.csv')
```

The software should provide results of Sound Pressure Levels at the receiver position, Reverberation time, Clarity and other energetic parameters in a pickle file called _resultsFVM.pkl_.

If using the following input data,

```
input_data = {
    "coord_source": [1.5, 1.5, 1.5], #source coordinates x,y,z
    "coord_rec": [2.0, 1.5, 1.5], #rec coordinates x,y,z
    "fc_low": 125, #lowest frequency
    "fc_high": 2000, #highest frequency
    "num_octave": 1, # 1 or 3 depending on how many octave you want
    "dt": 1/20000, #time discretization
    "m_atm": 0, #air absorption coefficient [1/m]
    "th": 3, #int(input("Enter type Absortion conditions (option 1,2,3):")) # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
    "tcalc": "decay" #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
}
```

and absorption coefficients in the .csv file:

```
"absorption coefficient": [0.3, 0.33, 0.5, 0.53, 0.7], #Absorption coefficient for all surfaces
```

test if the software provides the following results:

- Reverberation time (RT): [0.24, 0.22, 0.13, 0.12, 0.09] s; 
- Early Decay Time (EDT): [0.24, 0.22, 0.13, 0.13, 0.09] s;
- Clarity ($C_{80}$): [19.96, 22.25, 36.22, 38.89, 55.43] dB;
- Definition ($D_{50}$): [94.39, 95.95, 99.45, 99.63, 99.96] %
- Centre Time ($T_{s}$): [22.58, 20.75, 14.39, 13.67, 10.66] ms

The result file is a pickle file called _resultsFVM.pkl_. All the results are included in this file.