---
title: '`acousticDE`: A diffusion equation model package for room acoustics simulations'
tags:
  - python
  - diffusion equation
  - room acoustics
  - simulations
  - diffusion process
  - building physics
authors:
  - name: Ilaria Fichera
    orcid: 0000-0002-0097-1486
    #equal-contrib: true
    corresponding: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Cedric Van hoorickx
    orcid: 0000-0002-9671-5558
    #equal-contrib: false # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" 
  - name: Maarten Hornikx
    orcid: 0000-0002-8343-6613
    #equal-contrib: false # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" 
affiliations:
 - name: Department of Built Environment, Eindhoven University of Technology
   index: 1
   ror: 02c2kyt77
date: 24 September 2025
bibliography: paper.bib

---

# Summary
<!-- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience. -->

The acoustic properties of enclosed spaces are important for improving the quality of environments and to reduce noise induced health effects. With the rise of numerical room acoustics modelling [@Kuttruff2019RoomAcoustics], tools are needed that are both computationally efficient, accurate and versatile in capturing sound behaviour in enclosed spaces. `acousticDE` is an open-source computer program using the diffusion equation approach, a modelling method that has gained considerable interest because it offers efficient calcualtions and is flexible with regards to the type of spaces to be addressed. To simulate the sound field in a room, `acousticDE` requires a digital 3D model of the geometry, as well as material properties describing the acoustic boundaries. It can incorporate a sound sources, compute the spatial and temporal distribution of sound energy, and output key room acoustics parameters such as reverberation time. The diffusion equation method, originally introduced by Ollendorff [@Ollendorff1969AbsorptionWaves] and later refined [@Picaut1997AEquation], overcomes some of the limitations of statistical room acoustics, which is grounded in the diffuse field assumption and the classical Sabine reverberation formula [@Sabine1922CollectedAcoustics].
<!-- The Sabine equation also assumes a diffuse field, so this last sentence reads strange (I would omit the sabine mention) -->
Compared to wave-based approaches [@Hamilton2017FDTDTime], [@Wang2019RoomMethod], ray- and beam-tracing algorithms [@Kuttruff2019RoomAcoustics], or the image source method [@Kuttruff2019RoomAcoustics], the DE method provides a computationally lightweight yet physically informed framework for room acoustics simulations, research and practice [@Valeau2006OnPrediction].
<!-- Is this what Valeau write? "computationally lightweight yet physically informed framework for room acoustics simulations" To 'physically informed' could trigger people in thinking about ML methods as phycially infomred ML is a type of method nowadays. Also 'physically informed framework' does not give the reader information on the quality of the DE method. I would try to explain that the method is more accurate than statistical methods and suitable for quick estimates of room acoustic parameters  -->



# Statement of need
<!-- A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work. -->

Room acoustics simulation is essential for research and for consultancies to accurately predict and design the acoustics inside a room. However, existing methods often face a trade off between accuracy, transparency in the calculation method and software and computational efficiency. Traditional statistical methods, such as Sabine and Eyring equations, provide a fast global estimate of the acoustic reverberation time, however these are based on idealized assumptions (diffuse field environment) [@Kuttruff2019RoomAcoustics]. Some more accurate approaches are not commercially available and therefore it is very difficult to see how the calculation is actually done inside the software; this reduces transparency in the method of calculation. 
<!-- Do you mean with the last sentence that more accurate approaches are commercially available (you write: are not commercially available). If you mean that, they are also public, I mean gemoetrical acoustic open source codes.-->
Additionally, full wave-based approaches offer higher accuracy at the cost of significant computational resources, making them impractical for larger spaces. 
<!-- I don't see how this sentence connects to the previous one, as this one does not relate to commercial/open avaialability -->
The diffusion equation model has grown to be an effective compromise, providing physically meaningful spatial and temporal distributions of acoustic energy while remaining computationally efficient in the high-frequency range [@Foy2016IncludingApproach], [@Mou2023AnSpaces], [@SuGul2020ComparativeStations]. However, despite its demonstrated usefulness in the literature, accessible and user-friendly open-source software implementations of diffusion equation model have been limited.
<!-- Do they exist at all, if so we need to refer to that. I think you can refer to our InterNoise 2024 paper on the open source packages in room aocustics, as no DE code was identified back then -->

`acousticDE` is an open source software package designed for the simulation of room acoustics using the diffusion equation model (DE). 
The diffusion equation software has two different numerical methods: the Finite Different Method (FDM) by Du Fort&Frankel [@Navarro2012ImplementationRooms] and the Finite Volume Method (FVM) [@PaganMunoz2019NumericalMethods]. The FDM can be used for parallelepiped shapes and the room is discretised into a 3D grid of points at which the energy density is calculated. The FVM can be used for any shape/geometry, since the discretization is volumetric and tetrahedrical. Beyond numerical simulation, `acousticDE` enables auralization, allowing users to generate perceptual audio realizations of the modelled spaces with its boundary conditions and source location. The Diffusion Equation method software excels for its easeness of use and understanding 
<!-- This last phrase reads like a commercial talk, I mean, many people have difficulties with really translating what the diffusion equation represents in terms of acoustic wave propagation. Do you mean that it is a simple equation that is solved? (Well, the wave equation does not look complicated too ...). I also miss that you generate an impulse response from DE, as that is needed to auralize-->
and for its computational speed, allowing researchers and engineers to simulate quickly a wide range of room acoustics problems with accuracy and computational efficiency. 
<!-- You write 'with accuracy', but that reads like this is the ideal model, while the accuracy cannot be the same of a geometrical acoustics or wave based method, so we need to be a bit more precise here. Somehowe it needs to be clear that we cannot be as accurate as GA or let alone WB methods, but that for some spaces it is very accurate (compact spaces), that improvemets have been made for long spaces and that for other spaces it is approximative-->
The software is written in Python language. In its simpler settings, the `acousticDE` package is available as a back-end acoustic simulator of the Community Hub for Open-source Room Acoustics Software [CHORAS](https://github.com/choras-org/CHORAS) [@Willemsen2025].
<!-- "In its simpler settings", do you mean "With limited settings"? -->

The software has already been applied in a number of scientific publications [@Fichera2024DeterminationApproach], [@Fichera2025] and graduate-level teaching at the Eindhoven University of Technology (TU/e), where it has proven valuable for understanding the impact of room geometry and material absorption on acoustic behavior. By combining ease of use, computational efficiency, and research-grade accuracy, 
<!-- Again, we need to be clear if we write qualifications as "research-grade accuracy", as that leaves room for interpretation for the reader. To me it reads like it is as accurate as a wave-based method -->

`acousticDE` fills a critical need for researchers, consultants, and educators seeking a practical and transparent tool for room acoustics simulation.

# Method

The diffusion equation method aims to find the temporal and spatial distribution of acoustic energy within a specific room.
<!-- It feels like this is the 4rd time that I read a sentence like this now -->
The modelling method is based on solving the partial differential Diffusion Equation, which is applicable in the high-frequency range and it assumes only diffusely reflecting boundaries following the  Lambert's law [@Valeau2006OnPrediction]. The diffusion equation together with its boundary condition is as follow:

\begin{equation}
    \begin{cases}
      \frac{\partial{w(\mathbf{r}, t)}} {\partial{t}} = - D \mathbf{\nabla}^2 w(\mathbf{r}, t) - m c w(\mathbf{r}, t) + q(\mathbf{r}, t) \\
      -D \dfrac{\partial w}{\partial n} = h w \\
    \end{cases}  
    \label{eq:diffusionequationwithRobin}
\end{equation}

where $w(\mathbf{r}, t)$ is the energy density at each position $\mathbf{r}$ and at time $t$, $m$ is the atmospheric attenuation coefficient, $c$ is the speed of sound, $q(\mathbf{r}, t)$ is the source term and $h$ is the absorption term. Several methodologies propose expressions for $h$ based on different assumptions on how to treat the suface absorption coefficients of the room [@Picaut2002NumericalProcess], [@Jing2007APredication], [@Jing2008OnExperiments]. The diffusion coefficient $D$ is a term in $m^2/s$ that scales the energy density and indicates how quickly the sound diffuses around the room. This depends on the mean free path of the room.  

This method allows for an efficient way to calculate acoustics parameters in the room (e.g. sound pressure level and reverberation time) and for more accurate solutions of these parameters compared to the statistical methods of Sabine and Eyring, with the advantage of spatial and temporal distributions of the energies in the rooms while keeping the computational cost reasonable. 
<!-- It feels like this is the 4rd time that I read a sentence like this now -->
<!-- What about the spatially dependent D? -->

# Acknowledgements and Fundings

Ilaria Fichera acknowledges contributions from Silvin Willemsen, Marco Berzborn, Hassan Teymoori, Felipe Raymann, and Radovan Bast for their useful tips. 

This research is founded by the Dutch Research Council (NWO), Applied and Engineering Sciences (AES) under grant agreement No. 19430, with project title 'A new era of room acoustics simulation software: from academic advances to a sustainable open-source project and community'.

# References
<!-- A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline. -->
