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
  - name: Cédric Van hoorickx
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
date: 24 December 2025
bibliography: paper.bib

---

# Summary
<!-- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience. -->

Studying the acoustic properties of enclosed spaces is important for improving the quality of environments and reducing noise induced health effects. With the rise of numerical room acoustics modelling [@Kuttruff2019RoomAcoustics], tools have been developed that are both computationally efficient and versatile in capturing sound behaviour in enclosed spaces. `acousticDE` is an open-source software for simulating room acoustics using the diffusion equation approach. The software is designed to be easy to use and suitable not only for research applications but also for educational purposes. In practice, it can be used for evaluating design choices during early architectural design stages, optimising acoustic treatment layouts or testing performance of room intended for speech, music or general occupancy.

The method implemented in `acousticDE` - the diffusion equation method - has gained considerable interest because it offers efficient calculations and is flexible with regards to the type of spaces to be addressed. To simulate the sound field in a room, `acousticDE` requires a digital 3D model of the geometry, as well as material properties describing the acoustic boundaries. Once configured, `acousticDE` can incorporate sound sources, compute the distribution of sound energy throughout the space, and provide key room acoustics parameters relevant to practical acoustic assessments. 


# Statement of need
<!-- A section that clearly illustrates the research purpose of the software and places it in the context of related work. This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work. -->

Room acoustics simulation is essential for researchers and for consultants to accurately predict and design the acoustics inside a room. However, existing methods often face a trade off between accuracy, transparency in the calculation method and software and computational efficiency. Traditional statistical methods, such as the Sabine and Eyring calculations, provide a fast global estimate of the acoustic parameters, but these are based on unrealistic assumptions that could be valid for certain rooms and for certain frequency ranges (e.g. diffuse field environment) [@Kuttruff2019RoomAcoustics]. More accurate approaches, including both commercial and open-source software are available. However, commercial tools often lack transparency, as it is difficult to understand, how the calculations are performed internally. These are mostly based on common approaches to run high frequency simultions, i.e. with a geometrical acoustics method. Additionally, full wave-based approaches can achieve even higher accuracy but require significant computational resources, making them impractical for larger spaces. This creates a gap between simple tools, highly detailed numerical simulations and open source software, particularly for calculation of spatial information inside a room with manageable computational costs.

`acousticDE` was developed to address this gap. It is an open-source software package designed for the simulation of room acoustics using the diffusion equation model (DE), targeting researchers, consultants and educators who require a transparent and efficient method to high-frequency room acoustics.


# State of the field
<!-- A description of how this software compares to other commonly-used packages in the research area. If related tools exist, provide a clear “build vs. contribute” justification explaining your unique scholarly contribution and why existing alternatives are insufficient. -->

The diffusion equation model has grown to be an effective compromise for room acoustics simulations in the high frequency range, providing physically meaningful spatial and temporal distributions of acoustic energy while remaining computationally efficient [@Foy2016IncludingApproach], [@Mou2023AnSpaces], [@SuGul2020ComparativeStations]. However, despite its demonstrated usefulness in the literature, accessible and user-friendly software implementations of the diffusion equation model have been non existent [@Hornikx2024] or are somewhat limited. In fact, most implementations of the diffusion equation embedded in commercial tools (e.g. COMSOL or Femlab toolbox for Matlab), where room acoustics simulation are not the primary focus. An additional open source implementation exists, but it is currently at its initial version [@univgustaveeiffel]. 

`acousticDE` distinguishs itself by providing two different numerical methods for implementing the diffusion equation: the Du Fort&Frankel Finite Difference Method (FDM) [@Navarro2012ImplementationRooms] and the Finite Volume Method (FVM) [@PaganMunoz2019NumericalMethods]. The FDM can be used for parallelepiped shapes and the room is discretised into a 3D grid of points where the energy density is calculated. The FVM can be used for any 3D geometry, since the discretization is with tetrahedrical elements. The dual approach allows the user to select the appropriate balance between flexibility and efficiency. Compared to other methods for the high frequency range, the DE model is easy to use and understand and has a low computational speed, allowing users to quickly simulate room acoustic problems. However, for complex geometries or highly absorptive spaces, the diffusion equation does not reach the level of accuracy of more detailed (but computationally more costly) acoustics methods (geometrical acoustics or wave-based method). It is instead accurate for compact spaces with high reflective surfaces. Beyond numerical simulation, `acousticDE` enables the creation of the impulse response from the numerical results of the diffusion equation and auralization, allowing users to generate perceptual audio demonstrations of acoustic conditions. 

# Software design
<!-- An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application. This should demonstrate meaningful design thinking beyond a superficial code structure description. -->

The structure architecture adapts to the numerical method selected by the user. In this section, the FVM and Auralization are used as a reference case, as it represents the most structurally demanding workflow. The overall design follows a sequential pipeline:
1. creation of the 3D geometry;
2. conversion of the geometry into mesh;
3. specification of inputs paramters;
4. execution of the acoustics diffusion equation simulator;
5. execution of the auralization module.

This pipeline was intentionally chosen to make each stage of the modelling process explicit and inspectable. Room acoustics simulators are highly sensitive to early-stage modeling errors, such as gaps in the geometry, incorrect surface definitions or inconsistent material properties. By structuring the software as a sequence of separated steps, the user is encouraged to verify the geometry, inspect the mesh and validate the physical inputs before committing to the simulation. This practice reduces the likelihood of discovering errors only after a simulation has already begun. The core data structure used troughtout the software are linear data structures, with Numpy arrays, matrices and vector. This choice follows the mathematical and numerical formulation of the diffusion equation as discretised fields and time-stepping schemes are naturally represented by dense or sparse arrays. These structures allow efficient numerical operations. Mesh information imported from Gmsh is also stored as arrays, maintaining the consistency. 

The main simulation is organised in a set of modular, sequential functions. Each function has its own task and returns data that are needed by the next function in the next stage of the simulation. This functional decomposition improves readibility and support maintainability. In some cases, Python dictionaries appears - e.g. grouping material properties - but do not affect the mainly linear data structure and neither the computational efficiency of the software.

Python is chosen as implementation language because it is open-source, widely accessible, compared to alternatives such as Matlab or C++, which are not free and more dificult to implement and to understand from a user perspective.


# Method

The diffusion equation method, originally introduced by Ollendorff [@Ollendorff1969AbsorptionWaves] and later refined [@Picaut1997AEquation], overcomes some of the limitations of statistical room acoustics, which is grounded in the diffuse field assumption [@Sabine1922CollectedAcoustics]. Compared to wave-based approaches [@Hamilton2017FDTDTime], [@Wang2019RoomMethod], ray- and beam-tracing algorithms [@Kuttruff2019RoomAcoustics], or the image source method [@Kuttruff2019RoomAcoustics], the DE method provides a computationally lightweight software for room acoustics simulations, research and practice. The method is more accurate than statistical methods and suitable for quick estimates of room acoustic parameters.

The diffusion equation method aims to estimate how acoustic energy is distributed over time and space within a specific room. The modelling method is based on solving the partial differential diffusion equation, which is applicable in the high-frequency range and it assumes only diffusely reflecting boundaries following Lambert's law [@Valeau2006OnPrediction]. The diffusion equation method is based on statistical and energetical approach, where sound has undergone multiple reflections and scattering in an environment. It represents the late stage of wave propagation, where the energy, after being reflected multiple times, is spread out across the space, following the diffusion process of gas and heat in a medium. The method has some limitations, e.g. diffraction is not taken into account, the early decay part of the decay is not included, surfaces are fully scattering and it is not accurate for average high absorption coefficients. The method includes, though, the basic and most physical phenomena for room acoustics: surface absorption, atmospheric absorption and diffusion.

The diffusion equation together with its boundary condition reads as follows:

$$
\frac{\partial{w(\mathbf{r}, t)}} {\partial{t}} = - D \boldsymbol{\nabla}^2 w(\mathbf{r}, t) - m c w(\mathbf{r}, t) + q(\mathbf{r}, t) \\
$$

$$
-D \dfrac{\partial w}{\partial n} = h w
$$

<!-- \begin{equation}
    \begin{cases}
      \frac{\partial{w(\mathbf{r}, t)}} {\partial{t}} = - D \mathbf{\nabla}^2 w(\mathbf{r}, t) - m c w(\mathbf{r}, t) + q(\mathbf{r}, t) \\
      -D \dfrac{\partial w}{\partial n} = h w \\
    \end{cases}  
\end{equation} -->

where $w(\mathbf{r}, t)$ is the energy density at each position $\mathbf{r}$ and at time $t$, $m$ is the atmospheric attenuation coefficient, $c$ is the speed of sound, $q(\mathbf{r}, t)$ is the source term and $h$ is the boundary absorption term. Several methodologies propose expressions for $h$ based on different assumptions on how to treat the suface absorption coefficients of the room [@Picaut2002NumericalProcess], [@Jing2007APredication], [@Jing2008OnExperiments]. The diffusion coefficient $D$ is a term in m$^2$/s that indicates how quickly the sound diffuses around the room. This depends on the mean free path of the room.

From the energy density, the reverberant field Sound Pressure Level (SPL) can be calculated as follows: 
$$
\text{SPL}(\mathbf{r},t) = 10  \log_{10} \left(\frac{\rho c^2 w(\mathbf{r},t)} {p_{\text{ref}}^2}\right) 
$$
where $\rho$ is the density of air in [kg/m$^3$] and $p_{\text{ref}}$ is equal to $2 \times 10^ {- 5}$ Pa [@Navarro2015].

This method allows for an efficient way to calculate acoustics parameters in the room (e.g. sound pressure level and reverberation time) and for more accurate solutions of these parameters compared to the statistical methods of Sabine and Eyring, while keeping the computational cost reasonable.


# Research impact assessment
<!-- Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals). The evidence should be compelling and specific, not aspirational. -->

The software can be used as a standalone simulation tool or as a back-end acoustic simulator of the Community Hub for Open-source Room Acoustics Software [CHORAS](https://github.com/choras-org/CHORAS) [@Willemsen2025]. The software has already been applied in a number of scientific publications [@Fichera2024DeterminationApproach], [@Fichera2025], [@Fichera2026AnRooms] and graduate-level teaching at the Eindhoven University of Technology (TU/e), where it has proven valuable for understanding the impact of room geometry and material absorption on acoustic behavior and for designing rooms from an acoustics point of view. By combining ease of use and computational efficiency, `acousticDE` fills a critical need for researchers, consultants, and educators seeking a practical, transparent and approximate tool for room acoustics simulation.

# AI usage disclosure
<!-- Transparent disclosure of any use of generative AI in the software creation, documentation, or paper authoring. If no AI tools were used, state this explicitly. If AI tools were used, describe how they were used and how the quality and correctness of AI-generated content was verified.--> 

No AI was used for the creation of `acousticDE`.

# Acknowledgements and Fundings

Ilaria Fichera acknowledges contributions from Silvin Willemsen, Marco Berzborn, Hassan Teymoori, Amin Livani, Felipe Raymann, and Radovan Bast for their useful tips. 

This research is funded by the Dutch Research Council (NWO), Applied and Engineering Sciences (AES) under grant agreement No. 19430, with project title 'A new era of room acoustics simulation software: from academic advances to a sustainable open-source project and community'.

# References
<!-- A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline. -->
