# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 16:12:40 2023

@author: Ilaria Fichera
"""

#%%
###############################################################################
#IMPORT LIBRARIES
###############################################################################
#Code developed by Ilaria Fichera for the analysis of the FVM method adapted solving the 3D diffusion equation with one intermittent omnidirectional sound source
#Import modules
import math
import pickle
import time
import types
from math import ceil
from math import log
import gmsh
import numpy as np
import json

from acousticDE.FiniteVolumeMethod.FunctionClarity import clarity
from acousticDE.FiniteVolumeMethod.FunctionDefinition import definition
from acousticDE.FiniteVolumeMethod.FunctionCentreTime import centretime
from acousticDE.FiniteVolumeMethod.FunctionRT import t60_decay

import logging

# Create logger for this module
logger = logging.getLogger(__name__)

def check_should_cancel(json_file_path_in):
    try:
        if json_file_path_in is not None:
            with open(json_file_path_in, "r") as json_file_to_check:
                data = json.load(json_file_to_check)

        # Update the specified field value
        if "should_cancel" in data:
            return data["should_cancel"]

    except Exception as e:
        print("check_should_cancel returned: " + str(e))

#%%
###############################################################################
#SURFACE MATERIALS FUNCTIONS
###############################################################################

def number_freq(num_octave,fc_high,fc_low):
    """
    Calculate the frequency array and the number of frequency bands.

    Parameters
    ----------
    num_octave : int
        Number of octaves to calculate (1 or 3 octaves).
    fc_high : int
        The highest frequency in the calculation.
    fc_low : int
        The lowest frequency in the calculation.

    Returns
    -------
    nBands : int
        Number of frequency bands.
    center_freq : list of float
        Array of all the frequencies to calculate.
    """
    x_frequencies  = num_octave * log(fc_high/fc_low) / log(2)
    nBands = int(num_octave * log(fc_high/fc_low) / log(2) + 1)
    center_freq = fc_low * np.power(2,((np.arange(0,x_frequencies+1) / num_octave)))
    return nBands, center_freq


# Absorption term for boundary conditions 
def abs_term(th,abscoeff_list,c0):
    """
    Calculate the absorption term (Sabine, Eyring or Modified)

    Parameters
    ----------
        th : int
            The options for the absorption term; Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
        abscoeff_list : list
            Absrption coefficient for each frequency
        c0 : int 
            Speed of sound 

    Returns
    -------
        Absx_array : array of floats
            Calculated absorption term for each absorption coefficient for each frequency
    """
    Absx_array = np.array([])
    for abs_coeff in abscoeff_list:
        #print(abs_coeff)
        if th == 1:
            Absx = (c0*abs_coeff)/4 #Sabine
        elif th == 2:
            Absx = (c0*(-log(1-abs_coeff)))/4 #Eyring
        elif th == 3:
            Absx = (c0*abs_coeff)/(2*(2-abs_coeff)) #Modified by Xiang
        Absx_array = np.append(Absx_array, Absx)
    return Absx_array

def create_vgroups_names(file_path, should_initialise_gmsh=True):
    """
    Create a list of the material names assigned in SketchUp

    Parameters
    ----------
        file_path : str
            Full path to the mesh file

    Returns
    -------
        vGroupsNames : list
            Names of the materials in the msh file (the material name are the same as the one assigned in the SketchUp file)
    """
    if should_initialise_gmsh:
        gmsh.initialize() #Initialize msh file
    mesh = gmsh.open(file_path) #open the file
    dim = -1 #dimensions of the entities, 0 for points, 1 for curves/edge/lines, 2 for surfaces, 3 for volumes, -1 for all the entities 
    tag = -1 #all the nodes of the room
    vGroups = gmsh.model.getPhysicalGroups(-1) #these are the entity tag and physical groups in the msh file. 
    vGroupsNames = [] #these are the entity tag and physical groups in the msh file + their names
    for iGroup in vGroups:
        dimGroup = iGroup[0]  #entity tag: 1 lines, 2 surfaces, 3 volumes (1D, 2D or 3D)
        tagGroup = iGroup[1]  #physical tag group (depending on material properties defined in SketchUp)
        namGroup = gmsh.model.getPhysicalName(dimGroup, tagGroup) #names of the physical groups defined in SketchUp   
        alist = [dimGroup,tagGroup,namGroup] #creates a list of the entity tag, physical tag group and name
        #print(alist)
        vGroupsNames.append(alist)

    return vGroupsNames

# Gives absorption coefficients to a material (group) and links it to the surfaces
def surface_materials(group, abscoeff, surface_absorption, absorption_coefficient_dict, nBands,th,c0):
    """
    Calculation of absorption term for each material

    Parameters
    ----------
        group : list
            List of element type (1 for lines, 2 for surface, 3 for volumes), element type number and material name in the msh file
        abscoeff : list of strings
            Absorption coefficient of the group (material name)
        surface_absorption : list of tuples
            List initialization for each frequency of a tuple including the surface number and the absorption term for that surface
        absorption_coefficient_dict : dict
            Dictionary initialization of absorption coefficients per each surface and per each frequency 
        nBands : int
            Number of frequency bands
        th : int
            Option for the absorption term; Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
        c0 : int
            Speed of sound 

    Raises
    ------
        Exception: If the number of absorption coefficient typed are higher than the number of frequency bands it will raise an error

    Returns
    -------
        absorption_coefficient_dict : dict
            absorption coefficients per each surface and per each frequency
        surface_absorption : list of tuples
            Absorption term for each surface and for each frequency
    """

    #abscoeff = abscoeff.split(",")
    if len(abscoeff) != nBands:
        logger.error("Number of absorption coefficients doesn't match the number of frequency bands")
        raise Exception("Number of absorption coefficients doesn't match the number of frequency bands")

    #abscoeff = [float(i) for i in abscoeff][-1] #for one frequency
    abscoeff_list = [float(i) for i in abscoeff] #for multiple frequencies
    
    physical_tag = group[1] #Get the physical group tag
    entities = gmsh.model.getEntitiesForPhysicalGroup(2, physical_tag) #Retrieve all the entities in this physical group (the entities are the number of walls in the physical group)

    Abs_term = abs_term(th, abscoeff_list,c0) #calculates the absorption term based on the type of boundary condition th
    for entity in entities:
        absorption_coefficient_dict[entity] = abscoeff_list
        surface_absorption.append((entity, Abs_term.copy())) #absorption term (alpha*surfaceofwall) for each wall of the room
        surface_absorption = sorted(surface_absorption, key=lambda x: x[0])
    
    return absorption_coefficient_dict, surface_absorption

#%%
###############################################################################
#GMSH GET NODES, VOLUME ELEMENTS AND BOUNDARY ELEMENTS
###############################################################################
def get_nodes_elem(dim,tag):
    """
    Calculates the number of volumetric and boundary elements in the mesh

    Parameters
    ----------
        dim : int 
            Dimensions of the entities, 0 for points, 1 for curves/edge/lines, 2 for surfaces, 3 for volumes, -1 for all the entities 
        tag : int
            Indication for nodes (-1 indicates all the nodes of the room) 

    Returns
    -------
        nodecoords : array of floats
            The coordinates of each node in the mesh
        node_indices : dict 
            Indices of all the nodes in the mesh
        bounEl : array of int 
            Indices of all the boundary surfaces in the mesh
        bounNode : array of int 
            Indices of all the boundary nodes per each boundary surface in the mesh
        voluEl : array of int
            Indices of all the volume elements (tetrahedra) in the mesh
        voluNode : array of int 
            Indices of all the volumetric nodes per each colume element in the mesh
        belemNodes : array of int
            Indices of all the boundary nodes per each boundary surface in the mesh
        velemNodes : array of int 
            Indices of all the volumetric nodes per each colume element in the mesh
        boundaryEl_dict : dict
            Dictionary with key the index of the boundary element (boudnary surface) and value the Indices of the nodes of the surface
        volumeEl_dict : dict 
            Dictionary with key the index of the volumetric element (tetrahedra) and value the Indices of the nodes of the tetrahedra
    """
    # Nodes
    #tag = -1  # all the nodes of the room
    nodeTags, coords, parametricCoord = gmsh.model.mesh.getNodes(dim,
                                                                    tag)  # gets the tags for each node and the coordinates of each node
    nodecoords = coords.reshape((-1, 3))  # coordinates reshaped in a matrix 3xnumber of nodes

    node_indices = {tag: index for (index, tag) in enumerate(nodeTags)}

    # Element Types
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
    # elemTypes = 1 for lines, 2 for surfaces, 4 for tetrahedron
    # elemTags =  list of list of lines, boundary elements (surfaces) and volume elements (tetrahedron)
    # elemNodeTags = did not understand this yet, probably a tag gives to each line, surface and tetrahedron
    for e_type in elemTypes:
        if e_type == 1:  # if the e_type = 1, then get all the elements (lines) with that e_type
            tag = -1
            edgeEl, edgeNodeTagstype = gmsh.model.mesh.getElementsByType(e_type,
                                                                            tag)  # get the edge element tags and the nodes; edgeEl are lines
            # edgeEl = numbered edge elements (lines)
        elif e_type == 2:  # if the e_type = 2, then get all the elements (surfaces) with that e_type
            tag = -1
            bounEl, bounNodeTagstype = gmsh.model.mesh.getElementsByType(e_type,
                                                                            tag)  # get the boundary/surface element tags and the nodes; bounEl are surfaces
            bounNode = bounNodeTagstype.reshape((-1, 3))
            # boundEl = numbered boundary elements (surfaces)
            # bounNode = nodes of the surfaces
        elif e_type == 4:
            tag = -1  # if the e_type = 4, then get all the elements (tetrahedron) with that e_type
            voluEl, voluNodeTagstype = gmsh.model.mesh.getElementsByType(e_type,
                                                                            tag)  # get the volume element tags and the nodes; voluEl are tetrahedrons
            voluNode = voluNodeTagstype.reshape((-1, 4))
            # voluEl = numbered volume elements (tetrahedron)
            # voluNode = nodes of the tetrahedon
        else:
            print("The procedure is not possible")

    velemNodes = elemNodeTags[2].reshape((-1, 4))  # nodes per each tetrahedron
    belemNodes = elemNodeTags[1].reshape((-1, 3))  # nodes per each surface boundary

    eelement = 0
    for elem in range(len(edgeEl)):
        eelement = eelement + 1  # scalar number of the edge elements

    # Volume Element dictionary + nodes of each volume elements (4 nodes per element)
    volumeEl_dict = {}  # Initialization Dictionary of volumelements + its nodes
    for i in range(len(voluEl)):
        # print(i)
        volumeEl_dict[voluEl[i]] = velemNodes[i]  # Dictionary of volumelements + its nodes

    # Boundary Element dictionary + node per each surface elements (3 nodes per element)
    boundaryEl_dict = {}  # Initialization Dictionary of boundary elements + its nodes
    for i in range(len(bounEl)):
        boundaryEl_dict[bounEl[i]] = belemNodes[i]  # Dictionary of boundary elements + its nodes

    return nodecoords, node_indices, bounEl, bounNode, voluEl, voluNode, belemNodes, velemNodes, boundaryEl_dict, volumeEl_dict

#%%
###############################################################################
#CHECKING IF POINT SOURCE AND RECEIVER ARE INSIDE ROOM
###############################################################################
def point_inside_room(coordinate, nodecoords, velemNodes):
    """
    Checks if the point is inside the room volume

    Parameters
    ----------
        coordinate : : list 
            Coordinates of the source or receiver position
        nodecoords : array of floats
            The coordinates of each node in the mesh
        velemNodes : array of int 
            Indices of all the volumetric nodes per each colume element in the mesh

    Returns
    -------
        False : bool
            If the point is outside the room    
    """
    inside = False
    for tet in velemNodes:  # each tetra = indices of its 4 corner nodes
        tetra_coords = nodecoords[tet-1]
        T = np.hstack((tetra_coords, np.ones((4,1))))
        v = np.append(coordinate, 1)
        bary = np.linalg.solve(T.T, v)
        if np.all(bary >= -1e-12) and np.all(bary <= 1 + 1e-12):
            inside = True
            break
    if not inside:
        raise ValueError("Source or receiver points are outside the room mesh.")
        return False


#%%
###############################################################################
#CALCULATION OF VOLUME CELLS AND CENTRE OF VOLUME
###############################################################################

def velem_volume_centre(volumeEl_dict,nodecoords,node_indices):
    """
    Calculation of the volume of each cell and its center

    Parameters
    ----------
        volumeEl_dict : dict 
            Dictionary with key the index of the volumetric element (tetrahedra) and value an array with the Indices of the nodes of the tetrahedra (maximum 4 Indices)
        nodecoords : array of floats 
            The coordinates of each node in the mesh
        node_indices : dict
            Indices of all the nodes in the mesh

    Returns
    -------
        cell_center : array of floats 
            The coordinates of the center of the cell element per each element
        cell_volume : array of floats 
            The volume of each element cell
    """
    #Calculation of volume cells and centre of volume    
    vcell_dict = {} #volume of each element tetrahedron initialization
    centre_cell = {} #centre of the element tetrahedron initialization
    for i in volumeEl_dict.keys():
        coord_centre_cell = np.zeros(3)
        centre_cell[i] = []
        #print(i)
        vc0 = nodecoords[node_indices[volumeEl_dict[i][0]],:] #coordinates of node 0
        #vc0 = gmsh.model.mesh.getNode(volumeEl_dict[i][0])[0] #Coordinates of the node number zero of the volume element i
        #print(nc0)
        vc1 = nodecoords[node_indices[volumeEl_dict[i][1]],:]
        vc2 = nodecoords[node_indices[volumeEl_dict[i][2]],:]
        vc3 = nodecoords[node_indices[volumeEl_dict[i][3]],:]
        #vc1 = gmsh.model.mesh.getNode(volumeEl_dict[i][1])[0] #Coordinates of the node number one of the volume element i
        #print(nc1)
        #vc2 = gmsh.model.mesh.getNode(volumeEl_dict[i][2])[0] #Coordinates of the node number two of the volume element i
        #print(nc2)
        #vc3 = gmsh.model.mesh.getNode(volumeEl_dict[i][3])[0] #Coordinates of the node number three of the volume element i
        #print(nc3)
        for j in range(3): #three coordinates per each node
            coord_centre_cell[j] = (vc0[j]+vc1[j]+vc2[j]+vc3[j])/4 #coordinates of the centre of each volume element
            centre_cell[i].append(coord_centre_cell[j])
        vcell_dict[i] = abs(np.dot(np.cross(vc1-vc3,vc2-vc3),vc0-vc3))/6 #volume of each volume element
    
    #The dictionary centre cell is modified into an array of floats (list)
    cell_center = np.array(()) #initialization of an array of centre cell coordinates from the centre_cell dictionary
    for key in centre_cell:
        cell_center = np.append(cell_center,centre_cell[key])
    cell_center = cell_center.reshape((-1,3))
    
    cell_volume = np.array(()) #initialization of an array of cell volumes from the vcell_dict dictionary
    for key in vcell_dict:
        cell_volume = np.append(cell_volume,vcell_dict[key])
    
    return cell_center, cell_volume

#%%
###############################################################################
#CALCULATION OF BOUNDARY ELEMENTS AND CENTRE OF AREA (might not need this function actually)
###############################################################################
#Calculation of boundary elements area and centre 
def belem_area_centre(boundaryEl_dict,nodecoords,node_indices):
    """
    Calculation of the area of each boundary surface and its center

    Parameters
    ----------
        boundaryEl_dict : dict
            Dictionary with key the index of the boundary element (boudnary surface) and value an array with the Indices of the nodes of the surface (maximum 3 Indices)
        nodecoords : array of floats
            The coordinates of each node in the mesh
        node_indices : dict
            Indices of all the nodes in the mesh

    Returns
    -------
        barea_dict : dict 
            Area of each boundary surface element
        centre_area : dict 
            Center point of each boundary surface element
    """
    barea_dict = {} #surface of each element boundary initialization
    centre_area = {} #centre of the element tetrahedron initialization
    for i in boundaryEl_dict.keys():
        coord_centre_area = np.zeros(3)
        centre_area[i] = []
        #print(i)
        bc0 = nodecoords[node_indices[boundaryEl_dict[i][0]],:]
        
        #bc0 = gmsh.model.mesh.getNode(boundaryEl_dict[i][0])[0]
        #bnodeCoord_dict[boundaryEl_dict[i][0]] #Coordinates of the node number zero of the volume element i
        #print(nc0)
        bc1 = nodecoords[node_indices[boundaryEl_dict[i][1]],:]
        #gmsh.model.mesh.getNode(boundaryEl_dict[i][1])[0] #Coordinates of the node number one of the volume element i
        #print(nc1)
        bc2 = nodecoords[node_indices[boundaryEl_dict[i][2]],:]
        #gmsh.model.mesh.getNode(boundaryEl_dict[i][2])[0] #Coordinates of the node number two of the volume element i
        #print(nc2)
        for j in range(3):
            coord_centre_area[j] = (bc0[j]+bc1[j]+bc2[j])/3 #coordinates of the centre of each volume element
            centre_area[i].append(coord_centre_area[j])
        barea_dict[i] = abs(sum(np.cross(bc2-bc1,bc1-bc0)))/2 #volume of each volume element
    
    return barea_dict, centre_area


#%%
###############################################################################
#CALCULATION OF NEIGHBOURS (might not need this function actually)
###############################################################################
#Neighbours calculation; What are the neighbours faces of each volume? 3 per each minimum?
def get_neighbour_faces(voluEl):
    """
    Calculation of nighbours volumetric elements and faces

    Parameters
    ----------
        voluEl : array of int 
            Indices of all the volume elements (tetrahedra) in the mesh

    Returns
    -------
        fxt : dict 
            Dictionary with keys as the nodes of each face and values the volume elements of which this face is neighbour
        txt : dict 
            Dictionary with keys as the tetrahedron tag and values as the tet that are neighbours (the tetrahedrons at the boundary are not counted)
        neighbourVolume : array of floats 
            Array with the neighbours tetrahedron per each tetrahedron in order from 0 to the number of tetrahedrons
    """
    facenodes = gmsh.model.mesh.getElementFaceNodes(4,
                                                    3)  # 4 is the element type (tetrahedron) and three are the nodes per each face #get all the face tags of all the faces of the tetrahedrons

    # Computing face x tetrahedon incidence
    faces = []  # initialization of list of tuples of the faces nodes
    fxt = {}  # dictionary with keys as the nodes of each face and values the volume elements of which this face is neighbour
    for i in range(0, len(facenodes), 3):  # per each element basically, goes trhough the nodes of each face 3by3
        # print(i)
        f = tuple(sorted(facenodes[i:i + 3]))  # nodes of each face put in a tuple from node i to node i plus 3
        faces.append(f)
        tet = voluEl[i // 12]  # volume element number at which the faces are associated?
        if not f in fxt:  # if the face f (with its node) is already in the dictionary, just append the volume element neighbour to
            fxt[f] = [tet]
        else:
            fxt[f].append(tet)

    # Computing neighbors by face
    txt = {}  # dictionary with keys as the tetrahedron tag and values as the tet that are neighbours (the tetrahedrons at the boundary are not counted)
    for i in range(0, len(faces)):
        # print(i)
        f = faces[i]  # f is a tuple of the nodes of the face into consideration
        tet = voluEl[i // 4]  # tetrahedron at which the face is neighbour
        if not tet in txt:  # if the tet is not in the dictionary, add it, otherwise append the new tt
            txt[tet] = []
        for tt in fxt[f]:
            if tt != tet:
                txt[tet].append(int(tt - (voluEl[0] - 1)))  # volumes neighbours to each volume
    for values in txt.values():
        while len(values) < 4:  # if there are less than 4 tetrahedron neighbours it means that the others are boundary tetrahedrons, therefore add a zero for each tetrahedron missing
            values.append(0)

    neighbourVolume = np.array(
        ())  # initialization of an array with the neighbours tetrahedron per each tetrahedron in order from 0 to the number of tetrahedrons
    for key in txt:
        neighbourVolume = np.append(neighbourVolume, txt[key])
    for item in neighbourVolume:
        item = int(item)

    neighbourVolume = neighbourVolume.reshape(
        (-1,
            4))  # reshape the array so that we have the 4 tetrahedrons neighbours of the tetrahedron in consideration

    return fxt, txt, neighbourVolume

#%%
###############################################################################
#CALCULATION OF INTERIOR TETRAHEDRONS
###############################################################################
# Interior Tetrahedrons calculations

def interior_tetra(voluEl,cell_center,velemNodes,nodecoords,node_indices):
    """
    Calculation of shared area divided shared distance between tetrahedrons

    Parameters
    ----------
        voluEl : array of int 
            Indices of all the volume elements (tetrahedra) in the mesh
        cell_center : array of floats
            The coordinates of the center of the cell element per each element
        velemNodes : array of int 
            Indices of all the volumetric nodes per each colume element in the mesh
        nodecoords : array of floats 
            The coordinates of each node in the mesh
        node_indices : dict 
            Indices of all the nodes in the mesh

    Returns
    -------
        interior_tet : array of floats
            Matrix of tetrahedron per tetrahedron of the division between shared area and shared distance
        interior_tet_sum : array of floats
            Sum of interior_tet per columns
    """
    interior_tet = np.zeros((len(voluEl), len(voluEl))) #initialization matrix of tetrahedron per tetrahedron
    
    for i in range(len(voluEl)): #for each tetrahedron, take its centre
        #print(i)
        cell_center_i = cell_center[i]
        for j in range(len(voluEl)): #for each tetrahedron, take its centre
            #cell_center_j = cell_center[j]
            #print(j)
            if i != j: #if the tetrahedrons are not the same one, then check if there are shared nodes in between the two tetrahedron i and j
                shared_nodes = np.intersect1d(velemNodes[i], velemNodes[j])
                #shared_nodes = []
                #count = 0
                #for node in velemNodes[i]: #for each node in tetrahedron i
                #    print(node)
                #    if node in velemNodes[j]: #if each node of the tetrahedron i is in nodelist of tetrahedron j
                #        count += 1
                #        shared_nodes.append(node) #append the node that it is in common
                if len(shared_nodes) == 3: #after have done this for all the nodes, if the cound is 3 then calculate the shared area between the tetrahedrons
                    sc0 = nodecoords[node_indices[shared_nodes[0]],:]
                    #sc0 = gmsh.model.mesh.getNode(shared_nodes[0])[0] #coordinates of node 0
                    sc1 = nodecoords[node_indices[shared_nodes[1]],:]
                    #sc1 = gmsh.model.mesh.getNode(shared_nodes[1])[0] #coordinates of node 1
                    sc2 = nodecoords[node_indices[shared_nodes[2]],:]
                    #sc2 = gmsh.model.mesh.getNode(shared_nodes[2])[0] #coordinates of node 2
                    shared_area = np.linalg.norm(np.cross(sc2-sc0,sc1-sc0))/2 #compute shared area
                    shared_distance = np.linalg.norm(cell_center_i - cell_center[j])
                        #sqrt((abs(cell_center_i[0] - cell_center_j[0]))**2 + (abs(cell_center_i[1] - cell_center_j[1]))**2 + (abs(cell_center_i[2] - cell_center_j[2]))**2) #distance between volume elements
                    interior_tet[i, j] = shared_area/shared_distance #division between shared area and shared distance
                else:
                    shared_area = 0
                    interior_tet[i, j] = shared_area
    
    interior_tet_sum = np.sum(interior_tet, axis=1) #sum of interior_tet per columns (so per i element)
    
    return interior_tet, interior_tet_sum

#%%
###############################################################################
#CALCULATION OF 
###############################################################################
# Calculation surface areas
def surface_absorption_fun(vGroupsNames,df_abs,nBands,th,c0):
    """
    Assign the absoprtion term to each triangle boundary face depending on the surface material and the frequency

    Parameters
    ----------
        vGroupsNames : list 
            List of the names of the materials in the msh file (the material name are the same as the one assigned in the SketchUp file)
        df_abs : Dataframe 
            Absoprtion coefficient for each material assigned
        nBands : int 
            Number of frequency bands
        th : int
            Option for the absorption term; Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
        c0 : int 
            Speed of sound 

    Raises
    ------
        ValueError: If there are not values in the csv file, it will raise it as an error.

    Returns
    -------
        surface_absorption : list of tuples
            List of absorption term for each surface and for each frequency
        triangle_face_absorption : list of arrays 
            Absorption term value for each triangle face at the boundary
        absorption_coefficient_dict : dict
            Dictionary initialization of absorption coefficients per each surface and per each frequency 
    """
    #Initialize a list to store surface tags and their absorption coefficients
    surface_absorption = [] #initialization absorption term (alpha*surfaceofwall) for each wall of the room
    triangle_face_absorption = [] #initialization absorption term for each triangle face at the boundary and per each wall
    absorption_coefficient_dict = {}

    for group in vGroupsNames:
        if group[0] != 2:
            continue

        name_group = group[2]
        row = df_abs[df_abs["Material"] == name_group]
        if row.empty:
            raise ValueError(f"No absorption data found for surface: {name_group}")

        abscoeff = row.iloc[0, 1:].values.astype(float).tolist() #input(f"Enter absorption coefficient for frequency {fc_low} to {fc_high} for {name_abs_coeff}:") 
        absorption_coefficient_dict, surface_absorption = surface_materials (group, abscoeff, surface_absorption, absorption_coefficient_dict,nBands,th,c0)

    for entity, Abs_term in surface_absorption:
        triangle_faces, _ = gmsh.model.mesh.getElementsByType(2, entity) #Get all the triangle faces for the current surface
        triangle_face_absorption.extend([Abs_term] * len(triangle_faces)) #Append the Abs_term value for each triangle face

    return surface_absorption, triangle_face_absorption, absorption_coefficient_dict




#%%
###############################################################################
#CALCULATION OF ENTIRE SURFACE AREA PER MATERIAL
###############################################################################
# Calculation surface areas
def surface_area(surface_absorption, nodecoords, node_indices):
    """
    Calculation of total surface area for each model surface

    Parameters
    ----------
        surface_absorption : list of tuples 
            List of absorption term for each surface and for each frequency
        nodecoords : array of floats 
            The coordinates of each node in the mesh
        node_indices :dict 
            Indices of all the nodes in the mesh

    Returns
    -------
        surface_areas : dict
            Surface are for each surface of the model
    """
    surface_areas = {}   
    for entity, Abs_term in surface_absorption:   
        face_nodes_per_entity= gmsh.model.mesh.getElementFaceNodes(2, 3, tag=entity)
        surf_area_tot = 0
        for i in range(0, len(face_nodes_per_entity), 3): # per each element basically, goes trhough the nodes of each face 3by3
            #print(i)
            f = tuple(sorted(face_nodes_per_entity[i:i + 3])) 
            fc0 = nodecoords[node_indices[f[0]],:]
            #fc0 = gmsh.model.mesh.getNode(f[0])[0] #coordinates of vertix 0
            fc1 = nodecoords[node_indices[f[1]],:]
            fc2 = nodecoords[node_indices[f[2]],:]
            #fc1 = gmsh.model.mesh.getNode(f[1])[0] #coordinates of vertix 1
            #fc2 = gmsh.model.mesh.getNode(f[2])[0] #coordinates of vertix 2
            face_area = 0.5 * np.linalg.norm(np.cross(fc1 - fc0, fc2 - fc0)) #Compute the area using half of the cross product's magnitude
            surf_area_tot += face_area
            surface_areas[entity] = surf_area_tot
    return surface_areas

#%%
###############################################################################
#CALCULATION OF BOUNDARY ELEMENTS
###############################################################################
#FACE AREA & boundary_areas
def boundary_triang(velemNodes, nBands, bounNode, nodecoords, node_indices, triangle_face_absorption):
    """
    Assign the correct equaivalent area to each surface to each tetrahedra

    Parameters
    ----------
        velemNodes : array of int
            Indices of all the volumetric nodes per each colume element in the mesh
        nBands : int 
            Number of frequency bands
        bounNode : array of int 
            Indices of all the boundary nodes per each boundary surface in the mesh
        nodecoords : array of floats
            The coordinates of each node in the mesh
        node_indices : dict 
            Indices of all the nodes in the mesh
        triangle_face_absorption : list of arrays
            Absorption term value for each triangle face at the boundary

    Returns
    -------
        boundary_areas : array of floats 
            Product between the area and the correspondent absorption term for each surface for each tetrahedron
        total_boundArea : float
            Total surface area of the room
    """
    total_boundArea = 0  # initialization of total surface area of the room
    boundary_areas = []  # Initialize a list to store boundary_areas values for each tetrahedron
    import itertools
    face_areas = np.zeros(
        len(velemNodes))  # Per each tetrahedron, if there is a face that is on the boundary, include the area, otehrwise zero
    for idx, element in enumerate(velemNodes):  # for index and element in the number of tetrahedrons
        # if idx == 491:
        tetrahedron_boundary_areas = 0  # initialization tetrahedron face on boundary*its absorption term
        total_tetrahedron_boundary_areas = np.zeros(nBands)  # initialization total tetrahedron face on boundary*its absorption term if there are more than one face in the tetrahedron that is on the boundary
        # print(idx)
        node_combinations = [list(nodes) for nodes in itertools.combinations(element,3)]  # all possible combinations of the nodes of the tetrahedrons (it checks also for the order of the nodes in the same combination)
        # Check if the nodes are in any order in bounNode
        is_boundary = False  # variable to say that at the beginning the face in not on a boundary
        for nodes in node_combinations:  # for each node in each combination
            for surface_idx, surface in enumerate(bounNode):  # for index and surface in the number of nodes
                surface_set = sorted(set(surface))  # creates a set of the surface nodes
                surface_set_idx = surface_idx
                nodes_set = sorted(set(nodes))  # create a set of the node combination of the tetrahedron into consideration
                # surface_list = list(surface)
                if nodes_set == surface_set:  # if these are equal, it means that the tetrahedron into consideration has a surface in the boundary and therefore is_boundary gets the value of True.
                    # print(surface_set)
                    # print(surface_list)
                    is_boundary = True
                    if is_boundary:  # if the surface is at the boundary, then take the coordinates of each vertix
                        # Convert the vertices to NumPy arrays for vector operations
                        bc0 = nodecoords[node_indices[nodes[0]], :]
                        bc1 = nodecoords[node_indices[nodes[1]], :]
                        bc2 = nodecoords[node_indices[nodes[2]], :]

                        # bc0 = gmsh.model.mesh.getNode(nodes[0])[0] #coordinates of vertix 0
                        # bc1 = gmsh.model.mesh.getNode(nodes[1])[0] #coordinates of vertix 1
                        # bc2 = gmsh.model.mesh.getNode(nodes[2])[0] #coordinates of vertix 2

                        face_area = 0.5 * np.linalg.norm(np.cross(bc1 - bc0, bc2 - bc0))  # Compute the area using half of the cross product's magnitude
                        # print(face_area)

                        face_areas[idx] = face_area  # area of the surface that is on boundary per each tetrahedron
                        total_boundArea += face_area  # add to the total boundary area

                        if face_area > 0:
                            # Use the index to access the corresponding absorption area
                            face_absorption_product = face_area * triangle_face_absorption[surface_set_idx]  # calculate the product between the area*the correspondent absorption term
                            # print(face_absorption_product)

                            tetrahedron_boundary_areas += face_absorption_product  # add the calculation to the tetrahedron correspondent

                            total_tetrahedron_boundary_areas[:] = tetrahedron_boundary_areas  # if there are multiple surfaces on the boundary per each tetrahedron, then add also the second and the third one

        boundary_areas.append(np.array(
            total_tetrahedron_boundary_areas))  # Append the total boundary_areas for the tetrahedron to the list
    # print(total_tetrahedron_boundary_areas)
    boundary_areas = np.array(boundary_areas)
    boundary_areas = boundary_areas.T
    return boundary_areas, total_boundArea


#%%
###############################################################################
#CALCULATION OF EQUIVALENT ABSORPTION AREA
###############################################################################
def equiv_absorp_area(cell_volume, total_boundArea, surface_areas_in, absorption_coefficient_dict):
    """
    Calculation of the eqauivalent absorption area

    Parameters
    ----------
        cell_volume : array of floats
            The volume of each element cell
        total_boundArea : float 
            Total surface area of the room
        surface_areas_in : dict
            Surface are for each surface of the model
        absorption_coefficient_dict : dict 
            Dictionary of absorption coefficients per each surface and per each frequency

    Returns
    -------
        V : float
            Volume of the room
        S : float
            Total surface are of the room
        Eq_A : array of floats
            Equivalent absorption area
    """
    V = sum(cell_volume) #volume of the room
    S = total_boundArea #total surface area of the room
    
    sum_alpha_average = 0
    Eq_A = 0
    #Absorption parameters for room
    for entity in surface_areas_in:
        #print(entity)
        sum_alpha_average += np.multiply(absorption_coefficient_dict[entity],surface_areas_in[entity])
        Eq_A += np.multiply(absorption_coefficient_dict[entity],surface_areas_in[entity])
    #alpha_average = sum_alpha_average/S #average absorption
    
    return V,S,Eq_A


def calculation_sourceon_time(nBands, V, Eq_A):
    """
    Calculation of the time the sources stays on

    Parameters
    ----------
        nBands : int
            Number of frequency bands
        V : float
            Volume of the room
        Eq_A : array of floats
            Equivalent absorption

    Returns
    -------
        sourceon_time : float
            Time that the source stays on
    """
    RT_Sabine_band = []
    for iBand in range(nBands):
        #freq = center_freq[iBand]
        RT_Sabine = 0.16*V/Eq_A[iBand]
        RT_Sabine_band.append(RT_Sabine)
    
    sourceon_time = round(max(RT_Sabine_band),1)#time that the source is ON before interrupting [s]

    return sourceon_time


#%%
###############################################################################
#CALCULATION OF SABINE RT, TOTAL RECORDING TIME AND SOURCE ON TIME
###############################################################################

def calculation_rec_time(sourceon_time, dt, edt=-1, ir_length=-1.0):
    """
    Calculation of the simulation time run

    Parameters
    ----------
        sourceon_time : float
            Time that the source stays on
        dt : float
            Time step
        edt : int, optional
            Early decay time. Defaults to -1.
        ir_length : float, optional
            Impulse response length. Defaults to -1.

    Returns
    -------
        recording_time : float
            Length of the simulaton run
        t : array of floats
            Time array of time steps
        recording_steps : int
            Number of time steps
    """
    if (edt == -1) and (ir_length == -1.0):
        recording_time = sourceon_time * 2
    elif edt != -1:
        recording_time = sourceon_time + (sourceon_time / 60 * edt)
    elif ir_length != -1.0:
        recording_time = sourceon_time + ir_length

    # Time resolution
    t = np.arange(0, recording_time, dt)  # mesh point in time
    recording_steps = ceil(recording_time / dt)  # number of time steps to consider in the calculation
        
    return recording_time, t, recording_steps

#%%
###############################################################################
#CALCULATION OF DIFFUSION PARAMETERS
###############################################################################
#Diffusion parameters
def diffusion_coeff(V, S, c0):
    """
    Calculation of the theoretical diffusion coefficient

    Parameters
    ----------
        V : float
            Volume of the room
        S : float
            Total surface are of the room
        c0 : int
            Speed of sound 

    Returns
    -------
        Dx : float 
            Diffusion coefficient in the x direction (equal to the theoretical diffusion coefficient)
        Dy : float 
            Diffusion coefficient in the y direction (equal to the theoretical diffusion coefficient)
        Dz : float
            Diffusion coefficient in the z direction (equal to the theoretical diffusion coefficient)
    """
    mean_free_path = (4*V)/S #mean free path for 3D
    #mean_free_time = mean_free_path/c0 #mean free time for 3D
    #mean_free_time_step = int(mean_free_time/dt)
    Dx = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms x direction
    Dy = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms y direction
    Dz = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms z direction
    return Dx, Dy, Dz

#%%
###############################################################################
#CALCULATION OF SOURCE & RECEIVER DISTANCE
###############################################################################

def distance_source_receiver(coord_rec, coord_source):
    """
    Calculation of the distance between source and receiver positions

    Parameters
    ----------
        coord_rec : list 
            Coordinates of the receiver position
        coord_source : list
            Coordinates of the source position

    Returns
    -------
        dist_sr : float
            Distance between source and receiver position
    """
    #distance between source and receiver
    dist_sr = math.sqrt((abs(coord_rec[0] - coord_source[0]))**2 + (abs(coord_rec[1] - coord_source[1]))**2 + (abs(coord_rec[2] - coord_source[2]))**2) #distance between source and receiver
    return dist_sr

#%%
###############################################################################
#SOURCE INTERPOLATION
###############################################################################
def source_interp(cell_center, coord_source):
    """
    Interpolation of the source position

    Parameters
    ----------
        cell_center : array of floats
            The coordinates of the center of the cell element per each element
        coord_source : list
            Coordinates of the source position

    Returns
    -------
        cl_tet_s_keys : Dict keys
            Closest tetrahedrons Indices to the source
        total_weights_s : dict
            Weights for each 4 closest points to the source position for interpolation
    """
    
    #SOURCE INTERPOLATION CALCULATED WITHIN 4 CENTRE CELL SELECTED (TETRAHEDRON)
    #Position of source is the centre of a cell so the minimum distance with the centre of a cell has been calculated to understand which cell is the closest
    dist_source_cc_list = [] #initialise the list for all the distances between each cell centre and the source
    for i in range(len(cell_center)): #for each tetra
        dist_source_cc = math.sqrt(np.sum((cell_center[i] - coord_source)**2)) #calculate the distance between its centre cell and the source coordinate
        dist_source_cc_list.append(dist_source_cc) #append the distance in a list
    #source_idx = np.argmin(dist_source_cc_list) #take the minimum distance index; this is where the source will be positioned
    
    dist_source_cc_list_sorted = sorted(dist_source_cc_list) #sorted from the minimum to the maximum distance
    selected_source_cc_list = dist_source_cc_list_sorted[:4] #take only the first four element of the sorted list (take the first 4 cell centres closest to the source)
    
    dist_source_cc_list_sorted_indices = np.argsort(dist_source_cc_list)[:4] #takes the Indices of the minimum distances
    #selected_source_cc_list_indices = dist_source_cc_list_sorted_indices[:4] #does exactly the same as the previous line
    
    #cl_tet_s stands for cl=closest, tet=tetrahedron, s=to the source
    cl_tet_s = {} #initialise dictionary for closest tetrahedrons to the source
    for i in range(len(dist_source_cc_list_sorted_indices)): #for each of the 4 tetra closest to the source
        idx = selected_source_cc_list[i] #take its index
        cl_tet_s[dist_source_cc_list_sorted_indices[i]] = idx #put it in a dictionary
    
    total_weights_s = {} #initialise weights for each tetrahedron around the actual source position
    sum_weights_s = 0
    #Vs = 0
    for i, dist in cl_tet_s.items(): #for each key and value in the dictionary (so for each closest tetrahedron to the source)
        weights = np.divide(1.0 , dist)  #calculate the inverse distance weights, so closer to the point means higher weight
        #print(weights)
        sum_weights_s += weights
        #weights /= np.sum(weights)  # Normalize weights to sum to 1
        total_weights_s[i] = weights #put the wweigths (values) to the correspondent closest tetrahedron (keys)
        #Vs += cell_volume[i] #volume of the source calculated summing the volumes of all the tetrahedrons involved
    
    #total_weights_s_values = total_weights_s.values()
    for i,weight in total_weights_s.items():
        total_weights_s[i] = weight/sum_weights_s if sum_weights_s != 0 else 0
    
    cl_tet_s_keys = cl_tet_s.keys() #take only the keys of the cl_tet_s dictionary (so basically the Indices of the tetrahedrons)

    return cl_tet_s_keys, total_weights_s

#%%
###############################################################################
#CALCULATION OF SOURCE VOLUME
###############################################################################

def source_volume(velemNodes, nodecoords, coord_source, cell_volume):
    """
    Calculation of the volume of the source

    Parameters
    ----------
        velemNodes : array of int
            Indices of all the volumetric nodes per each colume element in the mesh
        nodecoords : array of floats
            The coordinates of each node in the mesh
        coord_source : list
            Coordinates of the source position
        cell_volume : array of floats
            The volume of each element cell

    Returns
    -------
        Vs : float
            Volume of the source
    """
    #To make sure that the source is in the correct tetrahedron position
    node_ids = velemNodes.T
    ori=nodecoords[node_ids[0,:]-1,:]
    #ori = nodecoords[node_indices[node_ids[:,0],:]]
    v_tet_s1=nodecoords[node_ids[1,:]-1,:]-ori
    v_tet_s2=nodecoords[node_ids[2,:]-1,:]-ori
    v_tet_s3=nodecoords[node_ids[3,:]-1,:]-ori
    n_tet=len(node_ids.T)
    v1s = v_tet_s1.T.reshape((3,1,n_tet))
    v2s = v_tet_s2.T.reshape((3,1,n_tet))
    v3s = v_tet_s3.T.reshape((3,1,n_tet))
    mat = np.concatenate((v1s,v2s,v3s), axis=1)
    inv_mat = np.linalg.inv(mat.T).T
    #if coord_source.size==3:  # to make rec has a dimension of (N_rec,3)
    #    rec=rec.reshape((1,3))
    coord_source_array = np.array(coord_source)
    if coord_source_array.size==3:  # to make rec has a dimension of (N_rec,3)
        coord_source_array=coord_source_array.reshape((1,3))
    N_sou=coord_source_array.shape[0]
    orir=np.repeat(ori[:,:,np.newaxis], N_sou, axis=2)
    newp=np.einsum('imk,kmj->kij',inv_mat,coord_source_array.T-orir)
    val=np.all(newp>=0, axis=1) & np.all(newp <=1, axis=1) & (np.sum(newp, axis=1)<=1)
    id_tet, id_p = np.nonzero(val)
    res = -np.ones(N_sou, dtype=id_tet.dtype) # Sentinel value
    res[id_p]=id_tet
        
    #VOLUME ORIGINAL
    Vs = cell_volume[res[0]] #volume of the source = to volume of cell where the source is 
    # Vs = 1
    
    return Vs

#%%
###############################################################################
#INITIAL CONDITIONS
###############################################################################
#Initial condition - Source Info (interrupted method)
def initial_cond(Ws, Vs, sourceon_time, dt, recording_steps):
    """
    Definition of initial conditions

    Parameters
    ----------
        Ws : float
            Power of the source
        Vs : float
            Volume of the source
        sourceon_time : float
            Time that the source stays on
        dt : float
            Time step
        recording_steps : int
            Number of time steps

    Returns
    -------
        source1 : array of floats
            Energy density of source number 1 at each time step position
        sourceon_steps : int
            Number of time steps while the source is on
    """
    w1=Ws/Vs #w1 = round(Ws/Vs,4) #power density of the source [Watts/(m^3))]
    sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
    s1 = np.multiply(w1,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position
    source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps
    return source1, sourceon_steps

#%%
###############################################################################
#DEFINITION OF SOURCE MATRIX
###############################################################################
def source_matrix(voluEl,cl_tet_s_keys, source1, total_weights_s):
    """
    Calculation of the source matrix

    Parameters
    ----------
        voluEl : array of int
            Indices of all the volume elements (tetrahedra) in the mesh
        cl_tet_s_keys : Dict keys
            Closest tetrahedrons Indices to the source
        source1 : array of floats
            Energy density of source number 1 at each time step position
        total_weights_s : dict
            Weights for each 4 closest points to the source position for interpolation

    Returns
    -------
        s : array of floats
            Matrix of tretrahedron central points inserting source energy
    """

    #INTERPOLATION WITH CELL CENTRES - SOURCE MATRIX
    s = np.zeros((len(voluEl))) #matrix of zeros for source
    for tet_s in cl_tet_s_keys:
        s[tet_s] = source1[0] *total_weights_s[tet_s]
    
    return s

#%%
###############################################################################
#RECEIVER INTERPOLATION
###############################################################################
def receiver_interp(cell_center, coord_rec):
    """
    Interpolation of the receiver position

    Parameters
    ----------
        cell_center : array of floats
            The coordinates of the center of the cell element per each element
        coord_rec : list
            Coordinates of the receiver position

    Returns
    -------
        cl_tet_r_keys : dict keys
            Closest tetrahedrons Indices to the receiver
        total_weights_r : dict
            Weights for each 4 closest points to the receiver position for interpolation
    """
    
    #RECEIVER INTERPOLATION CALCULATED WITHIN 4 CENTRE CELL SELECTED (TETRAHEDRON)
    #Position of receiver is the centre of a cell so the minimum distance with the centre of a cell has been calculated to understand which cell is the closest
    dist_rec_cc_list = [] #initialise the list for all the distances between each cell centre and the source
    for i in range(len(cell_center)): #for each tetra
        dist_rec_cc = math.sqrt(np.sum((cell_center[i] - coord_rec)**2)) #calculate the distance between its centre cell and the source coordinate
        dist_rec_cc_list.append(dist_rec_cc) #append the distance in a list
    #rec_idx = np.argmin(dist_rec_cc_list) #take the minimum distance index; this is where the source will be positioned
    
    dist_rec_cc_list_sorted = sorted(dist_rec_cc_list) #sorted from the minimum to the maximum distance
    selected_rec_cc_list = dist_rec_cc_list_sorted[:4] #take only the first four element of the sorted list (take the first 4 cell centres closest to the source)
    
    dist_rec_cc_list_sorted_indices = np.argsort(dist_rec_cc_list)[:4] #takes the Indices of the minimum distances
    #selected_rec_cc_list_indices = dist_rec_cc_list_sorted_indices[:4] #does exactly the same as the previous line
    
    #cl_tet_s stands for cl=closest, tet=tetrahedron, s=to the source
    cl_tet_r = {} #initialise dictionary for closest tetrahedrons to the source
    for i in range(len(dist_rec_cc_list_sorted_indices)): #for each of the 4 tetra closest to the source
        idx = selected_rec_cc_list[i] #take its index
        cl_tet_r[dist_rec_cc_list_sorted_indices[i]] = idx #put it in a dictionary
    
    total_weights_r = {} #initialise weights for each tetrahedron around the actual source position
    sum_weights_r = 0
    #Vs = 0
    for i, dist in cl_tet_r.items(): #for each key and value in the dictionary (so for each closest tetrahedron to the source)
        weights = np.divide(1.0 , dist)  #calculate the inverse distance weights, so closer to the point means higher weight
        #print(weights)
        sum_weights_r += weights
        #weights /= np.sum(weights)  # Normalize weights to sum to 1
        total_weights_r[i] = weights #put the wweigths (values) to the correspondent closest tetrahedron (keys)
        #Vs += cell_volume[i] #volume of the source calculated summing the volumes of all the tetrahedrons involved
    
    #total_weights_s_values = total_weights_s.values()
    for i,weight in total_weights_r.items():
        total_weights_r[i] = weight/sum_weights_r if sum_weights_r != 0 else 0
    
    cl_tet_r_keys = cl_tet_r.keys() #take only the keys of the cl_tet_s dictionary (so basically the Indices of the tetrahedrons)
    
    return cl_tet_r_keys, total_weights_r

#%%
###############################################################################
#CALCULATION OF LENGTH, WIDTH AND HEIGHT OF ROOM
###############################################################################

def room_dimensions(nodecoords):
    """
    Calculation of the room dimensions

    Parameters
    ----------
        nodecoords : array of floats
            The coordinates of each node in the mesh

    Returns
    -------
        room_length : float
            Length of the room
        room_width : float
            Width of the room
        room_height : float
            Height of the room
    """
    #Extract x-coordinates of all nodes LENGTH
    x_coordinates = nodecoords[:, 0]
    #Find the minimum and maximum x-coordinates to determine the length of the room
    min_x = np.min(x_coordinates)
    max_x = np.max(x_coordinates)
    #Calculate the length of the room
    room_length = max_x - min_x

    #Extract y-coordinates of all nodes WIDTH
    y_coordinates = nodecoords[:, 1]
    #Find the minimum and maximum x-coordinates to determine the width of the room
    min_y = np.min(y_coordinates)
    max_y = np.max(y_coordinates)
    #Calculate the width of the room
    room_width = max_y - min_y
    
    #Extract y-coordinates of all nodes HEIGHT
    z_coordinates = nodecoords[:, 2]
    #Find the minimum and maximum x-coordinates to determine the width of the room
    min_z = np.min(z_coordinates)
    max_z = np.max(z_coordinates)
    #Calculate the height of the room
    room_height = max_z - min_z
    
    return room_length, room_width, room_height
    
#%%
###############################################################################
#CALCULATION OF LINE RECEIVERS
###############################################################################
def line_receivers(room_length, room_width, coord_rec, coord_source, cell_center):
    """
    Definition of a line point receiver

    Parameters
    ----------
        room_length : float
            Length of the room
        room_width : float
            Width of the room
        coord_rec : list
            Coordinates of the receiver position
        coord_source : list
            Coordinates of the source position
        cell_center : array of floats
            The coordinates of the center of the cell element per each element

    Returns
    -------
        x_axis : array of floats
            Linspace on x_axis with distance dx
        y_axis : array of floats
            Linspace on y_axis with distance dy
        line_rec_x_idx_list : list
            Indices of center cells close to the line receiver position in the x axis
        dist_x : array of floats 
            Distance between each line receiver point cell and the source in the x axis
        line_rec_y_idx_list : list
            Indices of center cells close to the line receiver position in the y axis
        dist_y : array of floats 
            Distance between each line receiver point cell and the source in the y axis
    """
    #Arange linespace lines
    dx = 0.5
    x_axis = np.arange(0,room_length+dx,dx) #lispace on x_axis with distance dx
    y_axis = np.arange(0,room_width+dx,dx)

    #RECEIVERS IN A X LINE
    line_rec_x_idx_list = []
    dist_x = np.array([])
    for x_chang in x_axis:
        line_rec = [x_chang, coord_rec[1], coord_rec[2]]
        #Position of line_receiver is the centre of a cell
        dist_line_rec_x =  math.sqrt((abs(line_rec[0] - coord_source[0]))**2 + (abs(line_rec[1] - coord_source[1]))**2 + (abs(line_rec[2] - coord_source[2]))**2) #distance between source and line_receiver
        dist_x = np.append(dist_x, dist_line_rec_x)  # Append to the NumPy array
        dist_line_rec_x_cc_list = []
        for i in range(len(cell_center)):
            dist_line_rec_x_cc = math.sqrt(np.sum((cell_center[i] - line_rec)**2))
            dist_line_rec_x_cc_list.append(dist_line_rec_x_cc)
        line_rec_x_idx = np.argmin(dist_line_rec_x_cc_list)
        line_rec_x_idx_list.append(line_rec_x_idx)   
    
    #RECEIVERS IN A Y LINE
    line_rec_y_idx_list = []
    dist_y = np.array([])
    for y_chang in y_axis:
        line_rec = [coord_rec[0], y_chang, coord_rec[2]]
        #Position of line_receiver is the centre of a cell
        dist_line_rec_y =  math.sqrt((abs(line_rec[0] - coord_source[0]))**2 + (abs(line_rec[1] - coord_source[1]))**2 + (abs(line_rec[2] - coord_source[2]))**2) #distance between source and line_receiver
        dist_y = np.append(dist_y, dist_line_rec_y)  # Append to the NumPy array
        dist_line_rec_y_cc_list = []
        for i in range(len(cell_center)):
            dist_line_rec_y_cc = math.sqrt(np.sum((cell_center[i] - line_rec)**2))
            dist_line_rec_y_cc_list.append(dist_line_rec_y_cc)
        line_rec_y_idx = np.argmin(dist_line_rec_y_cc_list)
        line_rec_y_idx_list.append(line_rec_y_idx)  
    
    return x_axis, y_axis, line_rec_x_idx_list, dist_x, line_rec_y_idx_list, dist_y

#%%
###############################################################################
#CALCULATION OF BETA_ZERO
###############################################################################

def beta_zero_freq_fun(boundary_areas, dt, Dx, interior_tet_sum, cell_volume):
    """
    Calculation of factor beta_zero

    Parameters
    ----------
        boundary_areas : array of floats
            Product between the area and the correspondent absorption term for each surface for each tetrahedron
        dt : float
            Time step
        Dx : float
            Diffusion coefficient in the x direction (equal to the theoretical diffusion coefficient)
        interior_tet_sum : array of floats
            Sum of interior_tet per columns
        cell_volume : array of floats
            The volume of each element cell

    Returns
    -------
        beta_zero_freq : list
            Coefficient beta zero used in the calculation of the energy density per each frequency
    """
    beta_zero_freq = []
    for iBand in range(len(boundary_areas)):
        #print(iBand)
        #freq = center_freq[iBand]
        #print(boundary_areas[iBand])
        beta_zero_element = np.divide(dt*((Dx *interior_tet_sum) + boundary_areas[iBand]),cell_volume) #my interpretation of the beta_zero
        beta_zero_freq.append(beta_zero_element)
        
    return beta_zero_freq

#%%
###############################################################################
#MAIN CALCULATION - COMPUTING ENERGY DENSITY
############################################################################### 

def computing_energy_density(nBands, voluEl, recording_steps, beta_zero_freq, dt, c0, m_atm, Dx, interior_tet, cell_volume, s, cl_tet_r_keys, total_weights_r, tcalc, cl_tet_s_keys, source1, total_weights_s, t, sourceon_time, rho, json_file=None, use_gpu=False):
    """
    Computation of energy density

    Parameters
    ----------
        nBands : int
            Number of frequency bands
        voluEl : array of int
            Indices of all the volume elements (tetrahedra) in the mesh
        recording_steps : int 
            Number of time steps
        beta_zero_freq : list
            Coefficient beta zero used in the calculation of the energy density per each frequency
        dt : float
            Time step
        c0 : int 
            Speed of sound
        m_atm : float
            Air absorption coefficient
        Dx : float
            Diffusion coefficient in the x direction (equal to the theoretical diffusion coefficient)
        interior_tet : array of floats
            Matrix of tetrahedron per tetrahedron of the division between shared area and shared distance
        cell_volume : array of floats 
            The volume of each element cell
        s : array of floats 
            Matrix of tretrahedron central points inserting source energy
        cl_tet_r_keys : dict keys
            Closest tetrahedrons Indices to the receiver
        total_weights_r : dict
            Weights for each 4 closest points to the receiver position for interpolation
        tcalc : str
            Type of calculation; "decay" if the source switches off and "stationarysource" if the source is stationary
        cl_tet_s_keys : dict keys
            Closest tetrahedrons Indices to the source
        total_weights_s : dict
            Weights for each 4 closest points to the source position for interpolation
        source1 : array of floats
            Energy density of source number 1 at each time step position
        total_weights_s : dict
            Weights for each 4 closest points to the source position for interpolation
        t : array of floats
            Time array of time steps
        sourceon_time : float
            Time that the source stays on
        rho : float
            Density of air

    Returns
    -------
        w_new_band : list of arrays
            Energy density at the time step n+1 at each centre cell per each frequency band
        w_rec_band : list of arrays
            Energy density over time at the receiver position per each frequency band
        w_rec_off_band : list of arrays
            Energy density over time after the source is switched off at the receiver position per each frequency band
        w_rec_off_deriv_band : list of arrays
            Derivative of the energy density over time after the source is switched off at the receiver position per each frequency band
        p_rec_off_deriv_band : list of arrays
            Derivative of the pressure over time after the source is switched off at the receiver position per each frequency band
        idx_w_rec : int
            Time index at which the source is switched off
        t_off : array of floats
            Time array after the source is switched off
    """

    result_container = {}
    if json_file is not None:
        with open(json_file, "r") as f:
            result_container = json.load(f)

    called_from_choras = False
    if "results" in result_container:
        called_from_choras = True

    # ── CPU optimisation: sparse matvec ──
    # interior_tet is built as a dense (N, N) ndarray, but only ~5 entries per
    # row are nonzero (one per face-neighbour tet). Converting once to scipy
    # CSR turns the O(N**2) dense matvec into O(nnz). On meshes with 1000+
    # tets this is 50-200x faster per step and dominates the win on CPU.
    # Skip the conversion for tiny meshes (<= 64 tets) where dense BLAS wins.
    from scipy.sparse import issparse as _issparse, csr_matrix as _csr_matrix
    if not _issparse(interior_tet) and len(voluEl) > 64:
        nnz = int(np.count_nonzero(interior_tet))
        density = nnz / (len(voluEl) ** 2 + 1)
        if density < 0.10:
            interior_tet = _csr_matrix(interior_tet)
            print(f"acousticDE: interior_tet converted to CSR sparse "
                  f"(nnz={nnz}, density={density:.2%})")

    # CPU optimisation: pre-compute receiver / source index + weight arrays
    # once, replace per-step Python loops with one numpy gather + reduction.
    _r_idx_np = np.asarray(list(cl_tet_r_keys), dtype=np.int64)
    _r_w_np   = np.asarray([total_weights_r[k] for k in cl_tet_r_keys])
    _s_idx_np = np.asarray(list(cl_tet_s_keys), dtype=np.int64)
    _s_w_np   = np.asarray([total_weights_s[k] for k in cl_tet_s_keys])

    # ── Optional GPU path ──
    # When use_gpu=True we run the inner per-time-step math on the GPU via
    # CuPy. The CPU branch below is preserved verbatim so existing users get
    # bit-identical output. GPU output is float-equivalent (within reduction-
    # order epsilon, ~1e-12 relative for float64). If CuPy import fails (no
    # CuPy installed, no CUDA available), we transparently fall back to CPU.
    _xp = None
    _cusp = None
    _fused_kernel = None
    if use_gpu:
        try:
            import cupy as _xp                         # noqa: F401
            import cupyx.scipy.sparse as _cusp         # noqa: F401
            # Single-launch fused CUDA kernel that runs *all* time steps of
            # one frequency band in one kernel call. Eliminates the per-step
            # Python+CuPy dispatch overhead (~1 ms per step) that dominates a
            # naive drop-in port.
            #
            # Cooperative-groups variant: grid-stride loop across tets +
            # grid.sync() for the cross-block barrier on each step. Requires
            # GPU compute capability >= 6.0 (sm_60+); RTX 2060 is sm_75. The
            # kernel must be launched cooperatively (cuda.LaunchCooperative)
            # for grid.sync() to work. No 1024-tet cap.
            #
            # Memory: 4 small global arrays (w, w_old, s, w_rec) + the CSR
            # matrix. Receiver reduction uses an atomicAdd into w_rec[step]
            # instead of shared memory (because n_r is tiny, ~4, and we no
            # longer have a single-block shared-mem region across the grid).
            _fused_kernel = _xp.RawKernel(r"""
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

extern "C" __global__ void fvm_band_steps(
    double* __restrict__ w,
    double* __restrict__ w_old,
    double* __restrict__ s,
    double* __restrict__ w_rec,
    const double* __restrict__ M_data,
    const int*    __restrict__ M_indices,
    const int*    __restrict__ M_indptr,
    const double* __restrict__ cell_volume,
    double*       __restrict__ wTEMP_g,        // global scratch [N_tets]
    const double* __restrict__ source1,
    const int*    __restrict__ r_idx,
    const double* __restrict__ r_w,
    const int*    __restrict__ s_idx,
    const double* __restrict__ s_w,
    const double c1, const double c2, const double c3, const double c4,
    const int N_tets, const int N_steps,
    const int n_r, const int n_s,
    const int decay_mode
) {
    cg::grid_group grid = cg::this_grid();
    const int gid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x * blockDim.x;

    for (int step = 0; step < N_steps; ++step) {
        // 1) sparse matvec: grid-stride loop over rows
        for (int i = gid; i < N_tets; i += stride) {
            double acc = 0.0;
            const int rs = M_indptr[i];
            const int re = M_indptr[i + 1];
            for (int j = rs; j < re; ++j) {
                acc += M_data[j] * w[M_indices[j]];
            }
            wTEMP_g[i] = acc;
        }
        grid.sync();

        // 2) per-tet update: w_new = c1*w_old - c2*w + c3*M@w/V + c4*s
        // We compute new_val into a private register, sync, then write
        // both w_old <- w and w <- new_val in a second pass.
        // (We can't write w[i] = new_val directly because other threads
        // are still reading w[i] for receiver interpolation. So: 1st pass
        // computes; 2nd pass writes; the grid sync between separates them.)
        // We stash new_val in wTEMP_g[i] since wTEMP_g is no longer needed.
        for (int i = gid; i < N_tets; i += stride) {
            wTEMP_g[i] = c1 * w_old[i]
                       - c2 * w[i]
                       + c3 * (wTEMP_g[i] / cell_volume[i])
                       + c4 * s[i];
        }
        grid.sync();

        // 3) swap: w_old <- w; w <- new_val (stashed in wTEMP_g)
        for (int i = gid; i < N_tets; i += stride) {
            w_old[i] = w[i];
            w[i] = wTEMP_g[i];
        }
        grid.sync();

        // 4) receiver interpolation w_rec[step] = sum(w[r_idx[i]] * r_w[i])
        // n_r is small (typically 4). One atomicAdd per thread per step is
        // cheap. Zero the slot first via thread 0.
        if (gid == 0) {
            w_rec[step] = 0.0;
        }
        grid.sync();
        if (gid < n_r) {
            atomicAdd(&w_rec[step], w[r_idx[gid]] * r_w[gid]);
        }
        grid.sync();

        // 5) source update at source tets
        if (gid < n_s) {
            const double sv = decay_mode ? source1[step] : source1[0];
            s[s_idx[gid]] = sv * s_w[gid];
        }
        grid.sync();
    }
}
""", "fvm_band_steps", options=("--std=c++14",))
            from scipy.sparse import issparse as _issparse, csr_matrix as _csr_matrix
            print("acousticDE: computing_energy_density running on GPU (CuPy + fused kernel)")
        except Exception as _e:
            print(f"acousticDE: use_gpu=True but CuPy unavailable ({_e}); "
                  f"falling back to CPU")
            use_gpu = False
            _xp = None
            _cusp = None
            _fused_kernel = None

    w_new_band = []
    w_rec_band = []
    w_rec_off_band = []
    w_rec_off_deriv_band = []
    p_rec_off_deriv_band = []

    prev_percent_done = 0

    for iBand in range(nBands):

        w_new = np.zeros(len(voluEl)) #unknown w at new time level (n+1)
        #w_old = np.zeros(len(voluEl))
        w = w_new #w at n level
        w_old = w #w_old at n-1 level
        #w[source_idx] = w1 #w (m time step) at source position -> impulse source

        w_rec = np.zeros(recording_steps) #energy density at the receiver

        # Pre-compute the four time-invariant coefficients once per band (CPU
        # version computed these inside the per-step formula every step --
        # algebraically identical, micro-optimisation only).
        _beta = beta_zero_freq[iBand]
        _denom = 1.0 + _beta
        _c1 = (1.0 - _beta) / _denom
        _c2 = 2.0 * dt * c0 * m_atm / _denom
        _c3 = 2.0 * dt * Dx / _denom
        _c4 = 2.0 * dt / _denom

        fused_done = False
        if use_gpu:
            # One-time host->device transfers for this band. interior_tet may
            # be a dense ndarray or a scipy sparse matrix -- handle both.
            interior_tet_d = (_cusp.csr_matrix(interior_tet)
                              if _issparse(interior_tet)
                              else _xp.asarray(interior_tet))
            cell_volume_d = _xp.asarray(cell_volume)
            s_d           = _xp.asarray(s)
            # source1 is a 1-D time series of length recording_steps
            source1_d     = _xp.asarray(np.asarray(source1))
            # Receiver / source index lists + per-tet weights, on device once
            _r_keys = list(cl_tet_r_keys)
            _s_keys = list(cl_tet_s_keys)
            r_idx_d = _xp.asarray(_r_keys, dtype=_xp.int64)
            s_idx_d = _xp.asarray(_s_keys, dtype=_xp.int64)
            r_w_d   = _xp.asarray([total_weights_r[k] for k in _r_keys])
            s_w_d   = _xp.asarray([total_weights_s[k] for k in _s_keys])
            w_d     = _xp.asarray(w)
            w_old_d = _xp.asarray(w_old)

            # Fused-kernel fast path: one CUDA launch runs all recording_steps
            # of this band, eliminating per-step Python+CuPy dispatch overhead.
            # Requirements: CSR interior_tet (sparse) and N_tets <= 1024 (one
            # CUDA block). Otherwise fall back to the dispatch-bound per-step
            # CuPy path below, which is still correct.
            N_tets = len(voluEl)
            use_fused = (
                _fused_kernel is not None
                and _issparse(interior_tet)
            )
            if iBand == 0:
                print(f"acousticDE: GPU path = "
                      f"{'fused cooperative-launch kernel' if use_fused else 'per-step CuPy (dispatch-bound)'} "
                      f"(N_tets={N_tets}, sparse={_issparse(interior_tet)})")
            if use_fused:
                M_data_d    = interior_tet_d.data.astype(_xp.float64, copy=False)
                M_indices_d = interior_tet_d.indices.astype(_xp.int32, copy=False)
                M_indptr_d  = interior_tet_d.indptr.astype(_xp.int32, copy=False)
                r_idx_i32_d = r_idx_d.astype(_xp.int32, copy=False)
                s_idx_i32_d = s_idx_d.astype(_xp.int32, copy=False)
                w_rec_d     = _xp.zeros(recording_steps, dtype=_xp.float64)
                wTEMP_d     = _xp.empty(N_tets, dtype=_xp.float64)  # global scratch
                decay_mode  = 1 if tcalc == "decay" else 0

                # Cooperative launch: pick a (block_dim, grid_dim) sized to
                # cover N_tets with a grid-stride loop. The CUDA cooperative-
                # launch contract requires grid_dim * block_dim threads to be
                # co-resident; CuPy uses cuLaunchCooperativeKernel under the
                # hood when enable_cooperative_groups=True.
                block_dim = 256
                grid_dim  = (N_tets + block_dim - 1) // block_dim
                grid_dim  = max(grid_dim, 1)

                _fused_kernel(
                    (grid_dim,), (block_dim,),
                    (w_d, w_old_d, s_d, w_rec_d,
                     M_data_d, M_indices_d, M_indptr_d,
                     cell_volume_d, wTEMP_d, source1_d,
                     r_idx_i32_d, r_w_d, s_idx_i32_d, s_w_d,
                     np.float64(_c1), np.float64(_c2), np.float64(_c3), np.float64(_c4),
                     np.int32(N_tets), np.int32(recording_steps),
                     np.int32(len(_r_keys)), np.int32(len(_s_keys)),
                     np.int32(decay_mode)),
                    enable_cooperative_groups=True,
                )
                _xp.cuda.Stream.null.synchronize()
                w_rec[:] = w_rec_d.get()
                s[:]     = s_d.get()
                w_new    = w_d.get()
                fused_done = True

        # Iteration range: if the fused kernel already filled w_rec, just run
        # one "iteration" at steps = recording_steps - 1 to trigger the
        # post-step receiver-off / derivative computation (which only ever
        # used the FINAL iteration's values anyway).
        _step_iter = ([recording_steps - 1]
                      if (fused_done and recording_steps > 0)
                      else range(0, recording_steps))

        #Computing w;
        for steps in _step_iter:
            #Compute w at inner mesh points
            #time_steps = steps*dt #total time for the calculation

            #Computing w_new (w at n+1 time step)

            if fused_done:
                # Already populated by the single fused-kernel launch above.
                # Skip the inner math; the post-step receiver-off block runs
                # on the already-filled w_rec / w_new state.
                pass
            elif use_gpu:
                # Same mathematical update as the CPU branch, on device:
                #   w_new = c1*w_old - c2*w + c3*(M @ w)/V + c4*s
                wTEMP_d = interior_tet_d @ w_d
                w_new_d = _c1 * w_old_d - _c2 * w_d \
                          + _c3 * (wTEMP_d / cell_volume_d) + _c4 * s_d
                w_old_d = w_d
                w_d = w_new_d
                # Receiver interpolation: gather + dot, one device->host copy
                w_rec[steps] = float((w_new_d[r_idx_d] * r_w_d).sum().get())
                # Source update at the source tets
                if tcalc == "decay":
                    s_d[s_idx_d] = source1_d[steps] * s_w_d
                elif tcalc == "stationarysource":
                    s_d[s_idx_d] = source1_d[0] * s_w_d
                # Keep CPU-side w_new available for percentage / end-of-band
                # paths (only fetched once per band)
                w_new = None  # not used per-step on GPU path
            else:
                # CPU step. Same math as original, but:
                #   - per-band coefficients (_c1.._c4) precomputed outside loop
                #   - interior_tet is CSR sparse on non-tiny meshes (fast matvec)
                #   - receiver / source interpolation use numpy gather instead
                #     of Python for-loops over a few-element list
                wTEMP = interior_tet @ w
                w_new = _c1 * w_old - _c2 * w + _c3 * (wTEMP / cell_volume) + _c4 * s

                #Update w before next step
                w_old = w #The w at n step becomes the w at n-1 step
                w = w_new #The w at n+1 step becomes the w at n step

                # Vectorised receiver interpolation (replaces per-step Python loop)
                w_rec[steps] = float(np.dot(w_new[_r_idx_np], _r_w_np))

                # Vectorised source update at source tets
                if tcalc == "decay":
                    s[_s_idx_np] = source1[steps] * _s_w_np
                elif tcalc == "stationarysource":
                    s[_s_idx_np] = source1[0] * _s_w_np
            
            # Receiver-off / derivative block. The original code re-ran this
            # *every* time step but only the final iteration's values are
            # actually used (everything gets overwritten each step).
            # For recording_steps ~ 10000 and idx_w_rec slice length ~ 9000,
            # this was ~50k numpy ops per step -- the dominant CPU cost of
            # the loop, swamping the actual matvec. Guard so it only runs on
            # the last iteration; mathematically equivalent.
            if steps == recording_steps - 1:
                idx_w_rec = np.argmin(np.abs(t - sourceon_time)) #index at which the t array is equal to the sourceon_time; I want the RT to calculate from when the source stops.
                w_rec_off = w_rec[idx_w_rec:]
                p_rec_off = (w_rec[idx_w_rec:])*rho*c0**2
                t_off = t[idx_w_rec:]

                #Envelope of Impulse response from the energy density
                w_rec_off_deriv = w_rec_off #initialising an array of derivative equal to the w_rec_off -> this will be the impulse response after modifying it
                w_rec_off_deriv = np.delete(w_rec_off_deriv, 0) #delete the first element of the array -> this means shifting the array one step before and therefore do a derivation
                w_rec_off_deriv = np.append(w_rec_off_deriv,0) #add a zero in the last element of the array -> for derivation and to have the same length as previously
                #impulse = ((w_rec_off - w_rec_off_deriv))/dt#/(rho*c0**2) #This is the difference between the the energy density and the impulse response

                #Envelope of Impulse response from the pressure
                p_rec_off_transf = p_rec_off #initialising an array of derivative equal to the w_rec_off -> this will be the impulse response after modifying it
                p_rec_off_transf = np.delete(p_rec_off_transf, 0) #delete the first element of the array -> this means shifting the array one step before and therefore do a derivation
                p_rec_off_transf = np.append(p_rec_off_transf,0) #add a zero in the last element of the array -> for derivation and to have the same length as previously
                p_rec_off_deriv = ((p_rec_off - p_rec_off_transf))/dt

            #print(time_steps)
            percentDone = round(
                100 * (iBand / nBands + steps / recording_steps * 1 / nBands)
            )
            if percentDone > prev_percent_done:

                print(str(percentDone) + "% of main calculation completed")
                if called_from_choras:
                    
                    # Checking whether the user has cancelled the simulation (only one time per percentage increase)
                    if check_should_cancel(json_file):
                        print("breaking out of inner loop")
                        break

                    result_container["results"][0]["percentage"] = percentDone
                    with open(json_file, "w") as percentage_update:
                        percentage_update.write(
                            json.dumps(result_container, indent=4)
                        )

            prev_percent_done = percentDone

        if check_should_cancel(json_file):
            print("breaking out of outer loop")
            break

        # percentDone = round(100*iBand/nBands);
        #if (percentDone > curPercent):
        # print(str(percentDone) + "% of main calculation completed")
            #curPercent += 1;

        if use_gpu:
            # Bring the final per-band state back to CPU so the rest of the
            # pipeline (which is numpy-based) works unchanged. Also reflect
            # any source-mutation back into the CPU-side `s` so subsequent
            # bands start from the same state as the CPU branch would.
            w_new = w_d.get()
            s[:] = s_d.get()

        w_new_band.append(w_new)
        w_rec_band.append(w_rec)
        w_rec_off_band.append(w_rec_off)
        w_rec_off_deriv_band.append(w_rec_off_deriv)
        p_rec_off_deriv[-1] = 0 # making sure that unfinished impulse responses are not causing a click at the end because of the derivative
        p_rec_off_deriv_band.append(p_rec_off_deriv)
        
        import warnings
        warnings.filterwarnings("ignore")
    
    return w_new_band, w_rec_band, w_rec_off_band, w_rec_off_deriv_band, p_rec_off_deriv_band, idx_w_rec, t_off


#%%
###############################################################################
#POST-PROCESS
###############################################################################

def freq_parameters(nBands, line_rec_x_idx_list, w_new_band, line_rec_y_idx_list, rho, c0, Ws, dist_x, dist_y, pRef, w_rec_band, w_rec_off_band, tcalc, t, idx_w_rec, V, Eq_A, S, dist_sr):
    """
    Computation of sound presure level and reverberation time parameters

    Parameters
    ----------
        nBands : int
            Number of frequency bands
        line_rec_x_idx_list : list
            Indices of center cells close to the line receiver position in the x axis
        w_new_band : list of arrays
            Energy density at the time step n+1 at each centre cell per each frequency band
        line_rec_y_idx_list : list
            Indices of center cells close to the line receiver position in the y axis
        rho : float 
            Density of air
        c0 : int
            Speed of sound
        Ws : float
            Power of the source
        dist_x : array of floats
            Distance between each line receiver point cell and the source in the x axis
        dist_y : array of floats
            Distance between each line receiver point cell and the source in the y axis
        pRef : float
            Reference pressure
        w_rec_band : list of arrays
            Energy density over time at the receiver position per each frequency band
        w_rec_off_band : list of arrays
            Energy density over time after the source is switched off at the receiver position per each frequency band
        tcalc : str
            Type of calculation; "decay" if the source switches off and "stationarysource" if the source is stationary
        t : array of floats
            Time array of time steps
        idx_w_rec : int
            Time index at which the source is switched off
        V : float
            Volume of the room
        Eq_A : array of floats
            Equivalent absorption area
        S : float
            Total surface are of the room
        dist_sr : float
            Distance between source and receiver position

    Returns
    -------
        w_rec_x_band : list of arrays
            Energy density over time at the each line point receiver position in the x axis per each frequency band
        w_rec_y_band : list of arrays
            Energy density over time at the each line point receiver position in the y axis per each frequency band
        spl_stat_x_band : list of array
            Total sound pressure level (direct field plus reverberant) over time at the each line point receiver position in the x axis per each frequency band
        spl_stat_y_band : list of arrays
            Total sound pressure level (direct field plus reverberant) over time at the each line point receiver position in the y axis per each frequency band
        spl_r_band : list of arrays
            Sound pressure level over time at the receiver position per each frequency band
        spl_r_off_band : list of arrays
            Sound pressure level over time after the source is switched off at the receiver position per each frequency band
        spl_r_norm_band : list of arrays
            Sound pressure level over time at the receiver position per each frequency band normalised to its maximum level
        sch_db_band : list of arrays
            Energy density over time after the source is switched off at the receiver position per each frequency band
        spl_r_t0_band : list
            Energy density at the time when the source is switched off at the receiver position per each frequency band
        t30_band : array of floats
            Reverberation time T30 per each frequency band
        t20_band : array of floats
            Reverberation time T20 per each frequency band
        edt_band : array of floats
            Early decay time per each frequency band
        c80_band : array of floats 
            Clarity per each frequency band
        d50_band : array of floats 
            Definition per each frequency band
        ts_band : array of floats
            Centre time per each frequency band
    """
    w_rec_x_band = []
    w_rec_y_band = []
    spl_stat_x_band = []
    spl_stat_y_band = []
    spl_r_band = []
    spl_r_off_band = []
    spl_r_norm_band = []
    t30_band = []
    t20_band = []
    edt_band = []
    c80_band = []
    d50_band = []
    ts_band = []
    sch_db_band = []
    spl_r_t0_band = []
    
    enough_IR_for_parameters = True

    for iBand in range(nBands):
        
        w_rec_x_end = np.array([])
        for xr in line_rec_x_idx_list:
            w_rec_x = w_new_band[iBand][xr]
            w_rec_x_end = np.append(w_rec_x_end, w_rec_x)
            
        w_rec_y_end = np.array([])
        for yr in line_rec_y_idx_list:
            w_rec_y = w_new_band[iBand][yr]
            w_rec_y_end = np.append(w_rec_y_end, w_rec_y)
    
        spl_stat_x = 10*np.log10(rho*c0*(((Ws)/(4*math.pi*(dist_x**2))) + ((abs(w_rec_x_end)*c0)))/(pRef**2))
        spl_stat_y = 10*np.log10(rho*c0*(((Ws)/(4*math.pi*(dist_y**2))) + ((abs(w_rec_y_end)*c0)))/(pRef**2)) #It should be the spl stationary 
    
        #press_r = ((abs(w_rec_band[iBand]))*rho*(c0**2)) #pressure at the receiver
        spl_r = 10*np.log10(((abs(w_rec_band[iBand]))*rho*(c0**2))/(pRef**2)) #,where=press_r>0, sound pressure level at the receiver
        spl_r_off = 10*np.log10(((abs(w_rec_off_band[iBand]))*rho*(c0**2))/(pRef**2))
        
        spl_r_t0 = spl_r_off[0]
        
        spl_r_norm = 10*np.log10((((abs(w_rec_band[iBand]))*rho*(c0**2))/(pRef**2)) / np.max(((abs(w_rec_band[iBand]))*rho*(c0**2))/(pRef**2))) #normalised to maximum to 0dB

        schroeder = w_rec_off_band[iBand] #energy_r_rev_cum[::-1] #reverting the array again -> creating the schroder decay
        sch_db = 10.0 * np.log10(schroeder / max(schroeder)) #level of the array: schroeder decay
        
        if tcalc == "decay":
            try:
                t30 = t60_decay(t, sch_db, idx_w_rec, rt='t30') #called function for calculation of t60 [s]
                t20 = t60_decay(t, sch_db, idx_w_rec, rt='t20') #called function for calculation of t60 [s]
                edt = t60_decay(t, sch_db, idx_w_rec, rt='edt') #called function for calculation of edt [s]
                #Eq_A = 0.16*V/t60 #equivalent absorption area defined from the RT 
                c80 = clarity(t30, V, Eq_A[iBand], S, c0, dist_sr) #called function for calculation of c80 [dB]
                d50 = definition(t30, V, Eq_A[iBand], S, c0, dist_sr) #called function for calculation of d50 [%]
                ts = centretime(t30, Eq_A[iBand], S) #called function for calculation of ts [ms]
                
            except Exception as ex:
                t30 = -1
                t20 = -1
                edt = -1
                c80 = -1
                d50 = -1
                ts = -1
                print("Not enough impulse response data to caluclate parameters.")

            t30_band.append(t30)
            t20_band.append(t20)
            edt_band.append(edt)
            c80_band.append(c80)
            d50_band.append(d50)
            ts_band.append(ts)

        w_rec_x_band.append(w_rec_x_end)
        w_rec_y_band.append(w_rec_y_end)
        spl_stat_x_band.append(spl_stat_x)
        spl_stat_y_band.append(spl_stat_y)
        spl_r_band.append(spl_r)
        spl_r_off_band.append(spl_r_off)
        spl_r_norm_band.append(spl_r_norm)
        sch_db_band.append(sch_db)
        spl_r_t0_band.append(spl_r_t0)
    
    spl_r_off_band = np.array(spl_r_off_band)
    t30_band = np.array(t30_band)
    t20_band = np.array(t20_band)
    edt_band = np.array(edt_band)
    c80_band = np.array(c80_band)
    d50_band = np.array(d50_band)
    ts_band = np.array(ts_band)
    
    return w_rec_x_band, w_rec_y_band, spl_stat_x_band, spl_stat_y_band, spl_r_band, spl_r_off_band, spl_r_norm_band, sch_db_band, spl_r_t0_band, t30_band, t20_band, edt_band, c80_band, d50_band, ts_band

#%%
###############################################################################
#SAVING
###############################################################################

# Save all variables to a file
def save_fvm(filename,variables):
    """
    Saving of variables

    Parameters
    ----------
        filename : str
            Name of the file to save the results
        variables : dict
            Compilation of all the variables of the overall simulation
    """
    with open(filename, 'wb') as f:
        # Filter out modules, functions, and other unsupported types
        filtered_variables = {}
        for k, v in variables.items():
            try:
                # Check if the object can be pickled
                pickle.dumps(v)
                # Exclude some types explicitly known to cause issues
                if not k.startswith('__') and not isinstance(v, (types.ModuleType, types.FunctionType, types.BuiltinFunctionType, types.LambdaType, types.MethodType, types.MappingProxyType)):
                    filtered_variables[k] = v
            except Exception as e:
                print(f"Could not pickle {k}: {str(e)}")

        pickle.dump(filtered_variables, f)


if __name__ == '__main__':
    st = time.time() #start time of calculation

