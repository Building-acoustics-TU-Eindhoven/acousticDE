# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 15:23:01 2025

@author: 20225533
"""

import math

# %%
###############################################################################
# INPUT VARIABLES
###############################################################################
# The variables to be assigned by the user directly
# Only name of the file
name_of_obj_file = '3x3x3.obj'
name_of_geo_file = "3x3x3.geo"

###############################################################################

def obj_to_gmsh_geo(obj_file, geo_file, volume_name="RoomVolume", tol=1e-8):
    """
    Create the .geo file from an .obj file.

    Parameters
    ----------
    obj_file : str
        Name of the .obj file
    geo_file : str
        Name of the .geo file
    volume_name : str
        Name of the volume
    tol: float
        Tollerance in the defintinion of vertices

    Returns
    -------
    geo_file : .geo
        .geo file for the running of the code.
    """
    # Read OBJ file
    vertices = []          # initialization vertex coordinates
    faces = []             # initialization list of faces: list of vertex indices
    face_groups = []       # initialization group name for each face
    current_group = "default"

    with open(obj_file, "r") as f: #read the file
        for raw in f:
            line = raw.strip() #read line by line
            if not line:
                continue
            if line.startswith('v '): #checks for vertices
                parts = line.split()
                x, y, z = map(float, parts[1:4])
                # Convert from Blender (left-handed) to Gmsh (right-handed)
                # Flip the Z axis to correct mirroring
                coords = (x, -z, y)
                vertices.append(coords)
            elif line.startswith('g '): #checks for groups
                parts = line.split()[1:]  # skip 'g'
                # If first token starts with 'Mesh', skip it
                parts = [p for p in parts if not p.startswith("Mesh") and not p.startswith("Model")]
                current_group = parts[0] if parts else "default"
            elif line.startswith('usemtl'):  # Material-based grouping
                parts = line.split()[1:]
                # Materials also act as groups
                current_group = parts[0] if parts else "default"
            elif line.startswith('f '): #checks for faces
                parts = line.split()[1:]
                # face vertex indices (OBJ format may include v/vt/vn)
                idxs = [int(p.split('/')[0]) for p in parts]
                faces.append(idxs)
                face_groups.append(current_group)

    # Create duplicates of vertices to be assign to each face, preserving the order
    unique_vertices = []
    orig_to_unique = {}  # map from original index (1-based) -> unique index (1-based)
    for i, v in enumerate(vertices, start=1):
        # find existing within tolerance
        found = None
        for j, uv in enumerate(unique_vertices, start=1):
            if abs(uv[0] - v[0]) < tol and abs(uv[1] - v[1]) < tol and abs(uv[2] - v[2]) < tol:
                found = j
                break
        if found is None:
            unique_vertices.append(v)
            orig_to_unique[i] = len(unique_vertices)
        else:
            orig_to_unique[i] = found
            
    # Remap faces to unique vertex indices
    faces_mapped = [[orig_to_unique[i] for i in face] for face in faces]
    
    # Sort vertices deterministically 
    unique_vertices_sorted = sorted(
        enumerate(unique_vertices, start=1),
        key=lambda kv: (round(kv[1][0], 8), 
                        round(kv[1][1], 8),
                        round(kv[1][2], 8)) 
        ) 
    index_map = {old: new for new, (old, _) in enumerate(unique_vertices_sorted, start=1)} 
    unique_vertices = [v for _, v in unique_vertices_sorted] 
    faces_mapped = [[index_map[i] for i in face] for face in faces_mapped]


    # Merge pairs of triangles within same group into quads when needed
    nfaces = len(faces_mapped)
    merged_flag = [False] * nfaces
    merged_faces = []      # list of faces (each is 3 or 4 vertex indices)
    merged_groups = []

    # Get coords by unique index (1-based)
    def coords(idx):
        return unique_vertices[idx - 1]

    for i in range(nfaces):
        if merged_flag[i]:
            continue
        fi = faces_mapped[i]
        gi = face_groups[i]
        if len(fi) == 3:
            # try to find a partner triangle in same group sharing 2 vertices
            partner = None
            for j in range(i + 1, nfaces):
                if merged_flag[j]:
                    continue
                if face_groups[j] != gi:
                    continue
                fj = faces_mapped[j]
                if len(fj) != 3:
                    continue
                shared = set(fi) & set(fj)
                if len(shared) == 2:
                    partner = j
                    break
            if partner is not None:
                # build quad from union of vertices (4 vertices)
                union = list(dict.fromkeys(fi + faces_mapped[partner]))  # preserve order somewhat
                if len(union) == 4:
                    # order the 4 vertices into a planar loop consistently
                    pts = [coords(idx) for idx in union]
                    # compute plane normal using first triangle
                    v0 = pts[0]
                    v1 = pts[1]
                    v2 = pts[2]
                    nx = (v1[1] - v0[1]) * (v2[2] - v0[2]) - (v1[2] - v0[2]) * (v2[1] - v0[1])
                    ny = (v1[2] - v0[2]) * (v2[0] - v0[0]) - (v1[0] - v0[0]) * (v2[2] - v0[2])
                    nz = (v1[0] - v0[0]) * (v2[1] - v0[1]) - (v1[1] - v0[1]) * (v2[0] - v0[0])
                    an = (abs(nx), abs(ny), abs(nz))
                    # choose projection plane by largest normal component
                    if an[2] >= an[0] and an[2] >= an[1]:
                        # project to XY
                        proj = lambda p: (p[0], p[1])
                    elif an[1] >= an[0] and an[1] >= an[2]:
                        # project to XZ
                        proj = lambda p: (p[0], p[2])
                    else:
                        # project to YZ
                        proj = lambda p: (p[1], p[2])

                    uv = [proj(coords(idx)) for idx in union]
                    cx = sum(pt[0] for pt in uv) / 4.0
                    cy = sum(pt[1] for pt in uv) / 4.0
                    angles = [math.atan2(pt[1] - cy, pt[0] - cx) for pt in uv]
                    # sort union vertices by angle
                    union_ordered = [x for _, x in sorted(zip(angles, union))]
                    merged_faces.append(union_ordered)
                    merged_groups.append(gi)
                    merged_flag[i] = True
                    merged_flag[partner] = True
                    continue
                # if union not 4, fallthrough to keep triangle
            # no partner found => keep triangle
            merged_faces.append(fi)
            merged_groups.append(gi)
            merged_flag[i] = True
        else:
            # non-triangle face: keep as-is (maybe quad)
            merged_faces.append(fi)
            merged_groups.append(gi)
            merged_flag[i] = True

    # There may be faces leftover (if any not processed): ensure all covered
    for k in range(nfaces):
        if not merged_flag[k]:
            merged_faces.append(faces_mapped[k])
            merged_groups.append(face_groups[k])
            
    tag_to_surfaces = {}
    for sid, tag in enumerate(merged_groups, start=1):
        tag_to_surfaces.setdefault(tag, []).append(sid)

    room_center = tuple(
        sum(v[i] for v in unique_vertices) / len(unique_vertices)
        for i in range(3)
    )

    
    def is_outward_facing(face):
        """
        Checks if the normal of the face is facing outwards.

        Parameters
        ----------
        face : list
            List of vertices forming the face

        Returns
        -------
        dot : float
            If the float is negative, the normal points outward.
        """
        pts = [unique_vertices[i - 1] for i in face]
        v1, v2, v3 = pts[:3]
        # Face normal
        nx = (v2[1]-v1[1])*(v3[2]-v1[2]) - (v2[2]-v1[2])*(v3[1]-v1[1])
        ny = (v2[2]-v1[2])*(v3[0]-v1[0]) - (v2[0]-v1[0])*(v3[2]-v1[2])
        nz = (v2[0]-v1[0])*(v3[1]-v1[1]) - (v2[1]-v1[1])*(v3[0]-v1[0])
        normal = (nx, ny, nz)
        # Face centroid
        cx, cy, cz = tuple(sum(p[i] for p in pts) / len(pts) for i in range(3))
        # Vector from centroid to room center
        to_center = (
            room_center[0] - cx,
            room_center[1] - cy,
            room_center[2] - cz
        )
        # Dot product: if negative, normal points outward
        dot = sum(normal[i] * to_center[i] for i in range(3))
        return dot < 0
    
    for face in merged_faces:
        if not is_outward_facing(face):
            face.reverse()



    # Build unique edges (lines) with stable orientation
    edge_to_line = {}       # key = (min,max) -> line_id
    line_orientation = {}   # line_id -> (a,b) orientation used when created
    next_line_id = 1

    # collect edges from merged faces in consistent order
    face_line_loops = []  # list of lists of signed line indices (to write)
    for face in merged_faces:
        n = len(face)
        loop_line_ids = []
        for idx in range(n):
            a = face[idx]
            b = face[(idx + 1) % n]
            key = (a, b) if a < b else (b, a)
            if key not in edge_to_line:
                edge_to_line[key] = next_line_id
                # store orientation as the first encountered direction (a,b)
                if key == (a, b):
                    line_orientation[next_line_id] = (a, b)
                else:
                    line_orientation[next_line_id] = (b, a)
                next_line_id += 1
            lid = edge_to_line[key]
            # determine sign: +if orientation matches (a,b), - otherwise
            ori = line_orientation[lid]
            if ori == (a, b):
                loop_line_ids.append(lid)
            else:
                loop_line_ids.append(-lid)
        face_line_loops.append(loop_line_ids)

    # Write the GEO file
    with open(geo_file, "w") as g:
        # Points
        for i, v in enumerate(unique_vertices, start=1):
            g.write(f"Point({i}) = {{ {v[0]}, {v[1]}, {v[2]}, 1.0 }};\n")
        g.write("\n")

        # Lines (must write using stored orientation endpoints)
        # We need to output unique edge list, using the stored orientation endpoints
        # Build a mapping of line_id -> endpoints
        line_id_to_endpoints = {}
        for key, lid in edge_to_line.items():
            # endpoints should be line_orientation[lid]
            a, b = line_orientation[lid]
            line_id_to_endpoints[lid] = (a, b)

        # Write lines in increasing id order
        for lid in range(1, next_line_id):
            a, b = line_id_to_endpoints[lid]
            g.write(f"Line({lid}) = {{ {a}, {b} }};\n")
        g.write("\n")

        # Line Loops and Plane Surfaces
        for sid, loop in enumerate(face_line_loops, start=1):
            loop_str = ", ".join(str(x) for x in loop)
            g.write(f"Line Loop({sid}) = {{ {loop_str} }};\n")
            g.write(f"Plane Surface({sid}) = {{ {sid} }};\n")
        g.write("\n")

        # Physical Line (all lines)
        lines_all = ", ".join(str(i) for i in range(1, next_line_id))
        g.write(f'Physical Line("default") = {{ {lines_all} }};\n')

        # Physical Surface groups by OBJ group name
        unique_groups = []
        for grp in merged_groups:
            if grp not in unique_groups:
                unique_groups.append(grp)
        for grp in unique_groups:
            surf_ids = [str(i + 1) for i, gname in enumerate(merged_groups) if gname == grp]
            if surf_ids:
                g.write(f'Physical Surface("{grp}") = {{ {", ".join(surf_ids)} }};\n')

        # Surface Loop and Volume
        total_surfaces = len(face_line_loops)
        surf_list = ", ".join(str(i) for i in range(1, total_surfaces + 1))
        g.write(f"Surface Loop(1) = {{ {surf_list} }};\n")
        g.write("Volume(1) = { 1 };\n")
        g.write(f'Physical Volume("{volume_name}") = {{ 1 }};\n')
        g.write('Mesh.Algorithm = 6;\n')
        g.write('Mesh.Algorithm3D = 1; // Delaunay3D, works for boundary layer insertion.\n')
        g.write('Mesh.Optimize = 1; // Gmsh smoother, works with boundary layers (netgen version does not).\n')
        g.write('Mesh.CharacteristicLengthFromPoints = 1;\n')
        g.write('// Recombine Surface "*";\n')

    print(f"Wrote {geo_file}: {len(unique_vertices)} points, {next_line_id-1} lines, {len(face_line_loops)} surfaces.")



# Example
obj_to_gmsh_geo(name_of_obj_file, name_of_geo_file)


