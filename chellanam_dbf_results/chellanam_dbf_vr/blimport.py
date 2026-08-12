import bpy
import csv
import math
import os
import shutil

# **********************************************************************

# These parameters are automatically set by r.avaflow

pf = 'chellanam_dbf_' # prefix for objects
xdim = 168
ydim = 121
xmin = 640304.688
xmax = 648132.625
ymin = 1081277.125
ymax = 1086915.125
nfiles = 31
undef = -9999.000
layers = 0

# **********************************************************************

# These parameters have to be adapted manually, if required

wkdir = '' # location of data directory in file system (full path with normal slashes, trailing slash required)

impfiles = 1 # control for csv file import (1 = activated, 0 = deactivated)
standard = 0 # control for standard video (1 = activated, 0 = deactivated)
anaglyph = 0 # control for anaglyph video (1 = activated, 0 = deactivated)
stereo = 0 # control for stereo 3D video (1 = activated, 0 = deactivated)
animation = 0 # control for animation

vbackground = -9999.0 # background surface elevation value (excluded from meshes)
include_bg = 2 # control for including background (0 = no, 1 = only for impact area, 2 = yes)
subsurface = 0 # control for including phase 1 surface (only useful for layered multi-phase model: 0 = no, 1 = yes)
successive = 0 # control for successive simulation (0 = no, 1 = yes)

light = bpy.data.objects['Light'] # name of light object
light.location = (3, 4, 5) # light location
light.data.energy=200.0 # light energy

cam = bpy.data.objects['Camera'] # name of camera object
camx = 0.0 # x position of camera
camy = 0.0 # y position of camera
camz = 8.0 # z position of camera
cam.location = (camx, camy, camz) # camera location
cam.rotation_euler = (0.0 * math.pi/180, 0.0 * math.pi/180, 0.0 * math.pi/180) # camera rotation

camspeedx = 0.0 # camera speed in x direction during video (km per frame)
camspeedy = 0.0 # camera speed in y direction during video (km per frame)
camspeedz = 0.0 # camera speed in z direction during video (km per frame)
convdist = 3 # convergence distance of stereo camera
resx = 11 # control for x dimension of stereo 3D video (bits)
interval = 3 # interval between animation time steps (frames)

# **********************************************************************

# The following code should only be modified by experienced users

prefix = 'pv'
zref = 0.0

scn = bpy.context.scene
scn.camera = cam
scn.render.film_transparent=False
scn.render.image_settings.file_format='PNG'

collection = bpy.context.collection
meshes = set()

bpy.types.View3DShading.color = 'ATTRIBUTE'

if standard == 1 and os.path.exists(wkdir + pf + 'standard') and os.path.isdir(wkdir + pf + 'standard'):
    shutil.rmtree(wkdir + pf + 'standard')
    os.makedirs(wkdir + pf + 'standard')

if anaglyph == 1 and os.path.exists(wkdir + pf + 'anaglyph') and os.path.isdir(wkdir + pf + 'anaglyph'):
    shutil.rmtree(wkdir + pf + 'anaglyph')
    os.makedirs(wkdir + pf + 'anaglyph')

if stereo == 1 and os.path.exists(wkdir + pf + 'stereo') and os.path.isdir(wkdir + pf + 'stereo'):
    shutil.rmtree(wkdir + pf + 'stereo')
    os.makedirs(wkdir + pf + 'stereo')

ob = []
if layers > 0 and subsurface == 1: ob1 = []

if impfiles == 1:

    if successive == 0:

        for obj in [o for o in collection.objects if o.type == 'MESH']:

            meshes.add( obj.data )
            bpy.data.objects.remove( obj )

        for mesh in [m for m in meshes if m.users == 0]:
            bpy.data.meshes.remove( mesh )

        for material in bpy.data.materials: bpy.data.materials.remove(material)

        newmat = bpy.data.materials.new('VertCol')
        newmat.use_nodes = True
        
    else:
        
        newmat = bpy.data.materials['VertCol']

    if include_bg == 1:

        if nfiles < 10: filetext = '000' + str(nfiles)
        elif nfiles < 100: filetext = '00' + str(nfiles)
        elif nfiles < 1000: filetext = '0' + str(nfiles)
        else: filetext = str(file)

        path = wkdir + 'data/' + prefix + filetext + '.csv'

        imp = []

        with open(path, 'r') as csvfile:
            datareader = csv.reader(csvfile)
            next(datareader)

            i = 0

            for row in datareader:
        
                imp.append(int(row[8]))
           
                i += 1

    for file in range(0, nfiles + 1):
    
        if file < 10: filetext = '000' + str(file)
        elif file < 100: filetext = '00' + str(file)
        elif file < 1000: filetext = '0' + str(file)
        else: filetext = str(file)

        path = wkdir + 'data/' + prefix + filetext + '.csv'

        vertices = []
        edges = []
        faces = []
        colors = []

        if subsurface == 1 and file == 0: verticesb = []

        if layers > 0 and subsurface == 1:

            vertices1 = []
            edges1 = []
            faces1 = []
            colors1 = []

        if not include_bg == 1: imp = []

        xcotr = 0
        ycotr = 0

        with open(path, 'r') as csvfile:
            datareader = csv.reader(csvfile)
            next(datareader)

            i = 0

            for row in datareader:

                x = float(row[0]) * 0.001 - ( xmax + xmin ) * 0.5 * 0.001
                y = float(row[1]) * 0.001 - ( ymax + ymin ) * 0.5 * 0.001
                z0 = float(row[2])

                if z0 == vbackground: z = undef
                else: z = z0 * 0.001 - zref * 0.001

                if subsurface == 1 and file == 0:

                    zb0 = float(row[3])

                    if zb0 == vbackground: zb = undef
                    else: zb = zb0 * 0.001 - zref * 0.001

                if layers > 0 and subsurface == 1:

                    z10 = float(row[4])

                    if z10 == vbackground: z1 = undef
                    else: z1 = z10 * 0.001 - zref * 0.001

                r = float(row[5])
                g = float(row[6])
                b = float(row[7])
                if include_bg == 0: imp.append(int(row[8]))
                elif include_bg == 2: imp.append(1)
           
                vertices.append((x, y, z))
                if subsurface == 1 and file == 0: verticesb.append((x, y, zb))
                if layers > 0 and subsurface == 1: vertices1.append((x, y, z1))

                if not z0 == undef and xcotr > 0 and ycotr > 0:

                    i1 = i-1
                    i2 = i
                    i3 = i-xdim
                    i4 = i-xdim-1

                    if ( not vertices[i1][2] == undef and not vertices[i2][2] == undef 
                        and not vertices[i3][2] == undef and not vertices[i4][2] == undef and not imp[i] == 0 ):

                        faces.append((i1, i2, i3, i4))
                        colors.append((r, g, b, 1))

                    if ( layers > 0 and subsurface == 1 and not vertices1[i1][2] == undef and not vertices1[i2][2] == undef 
                        and not vertices1[i3][2] == undef and not vertices1[i4][2] == undef and not imp[i] == 0 
                        and ( vertices1[i1][2] < vertices[i1][2] or vertices1[i2][2] < vertices[i2][2] 
                        or vertices1[i3][2] < vertices[i3][2] or vertices1[i4][2] < vertices[i4][2] )):
                
                        faces1.append((i1, i2, i3, i4))

                if xcotr < xdim-1:
                    xcotr += 1
                else:
                    xcotr = 0
                    ycotr += 1

                i += 1

        mesh = bpy.data.meshes.new(pf + 'Mesh' + str(filetext))
        mesh.from_pydata(vertices, edges, faces)
        mesh.update()
        ob.append(bpy.data.objects.new(pf + 'Surface' + str(filetext), mesh))
        ob[file].data.update()

        if subsurface == 1 and file == 0:

            meshb = bpy.data.meshes.new(pf + 'Meshb')
            meshb.from_pydata(verticesb, edges, faces)
            meshb.update()
            obb = bpy.data.objects.new(pf + 'Surfaceb', meshb)
            obb.data.update()
            bpy.context.collection.objects.link(obb)

        if layers > 0 and subsurface == 1:

            mesh1 = bpy.data.meshes.new(pf + 'Mesh1' + str(filetext))
            mesh1.from_pydata(vertices1, edges1, faces1)
            mesh1.update()
            ob1.append(bpy.data.objects.new(pf + 'Surface1' + str(filetext), mesh1))
            ob1[file].data.update()

        if not mesh.vertex_colors:
            mesh.vertex_colors.new()

        color_layer = mesh.vertex_colors.active

        j=0
        k=0
        for face in mesh.polygons:

            for idx in face.loop_indices:
                color_layer.data[k].color = colors[j]
                k += 1

            j += 1

        if mesh.materials:
            mesh.materials[0] = newmat
        else:
            mesh.materials.append(newmat) 

        nodes = newmat.node_tree.nodes
        bsdf = nodes.get('Principled BSDF') 

        vertex_color_node = None
        if not 'Vertex Color' in [node.type for node in nodes]:
            vertex_color_node = nodes.new(type = 'ShaderNodeVertexColor')
        else:
            vertex_color_node = nodes.get('Vertex Color')

        vertex_color_node.layer_name = 'Col'

        links = newmat.node_tree.links
        link = links.new(vertex_color_node.outputs[0], bsdf.inputs[0])

        bpy.context.collection.objects.link(ob[file])
        if layers > 0 and subsurface == 1: bpy.context.collection.objects.link(ob1[file])

if standard == 1 or anaglyph == 1 or stereo == 1:

    if impfiles == 0:
        
        for file in range(0, nfiles + 1):

            if file < 10: filetext = '000' + str(file)
            elif file < 100: filetext = '00' + str(file)
            elif file < 1000: filetext = '0' + str(file)
            else: filetext = str(file)
            
            ob.append(bpy.data.objects[pf + 'Surface' + str(filetext)])
            ob[file].hide_render = True

    for file in range(0, nfiles + 1):

        if file < 10: filetext = '000' + str(file)
        elif file < 100: filetext = '00' + str(file)
        elif file < 1000: filetext = '0' + str(file)
        else: filetext = str(file)

        cam.location = (camx + float(file) * camspeedx, 
            camy + float(file) * camspeedy, 
            camz + float(file) * camspeedz)
    
        ob[file].hide_render = False

        if standard == 1:

            scn.render.use_multiview = False
            scn.render.resolution_x = 1920
            scn.render.resolution_y = 1080
            
            scn.render.filepath = wkdir + pf + 'standard/standard' + filetext + '.png'
            bpy.ops.render.render(write_still=1)

        if anaglyph == 1:

            scn.render.image_settings.views_format = 'STEREO_3D'
            scn.render.image_settings.stereo_3d_format.display_mode = 'ANAGLYPH'
            scn.render.use_multiview = True
            scn.render.resolution_x = 1920
            scn.render.resolution_y = 1080
            cam.data.stereo.convergence_distance = convdist
            cam.data.stereo.interocular_distance = convdist / 30.0
    
            scn.render.filepath = wkdir + pf + 'anaglyph/anaglyph' + filetext + '.png'
            bpy.ops.render.render(write_still=1)

        if stereo == 1:

            scn.render.image_settings.views_format = 'STEREO_3D'
            scn.render.image_settings.stereo_3d_format.display_mode = 'SIDEBYSIDE'
            scn.render.use_multiview = True
            scn.render.resolution_x = int(pow(2, resx))
            scn.render.resolution_y = int(pow(2, resx)) #int(scn.render.resolution_x / 2)
            cam.data.type = 'PANO'
            cam.data.stereo.use_spherical_stereo = True
            cam.data.stereo.convergence_mode = 'OFFAXIS'
            cam.data.stereo.convergence_distance = convdist
            cam.data.stereo.interocular_distance = convdist / 30.0
    
            scn.render.filepath = wkdir + pf + 'stereo/stereo' + filetext + '.png'
            bpy.ops.render.render(write_still=1)

        ob[file].hide_render = True

    for file in range(0, nfiles + 1):

        ob[file].hide_render = False

if animation == 1:

    arm = bpy.data.objects['Armature']
    arm.hide_render = True
    arm.hide_viewport = False

    for action in bpy.data.actions: bpy.data.actions.remove(action)
    
    bpy.context.view_layer.objects.active = arm
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    
    edit_bones = arm.data.edit_bones
    for bone in edit_bones:
        edit_bones.remove(bone)
    
    armt = []

    for file in range(0, nfiles + 1):

        armt.append(edit_bones.new('Bone' + str(file)))
        armt[file].head = (0.0, -1.0, 0.0)
        armt[file].tail = (0.0, 0.0, 0.0)

    bpy.ops.object.mode_set(mode='OBJECT')

    arm.animation_data_create()
    arm.animation_data.action = bpy.data.actions.new(name='AnimationAction')

    current_frame = 0
    pos = []
    scn.frame_start = 0
    scn.frame_end = ( nfiles + 1 ) * interval - 1

    for file in range(0, nfiles + 1):

        if file < 10: filetext = '000' + str(file)
        elif file < 100: filetext = '00' + str(file)
        elif file < 1000: filetext = '0' + str(file)
        else: filetext = str(file)

        if impfiles == 0: ob.append(bpy.data.objects[pf + 'Surface' + str(filetext)])

        ob[file].parent=arm
        ob[file].parent_type = 'BONE'
        ob[file].parent_bone = 'Bone' + str(file)

        if layers > 0 and subsurface == 1:

            ob1[file].parent=arm
            ob1[file].parent_type = 'BONE'
            ob1[file].parent_bone = 'Bone' + str(file)

        armf = arm.pose.bones['Bone' + str(file)]

        pos.append(armf.location)
        pos[file][0] += 1000
        pos[file][1] += 1000
        pos[file][2] += 1000
        armf.location = pos[file]
        ob[file].hide_render = False

    for file in range(0, nfiles + 1):

        armf = arm.pose.bones['Bone' + str(file)]

        fcurve = arm.animation_data.action.fcurves.new(data_path='location', index=file)
        k1 = fcurve.keyframe_points.insert(frame=current_frame, value=0)
        k1.interpolation = 'CONSTANT'
        
        if file > 0: armf.keyframe_insert(data_path = 'location', frame = current_frame - 1)

        pos[file][0] -= 1000
        pos[file][1] -= 1000
        pos[file][2] -= 1000
        armf.location = pos[file]
        armf.keyframe_insert(data_path = 'location', frame = current_frame)
        
        if file < nfiles: 
            
            armf.keyframe_insert(data_path = 'location', frame = current_frame + interval-1)
        
            pos[file][0] += 1000
            pos[file][1] += 1000
            pos[file][2] += 1000
            armf.location = pos[file]
            armf.keyframe_insert(data_path = 'location', frame = current_frame + interval)
            
        current_frame += interval
    
    bpy.ops.nla.bake(frame_start=0, frame_end=current_frame-interval, step=1, only_selected=False, 
        visual_keying=False, clear_constraints=False, clear_parents=False, 
        use_current_action=False, clean_curves=False, bake_types={'POSE'})
