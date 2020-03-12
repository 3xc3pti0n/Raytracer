# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:25:47 2019

@author: Julian
"""

import subprocess as sp
import shutil
import random
from pathlib import Path
import os
import Algorithmen
import numpy as np
import stl
import meshlabxml as mlx

#add meshlabserver to PATH
meshlabserver_path = 'C:\\Program Files\\VCG\\MeshLab'
os.environ['PATH'] = meshlabserver_path + os.pathsep + os.environ['PATH']

#Set mesh folders
new_mesh_folder = ("D:\Raytracer\\Remesh")
pointcloud_folder = ("D:\Raytracer\Pointclounds")
script_folder = ("D:\Raytracer\Meshlab_scripte")

#Set Raytracing folders
Raytracer_Path = 'D:\Raytracer'
#Raytracer_exe = os.path.join(Raytracer_Path,'The_Raytracer.jar')
Raytracer_exe = os.path.join(Raytracer_Path,'Raytracer_test_new.jar')
path_input_data = ("D:\Raytracer\Raytracer_local_input")
path_output_data = ("D:\Raytracer\Raytracer_local_output")
File_name = '200x200mm_20402tris_remesh_remesh'

n_pulses = 10
for i in range(1, n_pulses+1):
    print (i)
    # STL zu TRI konvertieren
    print('STL konvertieren')
    stl_file = os.path.join(path_input_data,File_name+".stl")
    tri_file = os.path.join(path_input_data,File_name+".tri")
    Algorithmen.stl_to_tri(stl_file, tri_file)
    #Algorithmen.STL2tri(stl_file, tri_file)
    
    
    #Raytracing calculation
    print('tracen')
    Parameterfile = os.path.join(Raytracer_Path, 'Parameter_Raytracer_neu_pulsed.txt')
    OUTfile = os.path.join(path_output_data,File_name)
    raytracer_arguments = [OUTfile, Parameterfile, tri_file] 
    cmd = "java -jar "+ Raytracer_exe +' '+ " ".join(raytracer_arguments)
    sp.call(cmd)
    print('finished tracing')
    
    #Excel File auslesen 
    f = open(os.path.join(path_output_data,File_name +'_output.csv'), 'r') # open the file for reading
    data_test = []
    for row_num, line in enumerate(f):
        # Remove the new line at the end and then split the string based on
        # tabs. This creates a python list of the values.
        values = line.strip().split('\t')
        if row_num == 0: # first line is the header
             header = values
        else:
            data_test.append([float(v) for v in values])
    basic_data = np.array(data_test)
    f.close()
    
    #Excel File auslesen nummer 2
    f = open(os.path.join(path_output_data,File_name +'_output.csv'), 'r') # open the file for reading
    data_test2 = []
    for row_num, line in enumerate(f):
        # Remove the new line at the end and then split the string based on
        # tabs. This creates a python list of the values.
        values = line.strip().split('\t')
        if row_num == 0: # first line is the header
             header = values
        else:
            data_test2.append([float(v) for v in values])
    basic_data2 = np.array(data_test2)
    f.close()
    
    #datahandler = open(os.path.join(path_output_data,File_name +'_out.csv'))
    #data = np.asanyarray([[float(x) for x in l.rstrip('\n').split('\t')] for l in datahandler.readlines()[1:]])
    SPs = basic_data[:,1:4] /1e3 # Schwerpunkte in cm
    n_points = len(SPs)
    Is = basic_data[:,10] #Intensitäten
    Nvecs = basic_data[:,4:7]
    Kvecs = basic_data[:,7:10] 
    
    
    #Parameter für Abtrag
    Phi_th = 0.0283 # Schwellfluenz des Materials in J/cm^2
    delta_opt = 6e-6 # optische eindringtiefe in um
    tau_pulse =8e-12 # Pulsdauer in s
    
    
    #Verschiebung der SP entlang NV und KV
    #k-Vektoren normieren
    norm_factor = np.linalg.norm(Kvecs,axis=1)
    rec_norm_k = np.where(norm_factor!=0,np.reciprocal(norm_factor),0.)
    Norm_k = Kvecs*rec_norm_k[:, np.newaxis]
    #Verschiebungsvektoren normieren
    shift_vector = (Nvecs+Norm_k)  #Verschiebungsvektor 
    norm_factor_shift = np.linalg.norm(shift_vector,axis=1)
    rec_norm_shift = np.where(norm_factor_shift!=0,np.reciprocal(norm_factor_shift),0.)
    shift_vector_norm = shift_vector*rec_norm_shift[:, np.newaxis]
    
    
    ablation_dist = Algorithmen.calc_ablat_dist(Is,Phi_th,tau_pulse,delta_opt)
    new_SPs = SPs - shift_vector_norm*ablation_dist[:, np.newaxis]
    #new_SPs = SPs  #testen ohne verschiebung
    comb_point_data = np.concatenate((new_SPs,shift_vector_norm),1) 
    #comb_point_data = np.concatenate((new_SPs,Nvecs),1) #testen ohne verschiebung
    
    
    
    #Neue Punktwolke in ply Datei schreiben
    print('Generate new Pointcloud file')
    #pointcloud_filename = 'Pointcloud_test' os.path.join(File_name + '_pointcloud')
    pointcloud_filename = os.path.join(File_name + '_pointcloud')
    pointcloud_file = open(os.path.join(pointcloud_filename + '.ply'),'w')
    #Header der Pointcloud Datei
    pointcloud_file.write("ply \n"
    "format ascii 1.0 \n"
    "comment VCGLIB generated \n"
    "element vertex " + str(n_points) + "\n"
    "property float x \n"
    "property float y \n"
    "property float z \n"
    "property float nx \n"
    "property float ny \n"
    "property float nz \n"
    "element face 0 \n"
    "property list uchar int vertex_indices \n"
    "end_header \n")
    for line in comb_point_data:
        pointcloud_file.write((str(line[0]) + ' '+str(line[1])+' '+ str(line[2])+ ' '+ str(line[3]) + ' '+str(line[4])+' '+ str(line[5])+'\n'))
    pointcloud_file.close()
    #"\nelement vertex" + str(len(normal_stl_file_new))
    
    #Meshlab zum neue Oberfläche erzeugen
    #https://github.com/3DLIRIOUS/MeshLabXML Meshlab XML package
    print('Mashlab remashing')
    meshlabserver_path = 'C:\\Program Files\\VCG\\MeshLab'
    os.environ['PATH'] = meshlabserver_path + os.pathsep + os.environ['PATH']
    #remesh_filename = 'Pointcloud_remesh' 
    remesh_filename = os.path.join(File_name + '_remesh') 
    
    #Pointcloud_test_out = mlx.FilterScript(file_in = 'Pointcloud_test.ply', file_out = 'Pointcloud_remesh.stl', ml_version='2016.12')
    Pointcloud_test_out = mlx.FilterScript(file_in = os.path.join(pointcloud_filename + '.ply'), file_out = os.path.join(remesh_filename + '.stl'), ml_version='2016.12')
    mlx.remesh.surface_poisson_screened(Pointcloud_test_out,depth=10,
                                 full_depth=5, cg_depth=0, scale=1,
                                 samples_per_node=1.5, point_weight=4.0,
                                 iterations=1, confidence=False, pre_clean=False)
    mlx.layers.delete(Pointcloud_test_out)
    mlx.remesh.simplify(Pointcloud_test_out, texture=False, faces=150000, target_perc=0.0,
                  quality_thr=0.95, preserve_boundary=False, boundary_weight=1.0,
                  optimal_placement=True, preserve_normal=True,
                  planar_quadric=False, selected=False, extra_tex_coord_weight=1.0,              
                  preserve_topology=True, quality_weight=False, autoclean=True)
    Pointcloud_test_out.run_script()
    
    #convert to ascii encoding
    print('converting...')
    i = stl.mesh.Mesh.from_file(os.path.join(Raytracer_Path,remesh_filename + '.stl'))
    i.save(os.path.join(Raytracer_Path,remesh_filename + '.stl'), mode=stl.Mode.ASCII)
    
    #Pointcloud file verschieben
    print('move new files to folders')
    shutil.move(os.path.join(Raytracer_Path,pointcloud_filename + '.ply'), pointcloud_folder)
    shutil.move(os.path.join(Raytracer_Path,remesh_filename + '.stl'), new_mesh_folder)



