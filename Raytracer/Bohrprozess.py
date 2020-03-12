# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:25:47 2019

@author: Julian
"""

import subprocess as sp
import shutil
import random
import time
from pathlib import Path
import os
import Algorithmen
import numpy as np
import stl
import meshlabxml as mlx
import openmesh

#add meshlabserver to PATH
meshlabserver_path = 'C:\\Program Files\\VCG\\MeshLab'
os.environ['PATH'] = meshlabserver_path + os.pathsep + os.environ['PATH']

#Set mesh folders
new_mesh_folder = ("D:\Raytracer_data\Remesh")
pointcloud_folder = ("D:\Raytracer_data\Pointclounds")

#Set Raytracing folders
Raytracer_Path = 'D:\Raytracer'
#Raytracer_exe = os.path.join(Raytracer_Path,'The_Raytracer.jar')
Raytracer_exe = os.path.join(Raytracer_Path,'The_Raytracer.jar')
path_input_data = ("D:\Raytracer_data\Raytracer_local_input")
path_output_data = ("D:\Raytracer_data\Raytracer_local_output")
path_progress_data = ("D:\Raytracer_data\Progress")
File_name = '200x200mm_hole'

#Set folder for third party libary
Path_poisson_recon = ("D:\Raytracer\Programme\PoissonRecon")

#File mit Simulations- und Remeshdaten (Performance)
time_log_filename = os.path.join(path_progress_data,File_name + '_logfile')
time_log_fil = open(os.path.join(time_log_filename + '.txt'),'w')
time_log_fil.write("Pulse Number; Tracing time; Remesh time; max z position; max ablation dist; min z position \n")
time_log_fil.close()

#File zum Speicher der Schnittprofile nach jedem Puls
profile_log_filename = os.path.join(path_progress_data,File_name + '_drilling_progress')
profile_log_file = open(os.path.join(profile_log_filename + '.txt'),'w')
profile_log_file.write("x and y profiles after single pulses \n")
profile_log_file.close()

n_pulses = 300
for i in range(1, n_pulses+1):
    if i%10 == 0:
        print (i)
    
    # STL zu TRI konvertieren
    print('STL konvertieren')
    if i==1:
        stl_file = os.path.join(path_input_data,File_name +".stl")
        tri_file = os.path.join(path_input_data,File_name +".tri")
    else:
        if 'remesh_filename' in locals():
            print(i)
            stl_file = os.path.join(path_input_data,remesh_filename +".stl")
            tri_file = os.path.join(path_input_data,remesh_filename +".tri")
        else:
            remesh_filename = os.path.join(File_name + '_pulse_' + str(i))
            remesh_filename_input = os.path.join(File_name + '_pulse_' + str(i-1))
            stl_file = os.path.join(path_input_data,remesh_filename_input +".stl")
            tri_file = os.path.join(path_input_data,remesh_filename_input +".tri")

    Algorithmen.stl_to_tri(stl_file, tri_file)
    #Algorithmen.STL2tri(stl_file, tri_file)
    
    
    #Raytracing calculation
    start_time_tracen = time.time()
    print('tracen')
    Parameterfile = os.path.join(Raytracer_Path, 'Parameter_Raytracer.txt')
    if i==1:
        OUTfile = os.path.join(path_output_data,File_name + '_pulse_' + str(i))
    else:
        OUTfile = os.path.join(path_output_data,File_name + '_pulse_' + str(i))
    raytracer_arguments = [OUTfile, Parameterfile, tri_file] 
    cmd = "java -jar "+ Raytracer_exe +' '+ " ".join(raytracer_arguments)
    sp.call(cmd)
    print('finished tracing')
    delta_time_tracen = round(time.time() - start_time_tracen, 6)
    
    
    #Excel File auslesen 
    f = open(os.path.join(path_output_data,File_name + '_pulse_' + str(i) + '_output.csv'), 'r') # open the file for reading
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
    
    
    #Excel File auslesen nummer 2 (buggy weil beim ersten auslesen manchmal murks rauskommen)
    f = open(os.path.join(path_output_data,File_name + '_pulse_' + str(i) + '_output.csv'), 'r') # open the file for reading
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
    SPs = basic_data[:,1:4] /1e3 # Schwerpunkte in mm
    n_points = len(SPs)
    Is = basic_data[:,10] #Intensitäten
    Nvecs = basic_data[:,4:7]
    Kvecs = basic_data[:,7:10] 
    
    
    #Parameter für Abtrag
    #Hier muss noch die Schwellfluenz in Abhöängigkeit der Pulsnummer berücksichtigt werden (Workaround)
    Phi_th = 0.0283 # Schwellfluenz des Materials in J/cm^2
    delta_opt = 1e-5 # optische eindringtiefe in mm
    tau_pulse = 8e-12 # Pulsdauer in s
    
    
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
    pointcloud_filename = os.path.join(File_name + '_pulse_' + str(i) + '_pointcloud')
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
    start_time_remesh = time.time()
    print('Mashlab remashing')
    remesh_filename = os.path.join(File_name + '_pulse_' + str(i)) 
    cmd_meshlab = "r"
    
    
    #Remesh ohne Meshlab (zum teil)
    #data_in = ("D:\Programmierung\Poisson_Recon\testing\hole.ply")
    data_in = os.path.join(pointcloud_filename + '.ply')
    #data_out = ("D:\Programmierung\Poisson_Recon\testing\holeneu.ply")
    data_out = os.path.join(pointcloud_filename + '_mesh.ply')
    depth = 10
    scale = 1.003
    samplePerNode = 3
    pointWeight = 2
    iters = 9
    cmd = Path_poisson_recon + " --in " + data_in + ' --out '+ data_out + ' --depth '+ str(depth) + ' --scale '+ str(scale) + ' --samplePerNode '+ str(samplePerNode) + ' --pointWeight '+ str(pointWeight) + ' --iters '+ str(iters)
    sp.call(cmd)
        
    #Pointcloud_test_out = mlx.FilterScript(file_in = 'Pointcloud_test.ply', file_out = 'Pointcloud_remesh.stl', ml_version='2016.12')
    Pointcloud_test_out = mlx.FilterScript(file_in = data_out, file_out = os.path.join(remesh_filename + '.stl'), ml_version='2016.12')
#    mlx.remesh.surface_poisson_screened(Pointcloud_test_out,depth=10,
#                                 full_depth=5, cg_depth=0, scale=1.0092,
#                                 samples_per_node=3, point_weight=2.0,
#                                 iterations=8, confidence=False, pre_clean=False)
#    mlx.layers.delete(Pointcloud_test_out)
    mlx.remesh.simplify(Pointcloud_test_out, texture=False, faces=10000, target_perc=0.0,
                  quality_thr=0.95, preserve_boundary=True, boundary_weight=10.0,
                  optimal_placement=True, preserve_normal=True,
                  planar_quadric=False, selected=False, extra_tex_coord_weight=1.0,              
                  preserve_topology=True, quality_weight=False, autoclean=True)
    Pointcloud_test_out.run_script()
    delta_time_remesh = round(time.time() - start_time_tracen, 6)
    os.remove(data_out)
    
    #convert to ascii encoding
    print('converting...')
    f = stl.mesh.Mesh.from_file(os.path.join(Raytracer_Path,remesh_filename + '.stl'))
    f.save(os.path.join(Raytracer_Path,remesh_filename + '.stl'), mode=stl.Mode.ASCII)
    
    #Schnittprofile in x/y-Richtung mit bestimmter breite bestimmen um fortschritt zu charackterisieren
    [x_profile_cut, y_profile_cut] = Algorithmen.get_profile_cut(new_SPs,0,0,0.005)
    print('Write progress to log file')
    with open(os.path.join(profile_log_filename + '.txt'),'a') as f:
        #Array in Liste umwandeln um Zeilenweise abzuspeichern
        np.savetxt(f, [x_profile_cut[:,0]], fmt='%2.10f', delimiter=';')
        np.savetxt(f, [x_profile_cut[:,1]], fmt='%2.10f', delimiter=';')
        np.savetxt(f, [y_profile_cut[:,0]], fmt='%2.10f', delimiter=';')
        np.savetxt(f, [y_profile_cut[:,1]], fmt='%2.10f', delimiter=';')
    
    
    #Logfile füllen
    print('Write to log-file')
    with open(os.path.join(time_log_filename + '.txt'),'a') as log:
        log.write(str(i) + ";" + str(delta_time_tracen) + ";" + str(delta_time_remesh) + ";" +str(np.amax(SPs,axis=0)[2]) + ";" +str(np.amax(ablation_dist)) + ";" +str(np.amin(SPs,axis=0)[2]) + "\n")
    
    #Pointcloud/remeshed file verschieben
    print('move new files to folders')
    shutil.move(os.path.join(Raytracer_Path,pointcloud_filename + '.ply'), pointcloud_folder)
    #shutil.move(os.path.join(Raytracer_Path,remesh_filename + '.stl'), new_mesh_folder)
    shutil.move(os.path.join(Raytracer_Path,remesh_filename + '.stl'), path_input_data)
    
    #Überflüssige Files löschen
    if (i>3 and i%50 != 0):
        if os.path.isfile(os.path.join(Raytracer_Path,path_input_data,File_name + '_pulse_' + str(i-2) + '.stl')):
            print('delete the fuck')
            prev_remesh_filename = os.path.join(File_name + '_pulse_' + str(i-2) + '.stl')
            prev_tri_filename = os.path.join(File_name + '_pulse_' + str(i-2) + '.tri')
            os.remove(os.path.join(Raytracer_Path,path_input_data,prev_remesh_filename )) #Remeshfile entfernen
            os.remove(os.path.join(Raytracer_Path,path_input_data,prev_tri_filename )) #Remeshfile entfernen (tri)
        excel_file = os.path.join(path_output_data,File_name + '_pulse_' + str(i) + '_output.csv')
        intensity_file = os.path.join(OUTfile + '_intensity_out.ply')
        settings_file = os.path.join(path_output_data,File_name + '_pulse_' + str(i) + '_settings.log')
        os.remove(os.path.join(pointcloud_folder,pointcloud_filename + '.ply')) #Pointcloundfile entfernen  
        os.remove(os.path.join(Raytracer_Path,path_output_data,settings_file)) #Settings file entfernen
        os.remove(os.path.join(Raytracer_Path,path_output_data,intensity_file)) #Intensitätsverteilungsfile entfernen
        os.remove(os.path.join(Raytracer_Path,path_output_data,excel_file))#excel file entfernen


