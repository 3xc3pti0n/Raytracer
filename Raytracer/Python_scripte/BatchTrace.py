# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 14:17:12 2019

@author: ifsw-ffetzer
"""

import os
import subprocess
import numpy as np
import stl
import shutil
import random

Raytracer_Path = '.'
Raytracer_Path_input = '.\\Raytracer_local_input'
Raytracer_Path_output = '.\\Raytracer_local_output'
Raytracer_exe = os.path.join(Raytracer_Path,'The_Raytracer.jar')

def STL2tri(STLfile, TRIfile):
    tristring = ''
    counter_triangle = -1
    f = open(STLfile)
    lines = f.readlines()
    f.close()
    for l in lines:
        if 'facet normal' in l:
            #print(counter_triangle)
            counter_triangle = counter_triangle + 1
            tristring = tristring + 'TRIANGLE ' + str(int(counter_triangle)) + '\n'
        if 'vertex' in l:
            tristring = tristring + l.split()[1] + ',    ' + l.split()[2] + ',    ' + l.split()[3] + '\n'
    f = open(TRIfile,'w')
    f.write(tristring)
    f.close()

akt_kap_name = 'act_cap.stl'
skript_mesh = 'points_to_mesh.mlx'
skript_simply = 'mesh_simply.mlx'
times= np.linspace(0,1,50)
n_points=-1
for t in times:
    print("Now at "+ str(t))
#    # STL to TRI
    STLfile = os.path.join(Raytracer_Path_input,akt_kap_name)
    TRIfile = os.path.join(Raytracer_Path_input,akt_kap_name[:-4] + '.tri')
    STL2tri(STLfile, TRIfile)
    TRIfile = os.path.join(Raytracer_Path_input, akt_kap_name[:-4] + '.tri')
#    # tracen
    print('tracen')
    Parameterfile = os.path.join(Raytracer_Path, 'Parameter_Raytracer.txt')
    OUTfile = os.path.join(Raytracer_Path_output, 'cap_'+str(t))
    raytracer_arguments = [OUTfile, Parameterfile, TRIfile] # out, Parameter, in (tri)
    cmd = "java -jar "+ Raytracer_exe +' '+ " ".join(raytracer_arguments)
    subprocess.run(cmd.split())
    print('finished tracing')
   # verschiebung berechnen    
    # get SP
    f = open(os.path.join(Raytracer_Path_output,'cap_'+str(t)+'.csv'))
    data = np.asanyarray([[float(x) for x in l.rstrip('\n').split('\t')] for l in f.readlines()[1:]])
    SPs = data[:,1:4]/1000. # get SPs und auf mm umrechnen
    if n_points == -1: # only in first iteration nr points is set
        n_points = len(SPs)
    Is = data[:,4] # get Intensities
    # get normal
    i = stl.mesh.Mesh.from_file(STLfile)
    norms = np.linalg.norm(i.normals,2,1)
    normals = i.normals/norms.reshape((len(norms),1))
    #normals = np.multiply(i.normals, np.linalg.norm(i.normals,2,1)) # auf 1 normiert
    verschiebungen = Is/np.max(Is)/10. # einfach um maximal 0.2mm verschieben.
    verschiebungen = verschiebungen.reshape((len(verschiebungen),1))
    # verschiebe punkte
    SPsnew = SPs - normals*verschiebungen
     # schreibe punkte in xyz file
    f = open('new_Sps.xyz','w')
    print('Nr. points before: ' + str(len(SPsnew)))
    sample = random.sample(range(len(SPsnew)), n_points)
    SPswrite = SPsnew[sample]
    print('Nr. points before: ' + str(len(SPswrite)))
    for s in SPswrite:
        f.write(str(s[0]) + ' '+str(s[1])+' '+ str(s[2])+'\n')
    f.close()
#    # meshlab BPA und Taubin
    fname_in = os.path.join('.', "new_Sps.xyz")
    fname_out = os.path.join(Raytracer_Path_input, 'cap_' + str(t)+'.stl')
    print('meshlabbing..stl generation.')
    subprocess.run(["meshlabserver.exe", "-i",fname_in ,"-o",fname_out,  '-s', skript_mesh])
#    
#    print('meshlabbing..stl simplification.')
#    subprocess.run(["meshlabserver.exe", "-i",fname_out ,"-o",fname_out,  '-s', skript_simply])
#    
    print('converting...')
    i = stl.mesh.Mesh.from_file(fname_out)
    i.save(fname_out, mode=stl.Mode.ASCII)# save as ASCII
    shutil.copy(fname_out, STLfile) # copy to akt_STL file
#    # slt in ASCII
    print("finished " + str(t))