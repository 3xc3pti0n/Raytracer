# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 15:22:24 2019

@author: Julian
"""

import numpy as np

import STLpy
import os
import csv
import subprocess

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
    
def stl_to_tri(inputfile, outputfile):
    """
    Die Funktion konvertiert eine *.stl-Datei (ASCII) zu einer *.tri-Datei
    """
    points = STLpy.read_ASCII_STL(inputfile)

    counter=0
    with open(outputfile, "w") as f:
        for i in range(0, len(points), 3):
            f.write("TRIANGLE " + str(counter) + "\n")
            f.write(str(points[i, 0]) + ",\t" + str(points[i, 1]) + ",\t"
                    + str(points[i, 2]) + "\n")
            f.write(str(points[i+1, 0]) + ",\t" + str(points[i+1, 1]) + ",\t"
                    + str(points[i+1, 2]) + "\n")
            f.write(str(points[i+2, 0]) + ",\t" + str(points[i+2, 1]) + ",\t"
                    + str(points[i+2, 2]) + "\n")
            counter = counter + 1
   
         
def read_ascii_stl_as_tris(filename):
    """
    Liest eine ascii-Datei ein und gibt alle Dreiecke als (ntris,3,3) Array 
    zurück
    """
    f=open(filename)
    lines = f.readlines()
    lines = [l for l in lines if 'vertex' in l]
    data = np.asanyarray([[float(x) for x in [e for e in l.rstrip('\n').split(' ') if e != ''][1:]] for l in lines])
    tris = np.reshape(data, (int(len(data)/3), 3, 3))
    f.close()
    return tris


def read_ascii_stl_facet_normals(filename):
    """
     Liest eine ascii-Datei ein und gibt alle Normalen der Dreiecke als 
     Array zurück
    """
    f = open(filename)
    lines = f.readlines()
    lines = [l for l in lines if 'facet normal' in l]
    normals = np.asanyarray([[float (x) for x in [e for e in l.rstrip('\n').split(' ') if e != ''][2:]] for l in lines])
    f.close()
    return normals


def get_absorbed_intensity(path_data, length_of_tris):
    """
    liefert die absorbierten Intensität der Patches, die zuvor im Raytracer 
    berechnet wurden
    """
    # convert from W/((µm)²) to W/(m²)
    absorbed_intensity = read_raytracer_output(path_data, length_of_tris)*1.0e6
    return absorbed_intensity


def read_raytracer_output(csv_data, length_of_tris):
    """
    Parser für die vom Raytracer erstellte Datei kapillare_output.csv
    """
    absorbed_intensity = np.empty(length_of_tris)
    
    with open(csv_data, newline="\n") as csv_file:
        raytracer_output = csv.reader(csv_file, delimiter=";")
        for row in raytracer_output:
            
            # skip empty lines and the header
            if (len(row) > 0) and str.isdigit(row[0]):
                
                absorbed_intensity[int(row[0])] = float(row[4])

    return absorbed_intensity


def get_absorbed_power(absorbed_intensity, tris):
    """
    Berechnung der absorbierten Leistung pro Dreieck
    """
    absorbed_power = np.empty(len(tris))
    
    for i in range(0, len(tris)):
        #    Vektoren des Patches
        vector1 = tris[i, 1] - tris[i, 0]
        vector2 = tris[i, 2] - tris[i, 0]
    
        # Fläche des Patches
        A_patch = 0.5*np.linalg.norm(np.cross(vector1, vector2))/1.0e6
        absorbed_power[i] = absorbed_intensity[i] * A_patch
        
    return absorbed_power


def run_shell_command(cmd):
    """
    executes a given command in the shell
    """
    print("Starting " + cmd)
    subprocess.run([cmd], shell=True)    
    print("Finished " + cmd) 


def start_raytracer(run_path):
    """
    calls the raytracer
    """
    raytracer_arguments = ["", "", ""]
    # name of the output file
    raytracer_arguments[0] = str(run_path/"OUT/kapillare")
    # name of the parameter file
    raytracer_arguments[1] = str(run_path/ "Parameter_Raytracer.txt")
    # name of the input file
    raytracer_arguments[2] = str(run_path/"Surfaces/kapillare.tri")
    print("calling the raytracer")
    run_shell_command("java -jar raytracer.jar " + " ".join(raytracer_arguments))
    
    
def tris_center_of_mass (tris):
    """
    Berechnung des Masseschwerpunktes jedes Dreiecks
    """
    tris_com = np.zeros([len(tris),3])
    for i in range(0,len(tris)):
        tris_com[i,0] = 1/3 * (tris[i,0,0]+tris[i,0,1]+tris[i,0,2]) 
        tris_com[i,1] = 1/3 * (tris[i,1,0]+tris[i,1,1]+tris[i,1,2]) 
        tris_com[i,2] = 1/3 * (tris[i,2,0]+tris[i,2,1]+tris[i,2,2]) 
        
    return tris_com

def tris_area (tris):
    """
    Berechnung der Fläche jedes Dreiecks (Für eventuelle Teilung zu großer
    Dreiecke nach Ball-pivot Algorithmus)
    """
    tris_a = np.empty(len(tris))
    for i in range(0,len(tris)):
        vector1 = tris[i, 1] - tris[i, 0]
        vector2 = tris[i, 2] - tris[i, 0]
        tris_a[i] = 0.5*np.linalg.norm(np.cross(vector1, vector2))/1.0e6 
        
    return tris_a


def calc_ablat_dist (tris_int, Phi_th, tau_pulse, delta_opt):
    fluence = tris_int * 1e8 * tau_pulse #Umrechnung von W/Mikrometer^2 in W/cm^2
    ablat_z = np.zeros([len(fluence),1])
   # tris_phi = intensities*t_p
    
    ablat_z = np.where((fluence != 0) & (fluence >= Phi_th),delta_opt*np.log(fluence/Phi_th),0)
    
    return ablat_z

#Methode die ein Schnittprofil entlang einer gegebenen Achse in x- und y-Richtung zurückgibt
def get_profile_cut(points, x_line, y_line, width):
    
    x_SP = np.where((points[:,0] < (x_line + width)) & (points[:,0] > (x_line - width)),0,1)
    y_SP = np.where((points[:,1] < (y_line + width)) & (points[:,1] > (y_line - width)),0,1)
    x_range = np.where(x_SP)
    y_range = np.where(y_SP)
    x_profile = np.delete(points, x_range, axis=0)
    x_profile = np.delete(x_profile, 0, axis=1)
    y_profile = np.delete(points, y_range, axis=0)
    y_profile = np.delete(y_profile, 1, axis=1)
    
    x_uniques = np.unique(x_profile[:,0])
    y_uniques = np.unique(y_profile[:,0])
    
    mean_depth_x = np.zeros([len(x_uniques),2])
    mean_depth_y = np.zeros([len(y_uniques),2])
    
    mean_depth_x[:,0] = x_uniques
    mean_depth_y[:,0] = y_uniques
    
    number_x=0
    for i in x_uniques:
        faktor = sum(np.where(x_profile[:,0] == i,1,0))
        summe = sum(np.where(x_profile[:,0] == i,x_profile[:,1],0))
        mean_depth_x[number_x,1] = summe/faktor
        number_x += 1
        
    number_y=0
    for i in y_uniques:
        faktor = sum(np.where(y_profile[:,0] == i,1,0))
        summe = sum(np.where(y_profile[:,0] == i,y_profile[:,1],0))
        mean_depth_y[number_y,1] = summe/faktor
        number_y += 1

    return mean_depth_x, mean_depth_y



