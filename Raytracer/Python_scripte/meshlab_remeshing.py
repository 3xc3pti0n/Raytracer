# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:42:40 2019

@author: Julian
"""
import subprocess as sp

cmd0 = "C:"
cmd1 = "cd \Program Files\VCG\MeshLab" 
cmd2 = 'meshlabserver -i D:\Raytracer\Pointclounds\Pointcloud_test.ply -o D:\Raytracer\remesh\Pointcloud_test_remesh.stl  -s D:\Raytracer\Meshlab_scripte\poisson_surface_reconstruction.mlx'

run1 = sp.Popen(['runas', '/noprofile', '/user:Administrator',cmd0],stdin=sp.PIPE, shell=True)
run1.stdin.write(b'serv-MK3')
run1.communicate()
run2 = sp.Popen(['runas', '/noprofile', '/user:Administrator',cmd1],stdin=sp.PIPE, shell=True)
run2.stdin.write(b'serv-MK3')
run2.communicate()
run3 = sp.Popen(['runas', '/noprofile', '/user:Administrator',cmd2],stdin=sp.PIPE, shell=True)
run3.stdin.write(b'serv-MK3')
run3.communicate()