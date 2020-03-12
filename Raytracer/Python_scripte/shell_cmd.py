# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 14:32:49 2019

@author: Julian
"""
import subprocess

def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print (proc_stdout)

subprocess_cmd('echo a; echo b')