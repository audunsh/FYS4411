# -*- coding: utf-8 -*-
"""
Created on Fri May  9 08:54:46 2014

@author: goranbs

FYS4411 - Quantum Mechanical Systems - Spring 2014

This program should call the c++ Hartree-Fock code and run for different distances R for two atoms.
The output from the Hartree-Fock code is the energy of the ground state, given the distances R between
a system of two particles. 

You also give the number of protons in the core; Z
The total number of electrons;                   N
And the number of single particle states;        Ns

The default is to use a STO-3G basis.
"""
import os
import subprocess as sub
import numpy as np

n = 100             # number of Hartree-Fock calculations (# of core distances R)
Z = 2               # number of protons in core
N = 2               # number of electrons orbiting a core
Ns = 6              # number of single particle states
#----------------------------------------------------------------------------

current_path =  os.path.realpath("Potential_HF.py")
path_to_Release = os.path.abspath(os.path.join(current_path, "..", "..", "build-theHartreeFockProject-Desktop-Release"))

path_to_HartreeFock = os.path.abspath(os.path.join(path_to_Release, "src", "HartreeFock"))
path_to_libs = os.path.abspath(os.path.join(path_to_Release, "src", "libs"))

os.environ["LD_LIBRARY_PATH"] = path_to_libs


# Varyation of the core distances:
R = np.zeros((n,1))
E = np.zeros((n,1))
for i in range(len(R)):
    R[i] = (1+i)*0.1   # just changing the x-position
    
#    returnval = os.popen('LD_LIBRARY_PATH="/home/goranbs/goran/CompPhys/FYS4411 - CompPhys2/build-theHartreeFockProject-Desktop-Release/src/libs/" ~/goran/CompPhys/FYS4411\ -\ CompPhys2/build-theHartreeFockProject-Desktop-Release/src/HartreeFock/HartreeFock %g %g %f' % (Z,N,R[i,0]),"r")
#    returnval = os.popen('LD_LIBRARY_PATH="/home/goranbs/goran/CompPhys/FYS4411 - CompPhys2/build-theHartreeFockProject-Desktop-Debug/src/libs/" /home/goranbs/goran/CompPhys/FYS4411\ -\ CompPhys2/build-theHartreeFockProject-Desktop-Debug/src/HartreeFock/HartreeFock %g %g %f' % (Z,N,R[i,0]),"r")    
    
    args = ["%d" % Z , "%d" % N , "%f" % R[i,0]]
    #print("args = ", args)
    
    #En = sub.check_output(["/home/goranbs/goran/CompPhys/FYS4411 - CompPhys2/build-theHartreeFockProject-Desktop-Release/src/HartreeFock/HartreeFock", args[0], args[1], args[2]])
    En = sub.check_output(["./HartreeFock", args[0], args[1], args[2]], cwd=path_to_HartreeFock)
    E[i] = float(En)    
    #print En

#    while 1:
#        line = returnval.readline()
#        if not line: break
#        E[i] = float(line)
#        print line
        
#    p = sub.Popen('LD_LIBRARY_PATH="/home/goranbs/goran/CompPhys/FYS4411 - CompPhys2/build-theHartreeFockProject-Desktop-Release/src/libs/" ~/goran/CompPhys/FYS4411\ -\ CompPhys2/build-theHartreeFockProject-Desktop-Release/src/HartreeFock/HartreeFock %g %g %f' % (Z,N,R[i,0]),stdout=sub.PIPE,stderr=sub.PIPE)
#    output, errors = p.communicate()
#    print output
        
           
#----------------------------------------------------------------------------

dx = R[1]-R[0]
F = np.zeros((n,1))
for i in range(n-1):
    F[i] = (E[i+1] - E[i])/dx

lenF = len(F)

# The plotting:    
import matplotlib.pyplot as plt    

h = plt.figure()
plt.plot(R,E,'b-*')
plt.title('Potential distribution for the Hydrogen atom')
plt.legend('E(r)')
plt.xlabel('r [a.u.]')
plt.ylabel('Potential [a.u]')

hh = plt.figure()
plt.plot(R[0:lenF],F[0:lenF],'r-d')
plt.title('Force on Hydrogen Atom from another')
plt.legend('F(r))')
plt.xlabel('r [a.u.]')
plt.ylabel('Force [units?]')

plt.show(True)
