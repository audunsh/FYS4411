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
Z = 4               # number of protons in core
N = 4               # number of electrons orbiting a core
Ns = 6              # number of single particle states
#----------------------------------------------------------------------------

# Varyation of the core distances:
R = np.zeros((n,1))
E = np.zeros((n,1))
for i in range(len(R)):
    R[i] = (1+i)*0.1   # just changing the x-position
    
    #E[i] = os.system("path/to/exec.exe %var1 %var2 %var3" % (var1, var2, var3))
    #E[i] = os.system("LD_LIBRARY_PATH="~/goran/CompPhys/FYS4411 - CompPhys2/build-theHartreeFockProject-Desktop-Release/src/libs/" ~/goran/CompPhys/FYS4411\ -\ CompPhys2/build-theHartreeFockProject-Desktop-Release/src/HartreeFock/HartreeFock %g %g %f" % (Z,N,R[i,0]))
#    returnval = os.popen('python ~/goran/teaching/FYS-MEK1110-mechanics/TrappingAtoms.py',"r")
#    while 1:
#        line = returnval.readline()
#        if not line: break
#        print line
    returnval = os.popen('LD_LIBRARY_PATH="/home/goranbs/goran/CompPhys/FYS4411 - CompPhys2/build-theHartreeFockProject-Desktop-Release/src/libs/" ~/goran/CompPhys/FYS4411\ -\ CompPhys2/build-theHartreeFockProject-Desktop-Release/src/HartreeFock/HartreeFock %g %g %f' % (Z,N,R[i,0]),"r")
    while 1:
        line = returnval.readline()
        if not line: break
        print line
    
#----------------------------------------------------------------------------
# The plotting:
    
import matplotlib.pyplot as plt    

h = plt.figure()
plt.plot(R,E,'b-*')
plt.title('Potential distribution for Z=%g N=%g Ns=%g' % (Z,N,Ns))
plt.legend('E(R)')
plt.xlabel('R [a.u.]')
plt.ylabel('Potential [a.u]')
plt.show(True)
