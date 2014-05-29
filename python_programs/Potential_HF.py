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

# Default values:
n = 100             # number of Hartree-Fock calculations (# of core distances R)
Z = 8               # number of protons in core
N = 8               # number of electrons orbiting a core
Ns = 6              # number of single particle states

val1 = 0
val2 = 0
x_list = [2,3,8,10,12,14,16,18]

Atom = 'some_atom'
#----------------------------------------------------------------------------

if (Z == 1):
    atom = 'H'
    Atom = 'Hydrogene'
elif (Z == 2):
    atom = 'He'
    Atom = 'Helium'
elif (Z == 4):
    atom = 'Be'
    Atom = 'Beryllium'
elif (Z == 8):
    atom = 'O'
    Atom = 'Oxygene'
elif (Z == 10):
    atom = 'Ne'
    Atom = 'Neon'
else:
    atom = 'not-specified'
    Atom = 'not-specified'

print 'atom='
print Atom

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
    
    print i
    # args = [#protons, #electrons, #distance from atom]
    args = ["%d" % Z , "%d" % N , "%f" % R[i,0], atom]
    #print("args = ", args)
    En = sub.check_output(["./HartreeFock", args[0], args[1], args[2], args[3]], cwd=path_to_HartreeFock)
    #print En
    E[i] = float(En)    
           
#----------------------------------------------------------------------------

dx = R[1]-R[0]
F = np.zeros((n,1))
for i in range(n-1):
    F[i] = -(E[i+1] - E[i])/dx

lenF = len(F)

##############################################################################
# The plotting:    
import matplotlib.pyplot as plt    

#---------------------------------------------------------------------------
ShowPlot = raw_input('Show plots? (yes/no): ')
showplot = False
saveplot = False
if (ShowPlot == 'yes' or ShowPlot == 'Yes' or ShowPlot == 'Y' or ShowPlot == 'y' or ShowPlot == 'true' or ShowPlot == 'True'):
    showplot = True
    print '----------------------------------------------------------------'
    print 'Show plot = True'
    print '----------------------------------------------------------------'
else:
    'Plots are not showed'
    
SavePlot = raw_input('Save plots? (yes/no): ')
if (SavePlot == 'yes' or SavePlot == 'Yes' or SavePlot == 'y' or SavePlot == 'Y' or SavePlot == 'true' or SavePlot == 'True'):
    saveplot = True
    print 'Saving plots...'
    print '----------------------------------------------------------------'
else:
    print 'Plots are not saved'
    #---------------------------------------------------------------------------

h = plt.figure()
plt.plot(R,E,'b-*')
plt.hold(True)
plt.plot(R,np.zeros(np.size(R)),'r-')
plt.title('Potential distribution for the %s atom. Z=%g, N=%g' % (atom,Z,N))
plt.legend(('E(r)','0'))
plt.xlabel('r [a.u.]')
plt.ylabel('Potential [a.u]')
if (saveplot == True):
    name1 = 'potential_%s_.png' % atom
    plt.savefig(name1, format='png')
    

hh = plt.figure()
plt.plot(R[0:lenF],F[0:lenF],'r-d')
plt.title('Force on %s atom from another. Z=%g, N=%g' % (atom,Z,N))
plt.legend('F(r)')
plt.xlabel('r [a.u.]')
plt.ylabel('Force [a.u.]')
if (saveplot == True):
    name2 = 'force_%s_.png' % atom
    plt.savefig(name2, format='png')
    
plt.show(showplot)
