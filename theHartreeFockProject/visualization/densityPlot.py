# -*- coding: utf-8 -*-
"""
Created on Thu May 15 10:04:45 2014

@author: audunsh
"""

from numpy import *
import Image

def RAWsave(ray,filename):
    savaetxt(filename, ray, delimiter="   ")

def RAWload(filename):
    return loadtxt(filename,dtype=float,delimiter="   ")

def Imray(fname):
    #open image to array
    im = Image.open(fname).convert("L")
    ray = asarray(im)
    return ray

def Rayim(IM, fname):
    #save array as image 
    im = Image.fromarray(IM).convert("L")
    im.save(fname)

class cubed():
    #loads data from filename and returns coordinatewise data
    def __init__(self, filename, N):
        self.N = N
        self.dataset = fromfile("testmap2.dat", dtype = float32)
        self.cube = zeros((N,N,N), dtype = float)
        for x in range(self.N):
            for y in range(self.N):
                for z in range(self.N):
                    self.cube[x,y,z] = self.at(x,y,z)
    def at(self, x,y,z):
        return self.dataset[x*self.N*self.N + y*self.N + z]
    def size(self):
        return len(self.dataset)



#a = cubed("testmap3.dat",39)
#Rayim(a.cube[0,:,:], "tesmap3.png")

N = 1000 #dimensions of plot
a = loadtxt("slice6.dataset")
img = zeros((N,N), dtype = float)
#img = a.resize(1000,1000)

for i in range(N):
    for e in range(N):
        img[i,:] =a[i]

#normalizing
img *= 256/img.max()

Rayim(img, "slice6.png")

#print sqrt(size(a))

