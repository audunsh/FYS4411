# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:13:56 2014

#Script to read out matrix elements 
"""

from numpy import *


el = open('m_elems.dat')

table = zeros((4,4,4,4), dtype = float)

for i in el:
    elems = i.split('$\langle')
    for e in elems:
        if len(e) != 0:
            m_e = e[22:-3]
            #m_e = m_e.replace("(", "float(")
            m_e = m_e.replace("Z", "")
            m_e = m_e.replace('\sqrt{',"*sqrt(")
            m_e = m_e.replace('}', ")")
            m_e = m_e.replace("$", "")
            m_e = m_e.replace(" ", "")
            m_e = m_e.replace("/", "*1.0/")
            print eval(m_e), e[1], e[2], e[6], e[7]           
            table[int(e[1])-1, int(e[2])-1, int(e[6])-1, int(e[7])-1] = eval(m_e)

print table[2,2,3,1]

f = open("m_elements_c.dat", "w")

for i in range(3):
    for e in range(3):
        for u in range(3):
            for y in range(3):
                s = str(i)+" "+str(e)+" "+str(u)+" " +str(y) + " " + str(table[i,e,u,y]) + '\n'
                f.write(s)

f.close()
            