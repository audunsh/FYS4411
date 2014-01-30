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
            

            #print m_e, eval(m_e)                
            table[e[1], e[2], e[6], e[7]] = eval(m_e)

print table[2,2,3,1]
            