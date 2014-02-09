# -*- coding: utf-8 -*-
"""
First attempt at Hartree Fock using python
1. Import all matrix elements
2. Set up HF-matrix, initial guess for C = I
3. Solve eigenvalueproblem using eig
4. Repeat procedure using the eigenvectors from 3
5. ...until convergence in lowest eigenvalue.

"""

from numpy import *


#Reading matrix elements into array m_elems
f = open('m_elems.dat', 'r')
m_elems = zeros((3,3,3,3), dtype = float)
for l in f:
    e = l.split(' ')
    m_elems[int(e[0]),int(e[1]),int(e[2]),int(e[3])] = float(e[4])

#Initial C
C = zeros((6,6), dtype = float)
C[range(6), range(6)] = 1.0
#Setting up HF-matrix

def HF(C,m_elems,h0, Np):
    N = len(C) #number of states
    h_HF = zeros((N,N), dtype = float)
    for alpha in range(N):
        for gamma in range(N):
            interaction = 0
            for p in range(Np):
                for beta in range(N):
                    for delta in range(N):
                        interaction += C[p,beta]*C[p,delta]*m_elems[alpha,beta,gamma,delta]
            h_HF[alpha,gamma] = h0[alpha,gamma]+interaction
    return h_HF


                    