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

class elements():
    def __init__(self, filename, Z):
        self.v = zeros((3,3,3,3), dtype = float)
        self.h_HF = zeros((6,6), dtype = float)
        el = open(filename, 'r')
        for l in el:
            e = l.split(' ')
            e[4] = e[4][0:-2]
            self.v[int(e[0]),int(e[1]),int(e[2]),int(e[3])] = float(e[4])
        self.init_spin()
    def init_spin(self):
        #Translates the 
        self.V = zeros((6,6,6,6), dtype = float)
        for p in range(6):
            for q in range(6):
                for r in range(6):
                    for s in range(6):
                        D = self.v[p/2,q/2,r/2,s/2] #Direct term
                        V = self.v[p/2,q/2,s/2,r/2] #Exchange term
                        self.V[p,q,r,s] = self.state(p,q,r,s,D,V)
    def state(self,p,q,r,s,V,D):
        #Evaluates the symmtrization conditions
        s1,s2,s3,s4 = p%2, q%2, e%2, s%2
        if s1==s2:
            if s3==s4:
                if s1==s2:
                    S = D - V #both terms
            else:
                S = 0 #integral cancels due to spin-orthogonality
        if s1!=s2:
            if s3!=s4:
                if s1==s3:
                    S = D #Direct term only
                else:
                    S = V #Exchange term only
            else:
                S = 0 #integral cancels due to spin-orthogonality
        return S
    def e(self,p,q,r,s):
        #Returns properly symmetrized matrix elements
        return self.V[p,q,r,s]
    def HF(self, C):
        #Set up the HF matrix, using the coefficients in C
        #self.C = zeros((6,6), dtype = float)
        #self.C[range(6), range(6)] = 1.0
        self.h_HF *= 0
        for alpha in range(6):
            for gamma in range(6):
                interaction = 0
                for p in range(6):
                    for beta in range(6):
                        for delta in range(6):
                            interaction += C[p,beta]*C[p,delta]*self.V[alpha,beta,gamma,delta]
                self.h_HF[alpha,gamma] = self.h0(alpha,gamma)+interaction
    def h0(self,alpha,gamma):
        #One-body interaction
        h = 0
        if alpha==gamma:
            #Return one body energy in potential
            #THIS IS NOT YET PROPERLY IMPLEMENTED
            h = 1
        return h
    def HF_solve(self):
        self.C = zeros((6,6), dtype = float)
        self.C[range(6), range(6)] = 1.0
        e_v = zeros(6, dtype = float)
        e_v_prev = ones(6, dtype = float)
        tolerance = 10**-14
        iters = 0
        while abs(min(e_v)-min(e_v_prev))>tolerance:
            iters += 1
            e_v_prev = e_v
            self.HF(self.C)
            e_v, self.C = eig(self.h_HF)
        print iters, min(e_v)


b = elements('m_elements_c.dat', 1)
b.HF_solve()



                    