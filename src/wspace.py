#!/usr/bin/env python

import sys
import numpy as np
import math

def set_LR_matrix(L, R, dens, velx, pres):
    gamma=1.4
    cs = np.sqrt((gamma-1.0)*pres/dens)
    enth = 0.5*velx*velx + (gamma*pres)/((gamma-1.0)*dens)

    R[0][0] = 1.0
    R[0][1] = 1.0
    R[0][2] = 1.0
    R[1][0] = velx-cs
    R[1][1] = velx
    R[1][2] = velx+cs
    R[2][0] = enth-velx*cs
    R[2][1] = 0.5*velx*velx
    R[2][2] = enth+velx*cs

    b1 = 0.5*(gamma-1.0)*(velx*velx)/(cs*cs)
    b2 = (gamma-1.0)/(cs*cs)

    L[0][0] = 0.5*(b1 + velx/cs)
    L[0][1] = -0.5*(1.0/cs + b2*velx)
    L[0][2] = 0.5*b2
    L[1][0] = 1.0-b1
    L[1][1] = b2*velx
    L[1][2] = -b2
    L[2][0] = 0.5*(b1 - velx/cs)
    L[2][1] =0.5*(1.0/cs - b2*velx)
    L[2][2] = 0.5*b2

def convert_to_wspace(x, dens, velx, pres, w1, w2, w3):
    gamma = 1.4
    L =  np.zeros([3,3],dtype=np.float_)
    R =  np.zeros([3,3],dtype=np.float_)

    for i in range(x.size):
        momx = dens[i]*velx[i]
        eneg = 0.5*dens[i]*velx[i]*velx[i] + pres[i]/(gamma-1.0)
        set_LR_matrix(L, R, dens[i], velx[i], pres[i])
        w1[i] = L[0][0]*dens[i] + L[0][1]*momx + L[0][2]*eneg
        w2[i] = L[1][0]*dens[i] + L[1][1]*momx + L[1][2]*eneg
        w3[i] = L[2][0]*dens[i] + L[2][1]*momx + L[2][2]*eneg
    
def read_cfd_data(file, x, dens, velx, pres):
    indx = 0
    for line in file:
        indx += 1

    file.seek(0,0)
    
    x.resize(indx,refcheck=False)
    dens.resize(indx,refcheck=False)
    velx.resize(indx,refcheck=False)
    pres.resize(indx,refcheck=False)

    indx = 0
    for line in file:
        data = line.split()
        x[indx] = data[1]
        dens[indx] = data[2]
        velx[indx] = data[3]
        pres[indx] = data[4]
        indx += 1

file = open(sys.argv[1],'r')
x = np.empty(1)
dens = np.empty(1)
velx = np.empty(1)
pres = np.empty(1)
read_cfd_data(file, x, dens, velx, pres)
file.close()

w1 = np.zeros(x.size)
w2 = np.zeros(x.size)
w3 = np.zeros(x.size)

convert_to_wspace(x, dens, velx, pres, w1, w2, w3)

for i in range(x.size):
    print(x[i], w1[i], w2[i], w3[i])
