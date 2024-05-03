#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import math

def read_cfd_data(file, x, dens, velx, pres, mach, shock):
    indx = 0
    for line in file:
        indx += 1

    file.seek(0,0)
    
    x.resize(indx,refcheck=False)
    dens.resize(indx,refcheck=False)
    velx.resize(indx,refcheck=False)
    pres.resize(indx,refcheck=False)
    mach.resize(indx,refcheck=False)
    shock.resize(indx,refcheck=False)

    indx = 0
    for line in file:
        data = line.split()
        x[indx] = data[1]
        dens[indx] = data[2]
        velx[indx] = data[3]
        pres[indx] = data[4]
        mach[indx] = data[5]
        shock[indx] = data[6]
        indx += 1

print(len(sys.argv))

file = open(sys.argv[1],'r')
x = np.empty(1)
dens = np.empty(1)
velx = np.empty(1)
pres = np.empty(1)
mach = np.empty(1)
shock = np.empty(1)
M = np.empty(1)
read_cfd_data(file, x, dens, velx, pres, mach, shock)
file.close()

datasize=np.size(pres)
M.resize(datasize,refcheck=False)

M = np.sqrt(1.6666*pres/dens)

if len(sys.argv)==3 :
    file_ref = open(sys.argv[2],'r')
    x_ref = np.empty(1)
    dens_ref = np.empty(1)
    velx_ref = np.empty(1)
    pres_ref = np.empty(1)
    mach_ref = np.empty(1)
    shock_ref = np.empty(1)
    read_cfd_data(file_ref, x_ref, dens_ref, velx_ref, pres_ref, mach_ref, shock_ref)
    file_ref.close()

if len(sys.argv)==3 :
    dens_max = max(dens_ref)
    dens_min = min(dens_ref)
else:
    dens_max = max(dens)
    dens_min = min(dens)

if np.log10(dens_max/dens_min)>2.0 :
    log_dens = 1
    dens_max = np.log10(dens_max)
    dens_min = np.log10(dens_min)
    dens = np.log10(dens)
    dens_ref = np.log10(dens_ref)
else:
    log_dens = 0

dens_max += (dens_max-dens_min)*0.1
dens_min -= (dens_max-dens_min)*0.1

print(dens_max, dens_min);

if len(sys.argv)==3 :
    velx_max = max(max(velx_ref),max(M))
    velx_min = min(min(velx_ref),min(M))
else:
    velx_max = max(velx)
    velx_min = min(velx)

velx_max += (velx_max-velx_min)*0.1
velx_min -= (velx_max-velx_min)*0.1

print(velx_max, velx_min);

if len(sys.argv)==3 :
    pres_max = max(pres_ref)
    pres_min = min(pres_ref)
else:
    pres_max = max(pres)
    pres_min = min(pres)

if np.log10(pres_max/pres_min)>2.0 :
    log_pres = 1
    pres_max = np.log10(pres_max)
    pres_min = np.log10(pres_min)
    pres = np.log10(pres)
    pres_ref = np.log10(pres_ref)
else:
    log_pres = 0

pres_max += (pres_max-pres_min)*0.1
pres_min -= (pres_max-pres_min)*0.1


print(pres_max, pres_min);

xmin = min(x)
xmax = max(x)

plt.plot(x, dens, color='b', marker='s', ms=5, ls='None')
if len(sys.argv)==3 :
    plt.plot(x_ref, dens_ref, color='k')
plt.xlabel('$x$', fontsize=20)
plt.ylabel('density', fontsize=20)
xaxis=[xmin,xmax]
yaxis=[dens_min, dens_max]
plt.axis(xaxis+yaxis)
plt.grid(which='both', ls='dashed', c='0.75')
plt.axes().set_axisbelow(True)
plt.savefig('dens.png', bbox_inches='tight')

plt.clf()

plt.plot(x, velx, color='b', marker='s', ms=5, ls='None')
plt.plot(x, M, color='k', marker='s', ms=5, ls='None')
if len(sys.argv)==3 :
    plt.plot(x_ref, velx_ref, color='k')
plt.xlabel('$x$', fontsize=20)
plt.ylabel('velocity', fontsize=20)
xaxis=[xmin,xmax]
yaxis=[velx_min, velx_max]
plt.axis(xaxis+yaxis)
plt.grid(which='both', ls='dashed', c='0.75')
plt.axes().set_axisbelow(True)
plt.savefig('velx.png', bbox_inches='tight')

plt.clf()

plt.plot(x, pres, color='b', marker='s', ms=5, ls='None')
if len(sys.argv)==3 :
    plt.plot(x_ref, pres_ref, color='k')
plt.xlabel('$x$', fontsize=20)
plt.ylabel('pressure', fontsize=20)
xaxis=[xmin,xmax]
yaxis=[pres_min, pres_max]
plt.axis(xaxis+yaxis)
plt.grid(which='both', ls='dashed', c='0.75')
plt.axes().set_axisbelow(True)
plt.savefig('pres.png', bbox_inches='tight')


