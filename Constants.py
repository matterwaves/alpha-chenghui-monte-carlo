import numpy as np
import matplotlib.pyplot as plt
import os

c0 = 2.99792458e+8
mu0 = 4*np.pi*1e-7
epsilon0 = 8.854187871e-12
hbar = 1.054571596e-34
e = 1.602176462e-19
muB = 9.274009994e-24
me = 9.10938188e-31
kB = 1.3806503e-23
g0 = 9.79958
gamma0 = 2.25e-6
G = 6.67408e-11
M = 2.20694650e-25
f_ref = 3.517309021052292e+14
detuning = 14e+9
f_L = f_ref+detuning+80e+6+90e+6
k0 = 2*np.pi*f_L/c0
omega_r = hbar*k0**2/2/M
omega_r/2/np.pi
n = 5
T = 0.1
v0 = -n*k0*hbar/M
phi = 2*n*k0*g0*T**2*(1+gamma0*(7/12*g0*T**2-T*v0-n*k0*hbar/M*T))

def conversion(quantity, dimension, unit = np.array([2/k0, 1/omega_r])):
    for i in range(len(dimension)):
        quantity = quantity/unit[i]**dimension[i]
    return quantity

g = conversion(g0,[1,-2])
gamma = conversion(gamma0,[0,-2])
omega0 = conversion(2*np.pi*f_L,[0,-1])
c = conversion(c0, [1,-1])
