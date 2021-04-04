#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:16:32 2020

@author: oliverk
"""

import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

#Qobj([[1],[2],[3],[4],[5]]) creates a 5x1 ket vector 

#x = np.array([[1, 2, 3, 4, 5]])
# Qobj(x) creates a 1x5 bra vector

plus = qt.basis(2,0)
ket_plus = plus.dag()
minus = qt.basis(2,1)
ket_minus = minus.dag()
print(plus, minus)
print(ket_plus, ket_minus)

#qt.sigmaz * spin , acting on spin yields +1 eigenvalue for vector 'up' and -1 for 'down'

completeness = (plus*ket_plus) + (minus*ket_minus)
print(completeness) #testing completeness is satisfied

E1 = 3 # energy of state 1
E2 = 2 # energy of state 2

prob_plus = 1/3
prob_minus = 1 - prob_plus

C_plus = np.sqrt(prob_plus)
C_minus = np.sqrt(prob_minus)

from scipy.optimize import fsolve

def equation(p):
    x, y = p
    return (x+(y*(C_minus/C_plus)), x**2 + y**2 - 1)

x, y =  fsolve(equation, (1, 1))

A_plus = x
A_mins = y

modsquare_C_plus = (C_plus)**2
modsquare_A_plus = (A_plus)**2
hbar = 1.054571

one = 1
i = 1
z=complex(one,i)

t = np.linspace(0, 20, 200)

energy1 = np.exp(-(z-1)*E1*t/hbar)
energy1star = np.exp((z-1)*E1*t/hbar)
energy2 = np.exp(-(z-1)*E2*t/hbar)
energy2star = np.exp((z-1)*E2*t/hbar)

def Probability_plus(t):
    return ((modsquare_C_plus*energy1)+(modsquare_A_plus*energy2))*(modsquare_C_plus*energy1star+(modsquare_A_plus*energy2star))

plt.figure()
plt.title("Time evolution of the probability density of a two-state system")
plt.xlabel("Time (10E-34 s)")
plt.ylabel("Normalised Probability")
plt.plot(t, Probability_plus(t))
plt.savefig('Time_evolution_of_two_state_system.png')
plt.show()

# ALternatively one can express the evolution of the probability density as

w = (E1-E2)/hbar

def Probability_plus2(t):
    return (modsquare_C_plus**2) + (modsquare_A_plus**2) + (2*modsquare_C_plus*modsquare_A_plus*np.cos(w*t))

plt.figure()
plt.title("Time evolution of the probability density of a two-state system")
plt.xlabel("Time (10E-34 s)")
plt.ylabel("Normalised Probability")
plt.plot(t, Probability_plus2(t))
plt.show()
