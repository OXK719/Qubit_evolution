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
print('completeness:', completeness) #testing completeness is satisfied

i=complex(0,1)

E1 = 3 # energy of state 1
E2 = 2 # energy of state 2

prob_plus = 1/3  #initial probability of qubit being in the plus state
prob_minus = 1 - prob_plus

C_plus = np.sqrt(prob_plus)
C_minus = np.sqrt(prob_minus)

from scipy.optimize import fsolve

def equation(p):
    x, y = p
    return (x+(y*(C_minus/C_plus)), x**2 + y**2 - 1)

x, y =  fsolve(equation, (1, 1))

A_plus = x
A_minus = i*y

hbar = 1.054571

t = np.linspace(0, 20, 200)

print(A_plus, A_minus)
print(C_plus, C_minus)

def Energy(E):
    return np.exp(-i*E*t/hbar)

def Prob_plus(t):
    return abs((((C_plus)**2)*Energy(E1))+(((A_plus)**2)*Energy(E2)))**2

def Prob_minus(t):
    return abs((((C_minus)**2)*Energy(E1))+(((A_minus)**2)*Energy(E2)))**2

plt.figure()
plt.title("Probability of Qubit measured in plus state")
plt.xlabel("Time (10E-34 s)")
plt.ylabel("Normalised Probability")
plt.plot(t, Prob_plus(t))
plt.savefig('Time_evolution_of_plus_state.png')
plt.show()

plt.figure()
plt.title("Probability of Qubit measured in minus state")
plt.xlabel("Time (10E-34 s)")
plt.ylabel("Normalised Probability")
plt.plot(t, Prob_minus(t))
plt.savefig('Time_evolution_of_minus_state.png')
plt.show()

