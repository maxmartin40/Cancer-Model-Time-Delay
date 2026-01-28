from sympy import *
import numpy as np
import random
from scipy.optimize import root
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#c, d, l, w, s, f, h, p, j, m, q, k, r_1, r_2, u, B, g = symbols('c d l w s f h p j m q k r_1 r_2 u B g')

# Parameter values from Patient 10
c = 8.32e-8
d = 4.362
f = 0.0956
l = 1.81
g = 0.029
q = 3617
p = 8166
s = 0.512
h = 2.1e-11
j = 0.0578
k = 3.17e-10
r_1 = 1.43e-4
r_2 = 0.175
K_T = 2.088
K_L = 1.392
K_N = 1.392
K_C = 1.392
m = 21.16
B = 0.0186
gamma = 2.088
alpha = 5e8
u = 382.1
w = 1e10-8
C = 1/B
#D = (d*L**l)/(w+s*T**l+L**l)
# dT_dt = T*(1-T) - c*N*T - D*T
# dN_dt = C - f*N + (g*N*T**2)/(h+T**2) - p*N*T
# dL_dt = -m*L + (j*L*D**2*T**2)/(k+D**2*T**2) - q*L*T + (r_1*N+r_2*C)*T - u*N*L**2

'''
Method #1: scipy root
'''
# Define system of ODEs
def F(v):
    T, N, L = v
    dT_dt = T*(1-T) - c*N*T - ((d*L**l)/(w+s*T**l+L**l))*T
    dN_dt = C - f*N + (g*N*T**2)/(h+T**2) - p*N*T
    dL_dt = -m*L + (j*L*((d*L**l)/(w+s*T**l+L**l))**2*T**2)/(k+((d*L**l)/(w+s*T**l+L**l))**2*T**2) - q*L*T + (r_1*N+r_2*C)*T - u*N*L**2

    # out = np.array([dT_dt, dN_dt, dL_dt], dtype=float)
    # out = np.nan_to_num(out, nan=1e6, posinf=1e6, neginf=-1e6)
    # return out
    return np.array([dT_dt, dN_dt, dL_dt])

# Define equation to find equilibria
# Bounds is initial guess range for T, N, L respectively
def find_equilibria(num_starts=1000, bounds=[(3e-4,0.009),(0,1000),(0,1000)]):
    # Initialize list of roots
    roots = []
    # roots.append(["T","N","L"])
    # roots.append(["B","d","m"])
    # Generate random guesses for B, d, m from ranges in Wei
    # B = random.uniform(0.001,1)
    B = 0.0186
    # d = random.uniform(0,7)
    d = 4.36
    m = random.uniform(0,25)
    # Generate random guesses inside range from above
    for i in range(num_starts):
        guess = np.array([
            random.uniform(3e-4,0.009),
            random.uniform(*bounds[1]),
            random.uniform(*bounds[2])
        ])
        # Run a solution to the system
        sol = root(F, guess, tol=1e-8)
        # Determine if there is a duplicate root and if not, add to root list
        if sol.success:
            r = sol.x
            # if not any(np.allclose(r, s, atol=1e-6) for s in roots):
            #     roots.append(r)
            #     roots.append([B,d,m])
            roots.append(r)
            roots.append([B,d,m])
    return roots

# equilibria = find_equilibria()

# print("Equilibria points found:")
# for r in equilibria:
#     print(r)

'''
Use the least squared methods to try to find other equilibrium points.
The system is the same but the RHS has been squared.
'''
# Define system of ODEs
def F_2(v):
    T, N, L = v
    dT_dt = (T*(1-T) - c*N*T - ((d*L**l)/(w+s*T**l+L**l))*T)**2
    dN_dt = (C - f*N + (g*N*T**2)/(h+T**2) - p*N*T)**2
    dL_dt = (-m*L + (j*L*((d*L**l)/(w+s*T**l+L**l))**2*T**2)/(k+((d*L**l)/(w+s*T**l+L**l))**2*T**2) - q*L*T + (r_1*N+r_2*C)*T - u*N*L**2)**2

    # out = np.array([dT_dt, dN_dt, dL_dt], dtype=float)
    # out = np.nan_to_num(out, nan=1e6, posinf=1e6, neginf=-1e6)
    # return out
    return np.array([dT_dt, dN_dt, dL_dt])

# Define equation to find equilibria
# Bounds is initial guess range for T, N, L respectively
def find_equilibria_2(num_starts=1000, bounds=[(0.5,1),(0,1000),(0,1000)]):
    # Initialize list of roots
    roots = []
    T_y = []
    m_x = []
    # roots.append(["T","N","L"])
    # roots.append(["B","d","m"])
    # Generate random guesses for B, d, m from ranges in Wei
    # B = random.uniform(0.001,1)
    B = 0.0186
    # d = random.uniform(0,7)
    d = 4.36
    # m = random.uniform(0,25)
    # Generate random guesses inside range from above
    for i in range(num_starts):
        m = random.uniform(0,25)
        guess = np.array([
            random.uniform(*bounds[0]),
            random.uniform(*bounds[1]),
            random.uniform(*bounds[2])
        ])
        # Run a solution to the system
        sol = root(F_2, guess, tol=1e-8)
        # Determine if there is a duplicate root and if not, add to root list
        if sol.success:
            r = sol.x
            # if not any(np.allclose(r, s, atol=1e-6) for s in roots):
            #     roots.append(r)
            #     roots.append([B,d,m])
            roots.append(r)
            roots.append([B,d,m])
    for i in range(len(roots)):
        if i%2 == 0:
            T_y.append(roots[i][0])
        else:
            m_x.append(roots[i][2])
    return (m_x, T_y)

equilibria = find_equilibria_2()

# print("Equilibria points found:")
# for r in equilibria:
#     print(r)
print(max(equilibria[1]))


'''
Graphing equilibrium points to recreate Figure 2(c).
Takes equilibrium points found above and plots them on a T vs m graph to mirror Figure 2(c) in Wei's paper
'''

T_y = equilibria[1]
m_x = equilibria[0]
plt.plot(m_x,T_y, '.')
plt.xlabel("m")
plt.ylabel("T")
plt.xticks([0,10,20])
plt.yticks([0,6.3e-6,3e-4,1])
plt.title("Recreation of Figure 2(c)")
plt.grid(True)
plt.show()
