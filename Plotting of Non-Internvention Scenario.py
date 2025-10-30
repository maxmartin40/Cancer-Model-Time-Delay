# Plotting the non-intervention scenario

import numpy as np
from ddeint import ddeint
import matplotlib.pyplot as plt

# Define paramaters from Patient 9 in Wei (2018)
c = 1.25e-7
d = 5.429
f = 0.0956
l = 2.09
g = 0.029
q = 3230
p = 7779
s = 0.0839
h = 2.1e-11
j = 0.0578
k = 2.05e-10
r_1 = 2.14e-4
r_2 = 0.2624
K_T = 2.088
K_L = 1.392
K_N = 1.392
K_C = 1.392
m = 0.4733
beta = 0.0278
gamma = 2.088
alpha = 7.5e8
u = 573.1
w = 1e10-8

# Equilibrium point as starting value
T_0, N_0, L_0 = 0, 1/(f*beta), 0

# Define model
def model(Y, t):
    T, N, L = Y(t)
    D = (d*L**l)/(w+s*T**l+L**l)
    C = 1/beta
    dTdt = T*(1-T)-c*N*T-D*T
    dNdt = C-f*N+(g*T**2*N)/(h+T**2)-p*N*T 
    dLdt = -m*L+(L*j*D**2*T**2)/(k+D**2*T**2)-q*L*T+(r_1*N+r_2*C)*T-u*N*L**2
    return np.array([dTdt, dNdt, dLdt])

# Define time before t=0
def history(t):
    return np.array([T_0, N_0, L_0]) + eps

# Perturbation
eps = np.array([0.01, 0.0, 0.0])

# Evaluation period
t_eval = np.linspace(0, 3, 200)

# Solve the system
sol = ddeint(model, history, t_eval)

# Graph the system
plt.figure(figsize=(8,4))
plt.plot(t_eval, sol[:,0], label="T(t)")
plt.plot(t_eval, sol[:,1], label="L(t)")
plt.plot(t_eval, sol[:,2], label="N(t)")
plt.xlabel("t")
plt.ylabel("x(t), y(t)")
plt.title("Perturbed system near equilibrium")
plt.legend()
plt.grid(True)
plt.show()
