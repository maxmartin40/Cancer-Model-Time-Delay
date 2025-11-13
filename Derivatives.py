from sympy import *

# Define symbols
T_0, T_1, N_0, N_1, L_0, L_1, C_0, C_1, epsilon, c, d, beta  = symbols('T_0 T_1 N_0 N_1 L_0 L_1 C_0 C_1 epsilon c d beta')
f, g, q, k, j, m, r_1, r_2, u, w, s, C_t, C_1_t, h, p, l = symbols('f g q k j m r_1 r_2 u w s C_t C_1_t h p l')

# Define equations
D = (d*(L_0+epsilon*L_1)**l)/(w+s*(T_0+epsilon*T_1)**l+(L_0+epsilon*L_1)**l)
T = (T_0+epsilon*T_1)*(1-(T_0+epsilon*T_1))-c*(N_0+epsilon*N_1)*(T_0+epsilon*T_1)-D*(T_0+epsilon*T_1)
N = C_t + epsilon*C_1_t - f*(N_0+epsilon*N_1) + (g*(N_0+epsilon*N_1)*(T_0+epsilon*T_1)**2)/(h+(T_0+epsilon*T_1)**2) - p*(N_0+epsilon*N_1)*(T_0+epsilon*T_1)
L = -m*(L_0+epsilon*L_1)+(j*(L_0+epsilon*L_1)*D**2*(T_0+epsilon*T_1)**2)/(k+D**2*(T_0+epsilon*T_1)**2) - q*(L_0+epsilon*L_1)*(T_0+epsilon*T_1) +(r_1*(N_0+epsilon*N_1)+r_2*(C_0+epsilon*C_1))*(T_0+epsilon*T_1)-u*(N_0+epsilon*N_1)*(L_0+epsilon*L_1)**2
C = 1-beta*(C_0+epsilon*C_1)

# Find first term in Taylor Series
dT0dt = T.subs(epsilon,0)
dN0dt = N.subs(epsilon,0)
dL0dt = L.subs(epsilon,0)
dC0dt = C.subs(epsilon,0)

# Find second term in Taylor Series
dT0dtde = diff(T, epsilon)
dN0dtde = diff(N, epsilon)
dL0dtde = diff(L, epsilon)
dC0dtde = diff(C, epsilon)

dT1dt = dT0dtde.subs(epsilon, 0)
dN1dt = dN0dtde.subs(epsilon, 0)
dL1dt = dL0dtde.subs(epsilon, 0)
dC1dt = dC0dtde.subs(epsilon, 0)
