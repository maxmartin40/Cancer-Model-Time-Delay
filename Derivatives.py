from sympy import *

# Define symbols
T_0, T_1, N_0, N_1, L_0, L_1, C_0, C_1, epsilon, c, d, beta  = symbols('T_0 T_1 N_0 N_1 L_0 L_1 C_0 C_1 epsilon c d beta')
f, g, q, k, j, m, r_1, r_2, u, w, s, C_t, C_1_t, h, p, l, t_0, lamb = symbols('f g q k j m r_1 r_2 u w s C_t C_1_t h p l, t_0 lamb')

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

# Define symbols for matrix
A_11, A_13, A_21, A_22, A_31, A_32, A_33 = symbols('A_11 A_13 A_21 A_22 A_31 A_32 A_33')

# Define matrix entries
A_11_sub = ((-l*s*(T_0**l)*L_0)/(T_0*L_0))/(L_0*l +s*T_0**l +w)**2 - (d*L_0**l)/(L_0*l +s*T_0**l +w) -c*N_0 - T_0*(1-T_0)
A_13_sub = ((l*T_0*L_0**l)/(T_0*L_0))/(L_0*l +s*T_0**l +w)**2 - (d*l*T_0*L_0**l)/(L_0*l +s*T_0**l +w)
A_21_sub = (2*g*N_0*T_0**3)/(T_0+h)**2 + (2*g*N_0*T_0)/(T_0**2+h) - 1
A_22_sub = (2*g*T-0)/(T_0**2+h) - p*T_0 -f
A_31_sub = (L_0**(2*l+1)*T_0**2*d**2*j*((-2*t_0**l*s*l*L_0)/(T_0*L_0))/((L_0**(2*l)*T_0**2*d**2)/(L_0*l +s*T_0**l +w)**2 +k)*(L_0*l +s*T_0**l +w)**3) + (2*L_0**(2*l)*T_0*d**2*j)/(((L_0**(2*l)*T_0**2*d**2)/(L_0*l +s*T_0**l +w)**2 +k)*(L_0*l +s*T_0**l +w)**2) -L_0*q+C_0*r_2+N_0*r_1 + ((L_0**(4*l+1)*T_0**4*d**4*j*((2*T_0**l*l*s*L_0)/(T_0*L_0))/(L_0*l +s*T_0**l +w)**3))/(((L_0**(2*l)*T_0**2*d**2)/(L_0*l +s*T_0**l +w)**2 +k)**2*(L_0*l +s*T_0**l +w)**2) - (((2*L_0**(4*l)*T_0**3*d**4*j/(L_0*l +s*T_0**l +w)**2)/(((L_0**(2*l)*T_0**2*d**2)/(L_0*l +s*T_0**l +w)**2 +k)**2)*(L_0*l +s*T_0**l +w)**2))
A_32_sub = -u*L_0**2 + T_0*r_1
A_33_sub = (L_0**(2*l+1)*T_0**2*d**2*j*(-2*L_0**l*T_0*l)/(T_0*L_0))/((L_0**(2*l)*T_0**2*d**2)/((L_0*l +s*T_0**l +w)**2+k)*(L_0*l +s*T_0**l +w)**3) -2*L_0*N_0*u + (2*L_0**(2*l)*T_0**2*d**2*j*l)/(((L_0*l +s*T_0**l +w)**2+k)*(L_0*l +s*T_0**l +w)**2) + (L_0**(2*l)*T_0**2*d**2*j)/(((L_0*l +s*T_0**l +w)**2+k)*(L_0*l +s*T_0**l +w)**2) -T_0*q -m + ((L_0**(4*l+1)*T_0**4*d**4*j*(2*L_0**l*l*T_0)/(T_0*L_0))/(L_0*l +s*T_0**l +w)**3)/(((L_0*l +s*T_0**l +w)**2+k)**2*(L_0*l +s*T_0**l +w)**2) - ((2*L_0**(4*l+1)*T_0**4*d**4*l*j)/(L_0*(L_0*l +s*T_0**l +w)**2))/(((L_0*l +s*T_0**l +w)**2+k)**2*(L_0*l +s*T_0**l +w)**2)

A = Matrix([[A_11-lamb,c*T_0,A_13,0],[A_21,A_22-lamb,0,exp(-lamb*t_0)],[A_31,A_32,A_33-lamb,T_0*r_2],[0,0,0,-beta-lamb]])

# Print determinant
char_eq = A.det()
pprint(char_eq)

characteristic_equation = (-beta-lamb)*(A_11*A_22*A_33-A_11*A_22*lamb-A_11*A_33*lamb+A_11*lamb**2+A_13*A_21*A_32-A_13*A_22*A_31*lamb-A_21*A_33*T_0*c+A_21*T_0*c*lamb-A_22*A_33*lamb+A_22*lamb**2+A_33*lamb**2-lamb**3)
full_char_eq = characteristic_equation.subs([(A_11,A_11_sub),(A_13,A_13_sub),(A_21,A_21_sub),(A_22,A_22_sub),(A_31,A_31_sub),(A_32,A_32_sub),(A_33,A_33_sub)])
#pprint(full_char_eq)
