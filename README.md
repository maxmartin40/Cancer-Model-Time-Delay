# Cancer-Model-Time-Delay

Wei (2018) investigated a model of cancer tumors that concluded that applying chemotherapy to a microscopic tumor could lead the tumor to grow into a lethal size where it would have remained non-lethal in a non-intervention scenario. This research prompts the question that if a time-delay in any part of the treatment (natural immune system response, immunotherapy, or chemotherapy) could be advantageous in the treatment of cancerous tumors.

The system of ODEs that will be investigated are defined below as in Wei (2018):

$$
\dot{T} = T(1-T)-cNT-DT-K_T(1-e^{-M})T
$$
$$
\dot{N} = C-fN+\frac{gNT^2}{h+T^2}-pNT-K_N(1-e^{-M})N
$$
$$
\dot{L} = -mL+\frac{jLD^2T^2}{k+D^2T^2}-qLT+(r_1N+r_2C)T-uNL^2 -K_L(1-e^{-M})L+v_L(t)
$$
$$
\dot{C} = 1-\beta C-K_C(1-e^{-M})C
$$
$$
\dot{M} = -\gamma M+v_M(t)
$$
$$
D = \frac{dL^l}{w+sT^l+L^l}
$$

The code in this repository includes plotting the system numerically, finding equilibrium points, and bifurcation analysis.

## References
Please see citation file for all references.
