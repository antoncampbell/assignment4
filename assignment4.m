%% Assignment 4
% By: Anton Campbell, 100975168
%% Question 1 and 2 - PA 9 and Transient Circuit Simulation
%
% The simulation was run:
MNA_Q1
pause(2);
%%
% 
% *(a)*
%
% The values of the C and G matrices are shown above.
%
% The G matrix is:
%
% $$\begin{array} {cccccc}1 & 0 & 0 & 0 & 0 & 0 \\ -\frac{1}{R_{1}} & \left
% ( \frac{1}{R_{2}}+\frac{1}{R_{1}} \right ) & 0 & 0 & 0 & 1 \\ 0 & 0 & 
% \frac{1}{R_{3}} & 0 & 0 & -1 \\ 0 & 0 & -\frac{\alpha}{R_{3}} & 1 & 0 & 0 
% \\ 0 & 0 & 0 & -\frac{1}{R_{4}} & \left ( \frac{1}{R_{4}}+\frac{1}{R_{O}}
% \right ) & 0 \\ 0 & -1 & 1 & 0 & 0 & 0 \end{array}$$
%
% The C matrix is:
%
% $$\begin{array} {cccccc}0 & 0 & 0 & 0 & 0 & 0\\ 
% -C_{1} & C_{1} & 0 & 0 & 0 & 0\\ 
% 0 & 0 & 0 & 0 & 0 & 0\\ 
% 0 & 0 & 0 & 0 & 0 & 0\\ 
% 0 & 0 & 0 & 0 & 0 & 0\\ 
% 0 & 0 & 0 & 0 & 0 & L_{1} \end{array}$$
%
% The V matrix is:
%
% $$\begin{array}{cccccc} V_{1}  \\ V_{2}  \\V_{3} \\
% V_{4}  \\ V_{5}  \\I_{L} \end{array}$$
%
% *(b)*
%
% Figure 1 shows a plot of the DC sweep. The relationship of $V_{O}$ and 
% $V_{3}$ to $V_{i}$ is linear.
%
% *(c)*
%
% Figure 2 shows the output voltage as a function of the angular frequency.
% Figure 3 shows the gain as a function of the angular frequency. The
% ciruit acts like a low pass filter.
%
% Figure 4 shows a histogram of Gain for pertubation in C. It has a normal
% distribution just like the perturbation in C.
%
% *(d)*
% 
% Figure 5 shows the input and output voltage for a step input. The output
% voltage overshoots and then exponentialy approaches a final votage.
%
% Figure 6 shows the input and output voltage for a sinusoidal input. The
% output voltage is also a sinusoidal signal but it has a higher amplitude and
% a phase shift. There is some irregularity since the input voltage starts
% from a steady state of 0V.
%
% Figure 7 shows the input and output voltage for a Gaussian pulse input.
% The output voltage resembles the Gaussian input but overshoots below 0V
% before gradually returning to 0V.
%
% *(e)*
%
% Figure 8 shows the frequency domain for the step input. The input
% and output both resemble a dirac-delta centred at a frequency of 0Hz. The
% output is a slighty worse dirac-delta.
%
% Figure 9 shows the frequency domain for the sinsoidal input. The input
% and output both resemble two dirac-delta at -33.3Hz and 33.3Hz. The
% output is a slighty worse pair of dirac-deltas.
%
% Figure 10 shows the frequency domain for the Gaussian pulse input. The
% input and output both resemble Gaussian distributions. The output has
% larger magnitude.


%% Question 3 - Circuit with Noise
%
% The simulation was run:
MNA_Q3
pause(2);


%%
% 
% *(a)*
%
% The value of the new C matrix is shown above.
%
% The new C matrix is:
%
% $$\begin{array} {cccccc}0 & 0 & 0 & 0 & 0 & 0\\ 
% -C_{1} & C_{1} & 0 & 0 & 0 & 0\\ 
% 0 & 0 & C_{n} & 0 & 0 & 0\\ 
% 0 & 0 & 0 & 0 & 0 & 0\\ 
% 0 & 0 & 0 & 0 & 0 & 0\\ 
% 0 & 0 & 0 & 0 & 0 & L_{1} \end{array}$$
%
% *(b)*
%
% Figure 11 shows the input and output voltages for added noise. 
%
% *(c)*
%
% Figure 12 shows the frequency domain for added noise. 
%
% *(d)*
%
% Figure 13 shows the input and output voltages for different values of Cn.
% The initial Cn used was 10\mu F which has noise visible on the output
% signal. For the middle subplot, a larger Cn of 1mF which eliminated most
% of the noise. For the bottom subplot, an even larger Cn of 10mF resulted
% in the output signal oscillating.
%
% *(e)*
%
% Figure 14 shows the input and output voltages for different numbers of
% timesteps. The upper subplot has 1000 timesteps. The lower subplot has
% 10 000 timesteps. The subplot with more timestpes appears to be more
% noisy. This is due to the higher sampling rate.

%% Question 4 - Non-linearity
% 


%%
% 
% *(a)*
%
% In order to implement a non-linear output stage voltage, the B matrix
% would need to be used. The new MNA equation that should be used is:
% 
% $$C\frac{dV}{dt}+GV+B\left ( V \right )=F\left ( t \right )$$
%
% *(b)*
% The MNA equation should be rearranged to the form:
%
% $$\left ( \frac{C}{\Delta t} +G\right )V_{n}-\frac{C}{\Delta t}V_{n-1}+-B\left ( V \right )-F\left ( t \right )=0$$
%
% Newton Raphson would need to be solved iteratively. I would create a H
% matrix. This H matrix would be used to iteratively solve for the V
% matrix.
%

