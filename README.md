# Markov-Chain-model-for-infectious-disease
Markov Chain Monte Carlo simulation for infections disease SEIRDM model

This is a MCMC with semi random transition matrix P. The states are susceptible (S), exposed (E), infected (I), recovered (R), dead (D), and immune (M). The transition probabilites between states susceptible -> exposed, and exposed -> ifected are drawn from uniform distribution at each Monte Carlo iteration with user defined interval: P(S->E) ~ U[a,b] and P(E->I) ~ U[c,d]. Transition probability P(I->R) is from worldometer stats (should be updated). The transition probability P(R->S) is user defined and can / should be changed in order to see how immunity-rate affects to the active time of disease in the population (intuition: if once recovered does not get immunity, they are again susceptible).

Function mc.mc input:
b = single positive integer value: number of Monte Carlo iterations
srate = single numerical value between 0 and 1: probability P(R->S)
pii = size element vector: system state / distribution at the beginning. Must sum into 1. i.e., pii = c(1,0,0,0,0,0) indicates that the system is at susceptible (S) state.
x = four element vector: random limits [a,b] and [c,d], that are use to calculated P(S->E) ~ U[a,b] and P(E->I) ~ U[c,d]. Must be between 0 and 1 and a < b, and c < d.
k = two element vector: constant rates / probabilites for P(I->R) and P(R->S)

Function output:
b x 6 matrix where each row sums to 1 and represent the system probability distribution and Monte Carlo iteration b_i. The system should reach "steady" state in approximately 1000 iterations; however, due to semi random transition matrix P, the steady state does not occur per se; however, law of large numbers apply and the system dampens quite quickly.

The real-world applicability of this model is questionable and the randomizen transition matrix P is another exotic piece. Typically, the transition matrix P is estimated from statistics and different scenarios are estimated in order to calculate some type of confidence intervals. That is, this model here is curioisity and cannot be used to estimate any real-life situation.
