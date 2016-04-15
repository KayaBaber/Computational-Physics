'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 3
Problem 1

Hamiltonian Dynamics of a Nonlinear Pendulum
Consider a simple pendulum of length in
gravitational field g. The frequency in the limit of small angles is Ω_0 ≡ radical(g/l) , but do not assume the limit
of small angles for the following calculations.
(a) Start with the Hamiltonian and develop two first order equations for the angle θ and its conjugate
momentum p_θ .

    ((d^2)θ/d(t^2)) = - (g/l)sin(θ)
    
     θ_dot = P_θ/(ml)^2
     P_θ_dot = -mlsin(θ)
    
    
(b) Use a second-order leapfrog algorithm to compute the motion of the pendulum. If we choose a
computational unit of time [T ] = Ω_0^(−1) , then 2π computational time units equals one period in the limit of
small oscillations. Another way to think about it is that we can choose a set of units such that Ω_0 = 1.
Make a graph of phase space trajectories for a variety of initial conditions.




(c) Liouville’s Theorem states that the phase-space volume of a infinitesimally close ensemble of states is
conserved. Demonstrate Liouville’s Theorem by considering an ensemble of closely spaced initial conditions.



'''


