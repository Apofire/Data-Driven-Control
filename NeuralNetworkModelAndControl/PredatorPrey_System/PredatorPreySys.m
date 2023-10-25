function [t,y] = PredatorPreySys(t,X0,params)
% This function computes the numerical solution of the Predator-Prey
% system. 
% INPUTS  -> t (time vector of the form [t0:dt:tf], which is time of 
%               simulation)
%         -> X0 (vector of initial initial conditions X0 = [x0 y0 z0]')
%         -> params (structure of Lorenz system parameters (α,β,γ,δ) )
% OUTPUTS -> t (time vector of simulation)
%         -> y (solution of the Lorenz system for the given initial 
%               conditions)
% =========================================================================

% DEFINE PREDATOR PREY MODEL
PredatorPrey = @(t,x)[ (params.alpha - params.alpha*x(2))*x(1); ...
                       (-params.gamma + params.delta*x(1))*x(2)];

% PredatorPrey = @(t,x) [(2-.5*x(2))*x(1); (-1+.5*x(1))*x(2)];
% SET ODE OPTIONS FOR NUMERICAL INTEGRATION
% ODEInitialisation = odeset('RelTol',1e-8,'AbsTol',1e-10);

% RUN ODE45 TO SIMULATE LORENZ SYSTEM
[t,y] = ode45(PredatorPrey,t,X0);

end