function [t,y] = LorenzSys(t,X0,params)
% This function computes the numerical solution of the Lorenz dynamical
% system. 
% INPUTS  -> t (time vector of the form [t0:dt:tf], which is time of 
%               simulation)
%         -> X0 (vector of initial initial conditions X0 = [x0 y0 z0]')
%         -> params (structure of Lorenz system parameters (σ,ρ,β) )
% OUTPUTS -> t (time vector of simulation)
%         -> y (solution of the Lorenz system for the given initial 
%               conditions)
% =========================================================================

% DEFINE LORENZ DYNAMICAL SYSTEM
Lorenz = @(t,x) ([params.sigma*(x(2) - x(1))     ; ...
                  x(1)*(params.rho - x(3)) - x(2); ...
                  x(1)*x(2) - params.beta*x(3)     ...
                 ]);


% SET ODE OPTIONS FOR NUMERICAL INTEGRATION
ODEInitialisation = odeset('RelTol',1e-8,'AbsTol',1e-10);

% RUN ODE45 TO SIMULATE LORENZ SYSTEM
[t,y] = ode45(Lorenz,t,X0,ODEInitialisation);

end