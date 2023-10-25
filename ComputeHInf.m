function [Hinf] = ComputeHInf(sys)
% -------------------------------------------------------------------------
% This function computes the H-infinity norm of a system defined as a state
% space model. This function is based on the bisection method as proposed 
% in [1]. For this to work, A needs to be a stability matrix. 
% Inputs  -> sys (state space object of the system)
% Outputs -> Hinf (H-infinity norm of the system)

% [1] Bruinsma, N. A., & Steinbuch, M. (1990). A fast algorithm to compute 
% the H∞-norm of a transfer function matrix. Systems & Control Letters, 
% 14(4), 287-293.
% -------------------------------------------------------------------------

% Extract the system matrices
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% Check if A is stable. If not stable, the algorithm cannot compute
% H_infinity norm. For this, check if Re(λ(A)) < 0. Equivalently, if 
% max{Re(λ(Α))} > 0, then A is unstable
eigA = diag(reig(A));
if max(real(eigA)) > 0 
    error("The matrix A is not a stability/Hurwitz matrix! Cannot compute norm");
end

%% BISECTION ALGORITHM

%  ----- Find the upper and lower bounds for the bisection algorithm -----
[~,SD,~] = svd(D);        % Compute SVD of matrix D
s_max_D  = SD(1);         % Store the largest singular value
Wc       = gram(sys,'c'); % Compute controllability Gramian
Wo       = gram(sys,'o'); % Compute observability Gramian
Hankel   = Wc*Wo;         % Compute Hankel matrix
n        = size(Hankel,1);% Size of Hankel matrix
H1       = sqrt(trace(Hankel/n)); % Max singular value of Hankel
H_all    = sqrt(n*trace(Hankel)); % Sum of Hankel singular values
rLB = max(s_max_D, H1);   % Lower bound 
rUB = s_max_D + 2*H_all;  % Upper bound

threshold = 1e-6; % Threshold for algorithm

while 2*(rUB - rLB) > threshold  % Loop until threshold is achieved   
    r   = (rLB + rUB)/2;      % Currrent iterate

    % Compute the Hamiltonian matrix (M_r)
    R = (r^2)*eye(size(D'*D)) - D'*D; % Temporary variable
    I = eye(size(D*inv(R)*D'));       % Temporary variable

    M_r = [A + B*inv(R)*D'*C, B*inv(R)*B'; 
           -C'*(I + D*inv(R)*D')*C, -(A + B*inv(R)*D'*C)']; 

    [~,Evals] = eig(M_r);  % Compute the eigenvalues of M_r
    Evals = diag(Evals);   % Convert to column vector

    % Check if any eigenvalue is purely imaginary
    if min( abs(real(Evals)) ) < 1.0e-5 
        rLB = r ;  % Set lower bound to current iterate
    else 
        rUB = r ;  % Set upper bound to current iterate
    end
   
end

Hinf = r; % Return the current iterate
end