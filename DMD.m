function [Abar,ATilde] = DMD(X, Xprime, r)
% This function computes the Dynamic Mode Decompostion (DMD) for an
% unforced system. The estimate of the original system matrix A and the
% corresponding reduced order matrix estimate. The algorithm is followed
% from [1].
% Inputs  -> X (Data of states)
%         -> Xprime (Shifted version of X)
%         -> r (truncation parameter)
% Outputs -> Abar (Approximate estimate of A)
%         -> ATilde (Reduced/truncated matrix A of order r)
% 
% Proctor, J. L., Brunton, S. L., & Kutz, J. N. (2016). Dynamic mode 
% decomposition with control. SIAM Journal on Applied Dynamical Systems,
% 15(1), 142-161.
% =========================================================================

% SVD of X
[U,S,V] = svd(X);

% Extract first r components
U_r = U(:,1:r);
Sigma_r = S(1:r,1:r);
V_r = V(1:r,:);

% Compute the approximate A
Abar = Xprime*V_r'*inv(Sigma_r)*U_r';

% Choose the co-ordinate transformation Px = xTilde, where P = U_r'. Use
% this for the reduced system and extract the corresponding ATilde matrix
ATilde = U_r'*Xprime*V_r'*inv(Sigma_r);


end