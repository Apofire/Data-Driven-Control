function [sysR,ATilde,BTilde] = DMDc(X,Xprime,Upsilon,r,rprime,NumStates,NumInputs)

% Original system of equation for regression
% Xprime = A*X + B*Upsilon 
% New system of equations
% Xprime = [A B][X; Upsilon] = G*Omega

% Stack the matrices X and Upsilon to form Omega
Omega = vertcat(X,Upsilon); % dimension - (n x p) + m

if NumStates > 5*rprime 
    % Compute SVD of Omega
    [U,S,V] = svd(Omega);

    % Truncate using r
    UTilde = U(:,1:r);
    SigmaTilde = S(1:r,1:r);
    VTilde = V(1:r,:);

    % Compute SVD of Xprime
    [U,S,V] = svd(Xprime);

    % Truncate using rprime
    UHat = U(:,1:rprime);
    SigmaHat = S(1:r,1:r);
    VHat = V(1:r,:);

    % Decompose UTilde = [UTilde1' UTilde2']'
    % where dim(UTilde1) = (n x r) and dim(UTilde2) = (p x r)
    UTilde1 = UTilde(1:NumStates,:);
    UTilde2 = UTilde(NumStates+1:NumStates + NumInputs,:);

    % Compute approximation of G = [Abar, Bbar] from (28) in [1]
    Abar = Xprime*VTilde'*inv(SigmaTilde)*UTilde1';
    Bbar = Xprime*VTilde'*inv(SigmaTilde)*UTilde2';

    % Compute the reduced matrices ATilde BTilde from (33) and (34) in [1]
    ATilde = UHat'*Xprime*VTilde'*inv(SigmaTilde)*UTilde1'*UHat;
    BTilde = UHat'*Xprime*VTilde'*inv(SigmaTilde)*UTilde2';

    % Generate the reduced model as a state sspace object
    % As suggested in [1], CTilde = Uhat for comparison of performance
    sysR = ss(ATilde,BTilde,UHat,zeros(size(UHat,1),NumInputs));
else
    % Compute 'econ' SVD of Omega
    [UTilde,SigmaTilde,VTilde] = svd(Omega,'econ');

    % Compute 'econ' SVD of Xprime
    [UHat,SigmaHat,VHat] = svd(Xprime,'econ');

    % Decompose UTilde = [UTilde1' UTilde2']'
    % where dim(UTilde1) = (n x r) and dim(UTilde2) = (p x r)
    UTilde1 = UTilde(1:NumStates,:);
    UTilde2 = UTilde(NumStates+1:NumStates + NumInputs,:);

    % Compute approximation of G = [Abar, Bbar] from (28) in [1]
    Abar = Xprime*VTilde*inv(SigmaTilde)*UTilde1';
    Bbar = Xprime*VTilde*inv(SigmaTilde)*UTilde2';

    % Compute the reduced matrices ATilde BTilde from (33) and (34) in [1]
    ATilde = UHat'*Xprime*VTilde*inv(SigmaTilde)*UTilde1'*UHat;
    BTilde = UHat'*Xprime*VTilde*inv(SigmaTilde)*UTilde2';

    % Generate the reduced model as a state sspace object
    % As suggested in [1], CTilde = Uhat for comparison of performance
    sysR = ss(ATilde,BTilde,UHat,zeros(size(UHat,1),NumInputs));

end


end