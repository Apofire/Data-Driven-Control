function [sysR,HSV] = EigenRealAlg(data, NumInputs, NumOutputs,r)
% The function computes the reduced system using the eigen system 
% realisation algorithm. The data is converted to a Hankel matrix and a
% shifted version of it. The Hankel matrix is decomposed using SVD and the
% reduced system is computed from this decomposition.
% Inputs  -> data (data in the form of [Input:Output:Samples])
%         -> NumInputs (Number of Inputs)
%         -> NumOutputs (Number of Outputs)
%         -> r (order of the reduced model)
% Outputs -> sysR (Reduced State Space system [A_r,B_r,C_r,D_r])
%         -> HSV (singular values of the Hankel matrix)
% =========================================================================

% Fix the dimension of the Hankel matrix (square)
HankelRows = floor( (length(data) - 1)/2);
HankelCols = HankelRows;

% Separate the first entry for each input-output data and set it to the D
% matirx. The remaining data is used in the Hankel matrix
for i = 1:NumInputs
    for j = 1:NumOutputs
        D(i,j) = data(i,j,1);
        y(i,j,:) = data(i,j,2:end); 
    end
end

assert(length(y(:,1,1)) == NumOutputs);
assert(length(y(1,:,1)) == NumInputs);
assert(length(y(1,1,:)) >= HankelCols + HankelRows);

% Form the Hankel matrix and the shifted version of it
for i=1:HankelRows
    for j=1:HankelCols
        for q=1:NumOutputs
            for p=1:NumInputs
               H(NumOutputs*i-NumOutputs+q,NumInputs*j-NumInputs+p) = y(q,p,i+j-1);
               Hprime(NumOutputs*i-NumOutputs+q,NumInputs*j-NumInputs+p) = y(q,p,i+j);
            end
        end
    end
end


%% ERA ALGORITHM 
[U,SigmaH,V] = svd(H,'econ'); % Compute SVD of Hankel Matrix
SigmaH_r = SigmaH(1:r,1:r); % Reduced order Singular values
U_r = U(:,1:r); % Reduced U matrix
V_r = V(:,1:r); % Reduced V matrix

A_r = SigmaH_r^(-0.5)*U_r'*Hprime*V_r*SigmaH_r^(-0.5); % Reduced A matrix
B_r = SigmaH_r^(-0.5)*U_r'*H(:,1:NumInputs); % Reduced B matrix
C_r = H(1:NumOutputs,:)*V_r*SigmaH_r^(-0.5); % Reduced C matrix

% Form the reduced system
sysR = ss(A_r,B_r,C_r,D,-1);

% Extract the singular values of the Hankel matrix
HSV = diag(SigmaH);

end