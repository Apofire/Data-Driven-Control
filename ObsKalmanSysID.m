function [H,M] = ObsKalmanSysID(y,u,r)
% The function computes the Markov parameters from the Observer based 
% method as proposed in [1]. The function takes the input and output data 
% and a truncation parameter to compute the observer-based Markov
% parameters, and the observer gain. 
% Inputs  -> y (output data)
%         -> u (input data)
%         -> r (truncation parameter)
% Outputs -> H (Markov parameters)
%         -> M (Observer gain)
% 
% [1] Juang, J. N., Phan, M., Horta, L. G., & Longman, R. W. (1993). 
% Identification of observer/Kalman filter Markov parameters-Theory and 
% experiments. Journal of Guidance, Control, and Dynamics, 16(2), 320-329.
% =========================================================================

% Get the number of inputs and ouputs and total samples for each
[NumInputs,InSamples] = size(u); % (m inputs and l samples of each input)
[NumOutputs,OutSamples] = size(y); % (q outputs and l output samples)

% Insamples must equal Outsamples
if (InSamples ~= OutSamples)
    error ("Number of Input and Output Samples Must Match! ");
end

% Choose p >= 5*r (Truncation index)
p = 5*r;

%% FORM V MATRIX USING (7) IN [1] AND COMPUTE Ybar USING (10)

% ------- Form V matrix --------------------
V1 = zeros((NumOutputs+NumInputs)*p + NumInputs, InSamples); % Initialise
for i = 1:InSamples
    V1(1:NumInputs,i) = u(1:NumInputs,i); % First row is inputs
end
    
for i = 2:p + 1
    for j = 1:InSamples - i + 1
        nu = [u(:,j); y(:,j)]; % Form augmented input
        V1( NumInputs + (i - 2)*(NumOutputs + NumInputs) + 1:NumInputs + ...
           (i - 1)*(NumOutputs + NumInputs),j + i - 1) = nu;
    end
end

% Solve for Ybar (Equation 10)
Ybar = y*pinv(V1, 1e-4);


%% OBTAIN THE MARKOV PARAMETERS Y ( = HANKEL MATRIX)

D = Ybar(:,1:NumInputs); % First column of the Y matrix [Ybar_(-1)]

% For the other entries follow the recursion in (15)

% Extract Ybar^(1) and Ybar^(2) from (11)
% The specific form of Ybar1 and Ybar2 is to be able to form the Markov
% parameters that is suitable to be used in ERA
for i = 1:p
    Ybar1(1:NumOutputs,1:NumInputs,i) = Ybar(:,NumInputs + 1 +  ...
            (NumInputs + NumOutputs)*(i-1): NumInputs + ...
            (NumInputs + NumOutputs)*(i-1) + NumInputs);
    Ybar2(1:NumOutputs,1:NumInputs,i) = Ybar(:,NumInputs + 1 + ...
            (NumInputs + NumOutputs)*(i-1) + NumInputs: NumInputs + ...
            (NumInputs + NumOutputs)*i);
end

% First Entry of Y
Y(:,:,1) = Ybar1(:,:,1) + Ybar2(:,:,1)*D;

% Successive entries from (15)
for j = 2:p
    Y(:,:,j) = Ybar1(:,:,j) + Ybar2(:,:,j)*D;
    for i=1:j-1
        Y(:,:,j) = Y(:,:,j) + Ybar2(:,:,i)*Y(:,:,j-i);
    end
end

% Form the Markov Parameters 
H(:,:,1) = D;
for j = 2:p+1
    H(:,:,j) = Y(:,:,j-1);
end


% Need to compute the gain M
M = 0;

