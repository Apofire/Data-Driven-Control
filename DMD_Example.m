% DMD EXAMPLE
clear all;
clc;

%% DEFINE RANDOM SYSTEM

NumStates = 100; % Number of states
NumInputs = 1;  % Number of inputs
NumOutputs = 1; % Number of outputs
NumSamples = 100; % Number of measurement samples

sys = drss(NumStates,NumOutputs,NumInputs);

Upsilon = randn(NumInputs,NumSamples-1); % Random Input data matrix

%% COLLECT DATA

% Initialise data matrices for states 
X = zeros(NumStates,NumSamples-1);
Xprime = X;

% Generate data matrix X by propagating the states through the system 
% dynamics by forcing the system through random inputs
for idx = 1:NumSamples-1
    X(:,idx+1) = sys.A*X(:,idx) + sys.B*Upsilon(:,idx);
end

% Form the Xprime and X matrices
Xprime(:,1:end) = X(:,2:end); % Last m-1 terms of X = First m-1 terms of Xprime
X = X(:,1:end-1); 

%% PERFORM DMD

r = 10; % Truncation parameter

% Compute DMD 
[Abar, ATilde] = DMD(X,Xprime,r);

% sysR = ss(ATilde,sys.B,sys.C,sys.D);

%% PLOTS FOR VERIFICATION
% Plot the singular values in the frequency domain (Generalisation of
% Bode plot in MIMO system)
% bodemag(sys,sysR)
% grid on;
% legend('Original System','System after DMD')
% 

