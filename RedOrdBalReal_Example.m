% SAMPLE SCRIPT TO VERIFY THE PERFORMANCE OF BALANCED ORDER REDUCTION

%% DEFINE A CONTINUOUS TIME SYSTEM
% EXAMPLE 1
A = [-.75 1; -.3 -.75];
B = [2; 1];
C = [1 2];
D = 0;
sys_full = ss(A,B,C,D);

% EXAMPLE 2
% n = 30;        % Order of full system
% NumInputs = 2;  % Number of inputs
% NumOutputs = 2; % Number of outputs
% sys_full = rss(n,NumInputs,NumOutputs); % Generate random system

% Parameters for reduced order model
METHOD = 'SVD';      % Method for Balanced realisation
% r = 10;              % Order of reduced model
% TRUNCATION = 'True'; % Parameter for truncation 


%% FIND THE REDUCED ORDER SYSTEM

[~, ~, ~, Wo_hat, ~] = RedOrdBalReal(sys_full,METHOD,TRUNCATION,r); % Find reduced order model

%% IMPULSE RESPONSE SIMULATION 
% Ti = 0;     % Initial time
% Ts = 1e-2;  % Sampling time 
% Tf = 20;    % Final time
% 
% T = Ti:Ts:Tf; % Define the time vector
% 
% figure(2)
% impulse(sys_full,T);
% hold on
% impulse(sys_R,T);
% grid on;




