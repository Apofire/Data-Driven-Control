% SAMPLE SCRIPT TO VERIFY THE PERFORMANCE OF OBSERVER BASED KALMAN SYSTEM
% IDENTIFICATION
clear;

%% DEFINE A RANDOM SYSTEM
NumInputs  = 2;   % Number of Inputs
NumOutputs = 2;   % Number of Ouputs
NumStates  = 100; % Number of States

sys  = rss(NumStates, NumOutputs, NumInputs); % Random model
t    = 0:1:(100*5)+1;          % Time vector for impulse response
data = impulse(sys,t);         % Generate impulse response data
data = permute(data,[2,3,1]);  % Arrange the data

%% SETUP SIMULATION FOR OKID
uRandom = randn(NumInputs,length(t));       % Random Input (not an Impusle now)
yRandom = lsim(sys,uRandom,1:length(t))';   % Output from random input
yNoisy  = yRandom + 0.*rand(size(yRandom)); % Noisy Output

%% RUN THE OKID AND ERA 
r = 10;  % Truncation parameter
[H,M] = ObsKalmanSysID(yNoisy,uRandom,r); % Run OKID
[sysR,HSV] = EigenRealAlg(H,NumInputs,NumOutputs,r); % Run ERA 


%% PLOTS

[y2,t2] = impulse(sysR);
[y1,t1] = impulse(sys,t2);

figure(1)
stairs(y1(:,1,1),'LineWidth',2);
hold on
stairs(y2(:,1,1),'LineWidth',1.2);
grid on
legend('Full System','Reduced System')
title('Impulse Response from Input 1 -> Output 1')

figure(2)
stairs(y1(:,1,2),'LineWidth',2);
hold on
stairs(y2(:,1,2),'LineWidth',1.2);
grid on
legend('Full System','Reduced System')
title('Impulse Response from Input 2 -> Output 1')

figure(3)
stairs(y1(:,2,1),'LineWidth',2);
hold on
stairs(y2(:,2,1),'LineWidth',1.2);
grid on
legend('Full System','Reduced System')
title('Impulse Response from Input 1 -> Output 2')

figure(4)
stairs(y1(:,2,2),'LineWidth',2);
hold on
stairs(y2(:,2,2),'LineWidth',1.2);
grid on
legend('Full System','Reduced System')
title('Impulse Response from Input 2 -> Output 2')
