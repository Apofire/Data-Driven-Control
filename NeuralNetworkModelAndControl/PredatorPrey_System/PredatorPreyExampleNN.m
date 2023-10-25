% LORENZ SYSTEM IDENTIFICATION EXAMPLE USING NEURAL NETWORK (FEEDFORWARD)

clear;
clc;

if exist('PredatorPreyTrainedNeuralNet.mat')
    load PredatorPreyTrainedNeuralNet.mat; % Load trained Neural Network 
    flag = 1;
else 
    flag = 0;
end

%% SETUP SIMULATION OF PREDATOR SYSTEM
% PREDATOR SYSTEM PARAMETERS 
params.alpha = 2;    % α
params.beta  = -0.5; % β
params.gamma = 1;    % γ
params.delta = 0.5;  % δ

% SIMULATION PARAMETERS
t0 = 0;             % Initial time of simulation
dt = 0.1;          % Time step
tf = 100;             % Final time of simulation
t  = t0:dt:tf;      % Time vector for simulation
X0 = randi([1 10],2,1); % Random initial conditions

if flag == 0
    %% PREPARE TRAINING SET
    NumTrainingSets = 200; % Number of training sets of data
    XTrain = GenPredatorPreyTrainSet(t,X0,params,NumTrainingSets);

    %% SETUP AND TRAIN NEURAL NETWORK
    NumLayers = 1;  % Number of hidden layers
    SizeLayer = 25; % Size of each layer
    HiddenLayers = SizeLayer*ones(1,NumLayers);

    NeuralNet = feedforwardnet(HiddenLayers); % Define neural netowrk with hidden layers
    NeuralNet.layers{1}.transferFcn = 'poslin';  % RELU
%     NeuralNet.layers{2}.transferFcn = 'radbas';  % RBF
%      NeuralNet.layers{2}.transferFcn = 'softmax';
%     NeuralNet.layers{4}.transferFcn = 'softmax'; % Softmax


    NeuralNet = train(NeuralNet,XTrain.Inputs',XTrain.Outputs');

    save PredatorPreyTrainedNeuralNet NeuralNet;
else 
   disp('Neural Network has already been trained. ');
end


%% TEST THE NEURAL NETWORK 

eg = input("Which example do you wish to test? 1 or 2 \n");

if eg == 1
    X0Test = randi([1 4],2,1);
    [tNet,yTest] = PredatorPreySys(t(1:500),X0Test,params);  % True output

    y0 = zeros(length(X0Test),1);

    yNN(1,:) = X0Test;
    for i = 2:length(tNet)
        y0 = NeuralNet(X0Test);
        yNN(i,:) = y0';
        X0Test = y0;
    end
else
    X0Test = [6;2];
    [tNet,yTest] = PredatorPreySys(t(1:500),X0Test,params);  % True output

    y0 = zeros(length(X0Test),1);

    yNN(1,:) = X0Test;
    for i = 2:length(tNet)
        y0 = NeuralNet(X0Test);
        yNN(i,:) = y0';
        X0Test = y0;
    end
end

%% PLOTS

figure(1)
plot(yTest(:,1),yTest(:,2),'LineWidth',2);
% plot3(X0Test(1),X0Test(2),X0Test(3),'ro');
hold on
plot(yNN(:,1),yNN(:,2),'LineStyle',':','LineWidth',2);
grid on
legend('True Dynamics','NN Estimate')
xlabel('$x_{1}(t)$','Interpreter','latex');
ylabel('$x_{2}(t)$','Interpreter','latex')
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 

figure(2)
subplot(2,1,1)
plot(yTest(:,1),'LineWidth',2); hold on; plot(yNN(:,1),'LineWidth',2);
ylabel('$x_{1}(t)$','Interpreter','latex')
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 
subplot(2,1,2)
plot(yTest(:,2),'LineWidth',2); hold on; plot(yNN(:,2),'LineWidth',2);
ylabel('$x_{2}(t)$','Interpreter','latex')
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 


figure(3)
subplot(2,1,1)
plot(abs(yTest(:,1) - yNN(:,1)),'LineWidth',2);
ylabel('$e_{x_{1}}(t)$','Interpreter','latex','FontSize',20)
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 
subplot(2,1,2)
plot(abs(yTest(:,2) - yNN(:,2)),'LineWidth',2);
ylabel('$e_{x_{2}}(t)$','Interpreter','latex','FontSize',20)
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 

