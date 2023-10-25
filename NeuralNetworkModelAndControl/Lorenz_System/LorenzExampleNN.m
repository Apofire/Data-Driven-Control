% LORENZ SYSTEM IDENTIFICATION EXAMPLE USING NEURAL NETWORK (FEEDFORWARD)

clear;
clc;

if exist('LorenzTrainedNeuralNet.mat')
    load LorenzTrainedNeuralNet.mat; % Load trained Neural Network 
    flag = 1;
else 
    flag = 0;
end

%% SETUP SIMULATION OF LORENZ SYSTEM
% LORENZ SYSTEM PARAMETERS
params.sigma = 10;  % σ
params.rho   = 28;  % ρ
params.beta  = 8/3; % β

% SIMULATION PARAMETERS
t0 = 0;             % Initial time of simulation
dt = 0.01;          % Time step
tf = 8;             % Final time of simulation
t  = t0:dt:tf;      % Time vector for simulation
X0 = 30*(randn(3,1) - 0.5); % Random initial conditions

if flag == 0
    %% PREPARE TRAINING SET
    NumTrainingSets = 100; % Number of training sets of data
    XTrain = GenLorenzTrainSet(t,X0,params,NumTrainingSets);

    %% SETUP AND TRAIN NEURAL NETWORK
    NumLayers = 3;  % Number of hidden layers
    SizeLayer = 10; % Size of each layer
    HiddenLayers = SizeLayer*ones(1,NumLayers);

    NeuralNet = feedforwardnet(HiddenLayers); % Define neural netowrk with hidden layers
    NeuralNet.layers{1}.transferFcn = 'logsig';  % Lograithmic Sigmoid
    NeuralNet.layers{2}.transferFcn = 'radbas';  % RBF
    NeuralNet.layers{3}.transferFcn = 'purelin'; % Pure Linear


    NeuralNet = train(NeuralNet,XTrain.Inputs',XTrain.Outputs');

    save LorenzTrainedNeuralNet NeuralNet;
else 
   disp('Neural Network has already been trained. ');
end


%% TEST THE NEURAL NETWORK 

eg = input("Which example do you wish to test? 1 or 2 \n");

if eg == 1
    X0Test = 30*(randn(3,1) - 0.5);
    [t,yTest] = LorenzSys(t,X0Test,params);  % True output

    y0 = zeros(length(X0Test),1);

    yNN(1,:) = X0Test;
    for i = 2:length(t)
        y0 = NeuralNet(X0Test);
        yNN(i,:) = y0';
        X0Test = y0;
    end
else
    X0Test = 50*(randn(3,1) - 0.5);
    [t,yTest] = LorenzSys(t,X0Test,params);  % True output

    y0 = zeros(length(X0Test),1);

    yNN(1,:) = X0Test;
    for i = 2:length(t)
        y0 = NeuralNet(X0Test);
        yNN(i,:) = y0';
        X0Test = y0;
    end
end

%% PLOTS

figure(1)
plot3(yTest(:,1),yTest(:,2),yTest(:,3),'LineWidth',2);
% plot3(X0Test(1),X0Test(2),X0Test(3),'ro');
hold on
plot3(yNN(:,1),yNN(:,2),yNN(:,3),'LineStyle',':','LineWidth',2);
grid on
legend('True Dynamics','NN Estimate')
xlabel('x(t)','Interpreter','latex');
ylabel('y(t)','Interpreter','latex')
zlabel('z(t)','Interpreter','latex')
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 

figure(2)
subplot(3,1,1)
plot(t,yTest(:,1),'LineWidth',2); hold on; plot(t,yNN(:,1),'LineWidth',2);
ylabel('x(t)','Interpreter','latex')
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 
subplot(3,1,2)
plot(t,yTest(:,2),'LineWidth',2); hold on; plot(t,yNN(:,2),'LineWidth',2);
ylabel('y(t)','Interpreter','latex')
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 
subplot(3,1,3)
plot(t,yTest(:,3),'LineWidth',2); hold on; plot(t,yNN(:,3),'LineWidth',2);
ylabel('z(t)','Interpreter','latex')
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 

figure(3)
subplot(3,1,1)
plot(t,abs(yTest(:,1) - yNN(:,1)),'LineWidth',2);
ylabel('$e_{x}(t)$','Interpreter','latex','FontSize',20)
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 
subplot(3,1,2)
plot(t,abs(yTest(:,2) - yNN(:,2)),'LineWidth',2);
ylabel('$e_{y}(t)$','Interpreter','latex','FontSize',20)
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 
subplot(3,1,3)
plot(t,abs(yTest(:,3) - yNN(:,3)),'LineWidth',2);
ylabel('$e_{z}(t)$','Interpreter','latex','FontSize',20)
grid on
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 

