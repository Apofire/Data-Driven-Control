% SAMPLE SCRIPT TO VERIFY THE PERFORMANCE OF ERA
clear;

%% DEFINE A RANDOM SYSTEM

% load testSys
% [yFull,t] = impulse(sysFull,0:1:(r*5)+1);  
% y = permute(yFull,[2 3 1]);
NumInputs  = 1;   % Number of Inputs
NumOutputs = 1;   % Number of Ouputs
NumStates  = 100; % Number of States

sys  = rss(NumStates, NumOutputs, NumInputs); % Random model
data = impulse(sys,0:1:(100*5)+1);            % Generate impulse response data
data = permute(data,[2,3,1]);                 % Arrange the data


%% Compute reduced order system of order r = 10
r = 10;
[sysR, HSV] = EigenRealAlg(data,NumInputs,NumOutputs,r);

%% PLOTS

[y2,t2] = impulse(sysR);
[y1,t1] = impulse(sys,t2);

% figure(1)
% stairs(y1(:,1,1),'LineWidth',2);
% hold on
% stairs(y2(:,1,1),'LineWidth',1.2);
% grid on
% legend('Full System','Reduced System (order 10)')
% title('Impulse Response from Input 1 to Output 1',FontSize=24)
% xlabel('Time Index',FontSize=20) ; ylabel('Magnitude',FontSize=20)
% grid on
% ax = gca;
% ax.Box;
% ax.LineWidth = 2;
% ax.GridLineStyle = '--';
% ax.FontSize = 20; 
% 
% figure(2)
% stairs(y1(:,1,2),'LineWidth',2);
% hold on
% stairs(y2(:,1,2),'LineWidth',1.2);
% grid on
% legend('Full System','Reduced System (order 10)')
% title('Impulse Response from Input 2 to Output 1',FontSize=24)
% xlabel('Time Index',FontSize=20) ; ylabel('Magnitude',FontSize=20)
% grid on
% ax = gca;
% ax.Box;
% ax.LineWidth = 2;
% ax.GridLineStyle = '--';
% ax.FontSize = 20; 
% 
% figure(3)
% stairs(y1(:,2,1),'LineWidth',2);
% hold on
% stairs(y2(:,2,1),'LineWidth',1.2);
% grid on
% legend('Full System','Reduced System (order 10)')
% title('Impulse Response from Input 1 to Output 2',FontSize=24)
% xlabel('Time Index',FontSize=20) ; ylabel('Magnitude',FontSize=20)
% grid on
% ax = gca;
% ax.Box;
% ax.LineWidth = 2;
% ax.GridLineStyle = '--';
% ax.FontSize = 20; 
% 
% figure(4)
% stairs(y1(:,2,2),'LineWidth',2);
% hold on
% stairs(y2(:,2,2),'LineWidth',1.2);
% grid on
% legend('Full System','Reduced System')
% title('Impulse Response from Input 2 to Output 2',FontSize=24)
% xlabel('Time Index',FontSize=20) ; ylabel('Magnitude',FontSize=20)
% grid on
% ax = gca;
% ax.Box;
% ax.LineWidth = 2;
% ax.GridLineStyle = '--';
% ax.FontSize = 20; 
% 
% figure(5)
% semilogy(HSV,LineWidth=2);
% grid on;
% title('Singuar Values of the Hankel Matrix',FontSize=24)
% grid on
% ax = gca;
% ax.Box;
% ax.LineWidth = 2;
% ax.GridLineStyle = '--';
% ax.FontSize = 20; 
% 
% figure(6)
% plot(0:length(HSV),[0; cumsum(HSV)/sum(HSV)],'k','LineWidth',2)
% hold on, grid on
% plot(r,sum(HSV(1:r))/sum(HSV),'ro','LineWidth',2)
% title('Reduction order in terms of cumulative sum of \sigma',FontSize=24)
% grid on
% ax = gca;
% ax.Box;
% ax.LineWidth = 2;
% ax.GridLineStyle = '--';
% ax.FontSize = 20; 

figure(7)
bodemag(sys,sysR);
grid on;
legend('Original System','Reduced System after ERA')
title('Bode Magnitude Plot',FontSize=20)
xlabel('Frequency',FontSize=20) ; ylabel('Magnitude (dB)',FontSize=20)
ax = gca;
ax.Box;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.FontSize = 20; 
ax.XLim = [1e-3 1e1];