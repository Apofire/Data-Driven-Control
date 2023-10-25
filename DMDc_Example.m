%% DMDc EXAMPLE


%% DEFINE RANDOM SYSTEM

NumStates = 50; % Number of states
NumInputs = 3;  % Number of inputs
NumOutputs = 3; % Number of outputs
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

%% PERFORM DMDc

r = 10;
rprime = r;
[sysR,ATilde,BTilde] = DMDc(X,Xprime,Upsilon,r,rprime,NumStates,NumInputs);

%% PLOTS FOR VERIFICATION
if NumStates <= 5*rprime 
    % Form the original system using C as proposed in [1]
    A = sys.A;
    B = sys.B;
    C = sysR.C;
    D =  0;
    sys1 = ss(A,B,C,D);

    % Plot the singular values in the frequency domain (Generalisation of 
    % Bode plot in MIMO system)
    sigmaplot(sys1,sysR,'r--')
    grid on;
    legend('Original System','System after DMDc')
    xlabel('Frequency',FontSize=20) ; ylabel('Magnitude',FontSize=20)
    title('Singular Values',FontSize=20);
    ax = gca;
    ax.Box;
    ax.LineWidth = 2;
    ax.GridLineStyle = '--';
    ax.FontSize = 20;
    set(gco,'LineWidth',2) 
else
    
    % Plot the singular values in the frequency domain (Generalisation of 
    % Bode plot in MIMO system)
    sigma(sys,sysR)
    grid on;
    legend('Original System','System after DMDc')
    xlabel('Frequency',FontSize=20) ; ylabel('Magnitude',FontSize=20)
    title('Singular Values',FontSize=20);
    ax = gca;
    ax.Box;
    ax.LineWidth = 2;
    ax.GridLineStyle = '--';
    ax.FontSize = 20;

end







