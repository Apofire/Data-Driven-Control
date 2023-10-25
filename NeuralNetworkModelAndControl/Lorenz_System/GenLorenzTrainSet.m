function XTrainSet = GenLorenzTrainSet(t,X0,params,NumTrainingSets)

% Total number of inputs
TrainingSetLength = (length(t) - 1)*NumTrainingSets;

% Preallocate Input and Output 
Input = zeros( TrainingSetLength,length(X0)); 
Output = zeros( TrainingSetLength,length(X0));

% Simulate System for 'NumTrainingSets' times
for j = 1:NumTrainingSets
    [t,y]  = LorenzSys(t,X0,params);
    Input((length(t)-1)*(j-1) + 1:(length(t)-1)*j,:)  = y(1:end-1,:); % Input  = Lorenz states at t
    Output((length(t)-1)*(j-1) + 1:(length(t)-1)*j,:) = y(2:end,:);   % Output = Lorenz states at t + dt
end

% Form Training Datset
XTrainSet.Inputs = Input;
XTrainSet.Outputs = Output;

end