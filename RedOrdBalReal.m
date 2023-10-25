function [sys_R, sys_B, Wc_hat, Wo_hat, T] = RedOrdBalReal(sys_orig, METHOD, TRUNCATION, r)
% -------------------------------------------------------------------------
% This function computes the balanced realisation of a given system and
% also gives the corresponding truncated reduced order balanced realisation
% based on the rank of the reduced system, as defined by the user. If the
% user does not define a reduced order, then the function plots the
% singular value decomposition (SVD) of the balanced transformed Gramian,
% and asks the user to input a value depending on the plot.  
% Inputs  -> sys_orig (original state space, as state space object)
%            METHOD [optional argument] (specify the method for 
%                                        balanced realisation)
%                       METHOD = 'EIGEN' (Eigen Decomposition)
%                       METHOD = 'SVD'   (SVD of Gramian)
%            TRUNCATION [optional input] (specify if order reduction is 
%                                         required. Default - False)
%                       TRUNCATION = 'True'  : Perform truncation
%                       TRUNCATION = 'False' : Do not Perform truncation
%            r [optional input] (order of the reduced system)
% Outputs -> sys_B (Balanced realisation as a state space object)
%            sys_R (Reduced balanced realisation as a state space object)
%            Wc_hat (The balanced Controllability Gramian)
%            Wo_hat (The balanced Observability Gramian)
%            T (The transformation matrix)
% -------------------------------------------------------------------------

% Check if the METHOD argument has been specified. If not, default it to
% the SVD method.
if nargin < 2
    METHOD = 'SVD';
end

% Check number of inputs. If number of inputs are less than 3, set the
% TRUNCATION argument to 'False' so that balanced truncation is not done.
if nargin < 3
    TRUNCATION = 'False'; 
end

% If the TRUNCATION argument is 'True' but the order of reduction is not
% input, then set a flag. 
if strcmp(TRUNCATION,'True') && nargin < 4
    flag = 1;
else 
    flag = 0;
end

if strcmp(METHOD,'Eigen')
    %% METHOD 1 (USING EIGEN DECOMPOSITION)
    % Compute the controllability and observability Gramians
    Wc = gram(sys_orig,"c"); % Controllability Gramian
    Wo = gram(sys_orig,"o"); % Observability Gramian

    % Find the Eigen decomposition to get the unscaled (normalised) eigen
    % vectors and the desired Diagonal matrix (Sigma^2)
    [T_unscaled, Sigma2] = eig(Wc*Wo);

    % Arrange the eigen values and eigen vectors in decreasing order.
    [~,idx] = sort(diag(Sigma2),'descend');   % Extract indices
    Sigma2 = Sigma2(idx,idx);       % Sorted eigenvalues
    T_unscaled = T_unscaled(:,idx); % Sorted eigenvectors

    % Compute Unscaled Balanced System
    sysb = ss(inv(T_unscaled)*sys_orig.A*T_unscaled, inv(T_unscaled)*sys_orig.B, sys_orig.C*T_unscaled, sys_orig.D);

    % Compute Balanced Gramians
    Wc = gram(sysb,'c');
    Wo = gram(sysb,'o');

    % Compute the scaling
    Sigma_c = diag(Wc)./diag(Wo);
    T = T_unscaled*diag(Sigma_c.^(1/4));

    % % Check if Wc_hat = Wo_hat = Sigma (up to a given threshold)
    % threshold = 1e-3;
    % if abs(Wc_hat - Wo_hat) < threshold*ones(size(Wo_hat))
    %     disp("Succesful Balanced Transformation!")
    % else
    %     error("Something Went Wrong!")
    % end

    % Balanced Realisation
    Abal = inv(T)*sys_orig.A*T;
    Bbal = inv(T)*sys_orig.B;
    Cbal = sys_orig.C*T;
    Dbal = sys_orig.D;

    sys_B = ss(Abal, Bbal, Cbal, Dbal);

    % Compute Balanced Gramians
    Wc_hat = gram(sys_B,'c');
    Wo_hat = gram(sys_B,'o');
else 
    %% METHOD 2 (USING SVD ON GRAMIAN DIRECTLY)
    Wc = gram(sys_orig,"c"); % Controllability Gramian
    Wo = gram(sys_orig,"o"); % Observability Gramian

    [U,S,~] = svd(Wc);  % SVD of Wc

    P1 = U*(S.^(0.5));  % Define the similarity transformation

    Wo_transformed = P1'*Wo*P1;  % Transformed Observability Gramian

    [U2,S_wo,~] = svd(Wo_transformed); % SVD of Wo_transformed

    P2 = U2*(S_wo^(-0.25));  % Balanced transformation

    Wo_hat = P2'*Wo_transformed*P2;  % Compute the balanced realisation
    Wc_hat = Wo_hat;

    T = P1*P2; % Transformation matrix

   % Balanced Realisation
    Abal = inv(T)*sys_orig.A*T;
    Bbal = inv(T)*sys_orig.B;
    Cbal = sys_orig.C*T;
    Dbal = sys_orig.D;

    sys_B = ss(Abal, Bbal, Cbal, Dbal);
end

%% BALANCED TRUNCATION FOR REDUCED ORDER REALISATION
% Perform this only if TRUNCATION argument is 'True'.
% If the user did not enter the order for the reduced realisation (flag = 1)
% then compute the SVD of Sigma (= Wc_hat = Wo_hat), and plot the singular
% values. The user then enters the order based on the singular value spread
% and the algorithm proceeds. In case the user specified the order (flag =
% 0), no need to plot the SVD output.
if strcmp(TRUNCATION,'True')

    [~,S,~] = svd(Wc_hat); % Compute SVD, and extract the singular values

    if flag == 1 % User did not specify the order
        % Plot the singular values
        figure(1)
        plot(diag(S),'o'); xlabel("Index"); ylabel("Singular values"); grid on; box on;

        % Ask user to enter the order of the reduced system
        r = input('Enter the order for the reduced system: \n');
    end

    % Define the transformation matrices to generate the reduced order system
    P = T(:,1:r); % Extract the first r columns of T
    T_inv = inv(T);
    Q = T_inv(1:r,:); % Extract the first r rows of T^(-1)

    % Compute the transformed reduced order system matrices
    ATilde = Q*sys_orig.A*P;
    BTilde = Q*sys_orig.B;
    CTilde = sys_orig.C*P;

    %     % Convert values less than 1e-5 to 0
    %     ATilde(real(ATilde) <= 1e-5) = 0;  ATilde(imag(ATilde) <= 1e-5) = 0;
    %     BTilde(real(BTilde) <= 1e-5) = 0;  BTilde(imag(BTilde) <= 1e-5) = 0;
    %     CTilde(real(CTilde) <= 1e-5) = 0;  CTilde(imag(CTilde) <= 1e-5) = 0;

    % Form the state space object of the reduced order system
    sys_R = ss(ATilde,BTilde,CTilde,sys_orig.D);

else
    sys_R = ss(0,0,0,0); % Define a random state space to avoid error
end

end