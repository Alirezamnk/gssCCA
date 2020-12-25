function [X_train , X_test , Y_train , Y_test] = Do_Cross_Validation(X , Y , i_CV)

global n_CV_samples
global whole_sample

n = size(X , 1);

% for i_CV = 1 : n_CV
    temp_sample = whole_sample;
    testing_sample = whole_sample(((i_CV - 1) * n_CV_samples/2 + 1) : (i_CV * n_CV_samples/2));
    testing_sample = [testing_sample        testing_sample + n/2];

    temp_sample(testing_sample) = [];   
    training_sample = temp_sample;

    X_train = X(training_sample, :);
    X_test = X(testing_sample, :);

    Y_train = Y(training_sample, :);
    Y_test = Y(testing_sample, :);

%     [X_train , mu_X , C_X] = standardize(X_train);
%     [Y_train , mu_Y , C_Y] = standardize(Y_train);
% 
%     U = X_train * W1;
%     V = Y_train * W2;
%     cor1 = corr(U , V);
% 
%     X_test = standardize(X_test , mu_X , C_X);
%     Y_test = standardize(Y_test , mu_Y , C_Y);
% 
%     U_test = X_test * W1;
%     V_test = Y_test * W2;
% end
