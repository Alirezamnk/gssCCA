function [opt_a1, opt_a2] = Find_sub_opt_alphas(X, Y, indices, lam_1, lam_2, alpha_1, alpha_2, Omegu, Omegv, ...
                posu, posv, maxit, iterS)

global n_CV
global whole_sample

[X_temp] = standardize(X);
[Y_temp] = standardize(Y);
[startu, startv] = svd_using_QR(X_temp, Y_temp, 1);
clear X_temp
clear Y_temp

K = 1;
N_a1 = length(alpha_1); 
N_a2 = length(alpha_2); 

delta_corr = zeros(N_a1, N_a2);

Lu = zeros(N_a1, 1);
Lv = zeros(N_a2, 1);

for i_a1 = 1:N_a1
    for i_a2 = 1:N_a2

        for i_CV = 1 : n_CV             
            temp_sample = whole_sample;
            testing_sample = whole_sample(indices == i_CV); 
            temp_sample(testing_sample) = [];   
            training_sample = temp_sample;
    
            X_train = X(training_sample, :);
            X_test = X(testing_sample, :);

            Y_train = Y(training_sample, :);
            Y_test = Y(testing_sample, :);

            [X_train , mu_X , C_X] = standardize(X_train);
            [Y_train , mu_Y , C_Y] = standardize(Y_train);

            [W1, W2, ~, ~, ~, ~, ~, ~, ~ , Lu(i_a1) , Lv(i_a2)] = sfpca_nested_bic(X_train, Y_train, K, lam_1, lam_2, ...
                alpha_1(i_a1), alpha_2(i_a2), Omegu, Omegv, startu, startv, posu, posv, maxit, iterS , Lu(i_a1) , Lv(i_a2), 0);

            U = X_train * W1;
            V = Y_train * W2;
            cor_train = corr(U , V);

            X_test = standardize(X_test , mu_X , C_X);
            Y_test = standardize(Y_test , mu_Y , C_Y);

            U_test = X_test * W1;
            V_test = Y_test * W2;
            cor_test = corr(U_test , V_test);
            
            delta_corr(i_a1 , i_a2) = delta_corr(i_a1 , i_a2) + abs(cor_train - cor_test) / n_CV;
        end        
    end
end

Best_index = delta_corr == min(delta_corr(:));
[inda_1, inda_2] = ind2sub([N_a1 N_a2], find(Best_index == 1));

opt_a1 = alpha_1(inda_1)
opt_a2 = alpha_2(inda_2)
delta_corr
