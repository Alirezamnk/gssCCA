clear;
close all;
clc;

global n_CV
global whole_sample

Excel_monitoring = 0;                      %1: Monitor    0: Non-Monitor
Read_saved_monitoring = 1;                 %1: Read the previously saved indices    0: Re-arrenge indices
find_best_Groupness_parameters = 0;        %1: Find    0: No
find_best_sparsity_parameters = 0;         %1: Find    0: No
Do_Cross_Validation = 0;
Do_Permutation_test = 0;
Find_Result_for_All = 1;
NUM = 1;

%-------------------------------------------------------------------------------------
%                                       Excel setup                                          %
%-------------------------------------------------------------------------------------

fname = 'Book1.xlsx';
sname_1 = 'Sheet1';
sname_2 = 'Sheet2';
[~ , ~ , Data] = xlsread(fname , sname_1);      %read in the old data, text and all
nextRow = size(Data , 1) + 1;                 %get the row number of the end

range_1 = sprintf('%s%d : %s%d' , 'A' , nextRow , 'AZ' , nextRow);  %this tells excel where to stick it
range_2 = sprintf('%s%d : %s%d' , 'A' , nextRow , 'AG' , nextRow);  %this tells excel where to stick it
range_3 = sprintf('%s%d : %s%d' , 'A' , nextRow + 1 , 'AZ' , nextRow + 1);  %this tells excel where to stick it
range_4 = sprintf('%s%d : %s%d' , 'A' , nextRow + 1 , 'AG' , nextRow + 1);  %this tells excel where to stick it

%-------------------------------------------------------------------------------------
Data_DIR = 'Path to the pre-processed data';

% Read the input X and Y matrices
X = importdata(strcat(Data_DIR , 'Modality_X.mat'));
Y = importdata(strcat(Data_DIR , 'Modality_Y.mat'));

% When we convert a 3D volume image to a vector, we need to save the co-ordinates of each voxel
Tagged_X_DIR = strcat(Data_DIR , 'Tagged_Modality_X');
Tagged_Y_DIR = strcat(Data_DIR , 'Tagged_Modality_Y');

[n, p] = size(X);
[~, q] = size(Y);

% List of Subject IDs (We can prepare a text file and then read the list based on that txt file)
G_1 = [002, 010, 013, 014, 017, 018, 021, 022, 025, 026, 029, 030, 033, 037, 045, 049, 053, 057, 061, 003, ...
            019, 020, 059, 071, 031];
G_2 = [004, 007, 008, 011, 012, 015, 016, 023, 024, 027, 028, 039, 043, 047, 051, 055, 063, 067, 075, 001, ...
            005, 006, 009, 041];

G = [G_1   G_2];
N_G1 = size(G_1,2);
N_G2 = size(G_2,2);

% Read the Laplacian matrix of each modality
Laplac_DIR = 'Path to the Laplacian Matrices';
Omegu = importdata(strcat(Laplac_DIR , 'LX.mat'));
Omegv = importdata(strcat(Laplac_DIR , 'LY.mat'));

%-------------------------------------------------------------------------------------
%                               Parameter Setting                                        %
%-------------------------------------------------------------------------------------

% Rang of the hyper-parameters
opt_lam_1 = (1);
opt_lam_2 = (2);

opt_alpha_1 = (0.1);
opt_alpha_2 = (0.05);
%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------------------
%                Produce K-fold labels for all of the subjects                %
%-------------------------------------------------------------------------------------

n_CV = 5;                                        % Number of Cross-Validation
whole_sample = 1 : n;
posu = 1;       posv = 1;       maxit = 1e3;        iterS = 20;

for num = 1:NUM
    if (Read_saved_monitoring == 0)
        species_type = cell(size(G,2),1);
        labels = {'HC','SZ'};
        species_type(1 : size(G_1,2)) = labels(1);
        species_type(size(G_1,2) + 1 : end) = labels(2);

        indices = crossvalind('Kfold' , species_type , n_CV);
        save('indices.mat' , 'indices');
    else
        indices = importdata(strcat(pwd,'/indices.mat'));
    end

%------------------------------------------------------------------------------------%
%                               Do a sub-optimal K-fold CV to find the best Groupness parameters                          %
%------------------------------------------------------------------------------------%

    if (find_best_Groupness_parameters == 1)

        lam_1 = (0);
        lam_2 = (0);

        alpha_1 = (0.20 : 0.004 : 0.30);
        alpha_2 = (0.20 : 0.004 : 0.30);

        [opt_alpha_1, opt_alpha_2] = Find_sub_opt_alphas(X, Y, indices, lam_1, lam_2, alpha_1, alpha_2, Omegu, Omegv, ...
                                                                                            posu, posv, maxit, iterS);
    end

%------------------------------------------------------------------------------------%
%                           Do a sub-optimal K-fold CV to find the best Sparsity parameters                                   %
%------------------------------------------------------------------------------------%

    if (find_best_sparsity_parameters == 1)

        opt_alpha_1 = 0.3;
        opt_alpha_2 = 0.3;

        lam_1 = (1.0 : 0.1 : 2.0);
        lam_2 = (1.0 : 0.1 : 2.0);

        [opt_lam_1, opt_lam_2] = Find_sub_opt_lambdas(X, Y, indices, lam_1, lam_2, opt_alpha_1, opt_alpha_2, ...
                                                                                         Omegu, Omegv, posu, posv, maxit, iterS);
    end

%--------------------------------------------------------------------------------------------------------------

    if (Do_Cross_Validation == 1)
        Lu = 0;
        Lv = 0;
        for i_CV = 1 : n_CV
            temp_sample = whole_sample;
            testing_sample = whole_sample(indices == i_CV);
            temp_sample(testing_sample) = [];
            training_sample = temp_sample;

            N_train_G1 = sum(training_sample <= N_G1);
            N_train_G2 = sum(training_sample > (n - N_G2));

            N_test_G1 = N_G1 - N_train_G1;
            N_test_G2 = N_G2 - N_train_G2;

            X_train = X(training_sample, :);
            X_test = X(testing_sample, :);

            Y_train = Y(training_sample, :);
            Y_test = Y(testing_sample, :);

            [X_train , mu_X , C_X] = standardize(X_train);
            [Y_train , mu_Y , C_Y] = standardize(Y_train);

        %     Optimal_alpha_selection(X , Y , alphaus , alphavs , Omegu , Omegv , 0 , 0 , 10);

            K = 1;%Number of Canonical Variates
            U = zeros(N_train_G1 + N_train_G2 , K);
            V = zeros(N_train_G1 + N_train_G2 , K);
            [W1, W2, ~, ~, ~, ~, ~, ~, ~ , Lu , Lv] = sfpca_nested_bic(X_train, Y_train, K, opt_lam_1, opt_lam_2, ...
                    opt_alpha_1, opt_alpha_2, Omegu, Omegv, 0, 0, posu, posv, maxit, iterS , Lu , Lv, Do_Permutation_test);

            for i_K = 1 : K
                Show_Canonical_Coeffitient(W1(:,i_K) , strcat('\W1_CV_',num2str(i_CV),'.nii'));
                Show_Canonical_Coeffitient(W2(:,i_K) , strcat('\W2_CV_',num2str(i_CV),'.nii'));

                U(: , i_K) = X_train * W1(:,i_K);
                V(: , i_K) = Y_train * W2(:,i_K);
                cor1 = corr(U(: , i_K) , V(: , i_K));

                Show_Canonical_variates(U(:,i_K) , V(:,i_K) , cor1(1,1) , i_K , N_train_G1 , N_train_G2);
            end

            if Excel_monitoring == 1
                Write_CV_data_2_Excel(U , V , training_sample , N_G1 , N_G2 , N_train_G1 , i_CV , fname , sname_1 , sname_2 , ...
                        range_1 , range_2 , range_3 , range_4);
            end

            X_test = standardize(X_test , mu_X , C_X);
            Y_test = standardize(Y_test , mu_Y , C_Y);

            U_test = X_test * W1;
            V_test = Y_test * W2;
            cor1 = corr(U_test(: , i_K) , V_test(: , i_K))

            training_label_vector = [ones(N_train_G1,1)    ;    -ones(N_train_G2,1)];

            hold off;
            figure
            model = svmtrain(U , training_label_vector,'ShowPlot',true , 'kernel_function' , 'rbf' , 'rbf_sigma' , 1 , 'method' , 'QP');
            species = svmclassify(model , U_test , 'ShowPlot' , true);
            [sum(species(1 : N_test_G1) == 1)  sum(species(1 : N_test_G1) == -1)     sum(species(N_test_G1 + 1 : end) == -1)     sum(species(N_test_G1 + 1 : end) == 1)]
            hold off;

            figure
            model = svmtrain(V , training_label_vector,'ShowPlot',true , 'kernel_function' , 'rbf' , 'rbf_sigma' , 1 , 'method' , 'QP');
            species = svmclassify(model , V_test , 'ShowPlot' , true);
            [sum(species(1 : N_test_G1) == 1)  sum(species(1 : N_test_G1) == -1)     sum(species(N_test_G1 + 1 : end) == -1)     sum(species(N_test_G1 + 1 : end) == 1)]

            Plot_CV(U , U_test , G(training_sample) , G(testing_sample) , 'U' , N_train_G1 , N_test_G1);
            Plot_CV(V , V_test , G(training_sample) , G(testing_sample) , 'V' , N_train_G1 , N_test_G1);

            close all;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                       Result for All data                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (Find_Result_for_All == 1)
        Lu = 0;
        Lv = 0;

        [X] = standardize(X);
        [Y] = standardize(Y);

        K = 1;%Number of Canonical Variates
        U = zeros(n , K);
        V = zeros(n , K);

        [W1, W2, dh, optaus, optavs, optlus, optlvs, bicu, bicv , Lu , Lv] = sfpca_nested_bic(X, Y, K, opt_lam_1, opt_lam_2, ...
                              opt_alpha_1, opt_alpha_2, Omegu, Omegv, 0, 0, posu, posv, maxit, iterS , Lu , Lv, Do_Permutation_test);

        for i_K = 1 : K
            Show_Canonical_Coeffitient_X(W1(:,i_K) , strcat('/W1_ALL_',num2str(num),'.nii') , Tagged_X_DIR);
            Show_Canonical_Coeffitient_Y(W2(:,i_K) , strcat('/W2_ALL_',num2str(num),'.nii') , Tagged_Y_DIR);

            U(: , i_K) = X * W1(:,i_K);
            V(: , i_K) = Y * W2(:,i_K);
            cor1 = corr(U(: , i_K) , V(: , i_K));

            Show_Canonical_variates(U(:,i_K) , V(:,i_K) , cor1(1,1) , i_K , N_G1 , N_G2);
            [opt_lam_1 , opt_lam_2 ; opt_alpha_1 , opt_alpha_2 ; mean(U(1:N_G1,1)) - mean(U(N_G1 + 1 : end,1)) , mean(V(1:N_G1,1)) - mean(V(N_G1 + 1 : end,1))]

            test_permutation(U(:,i_K));
            test_permutation(V(:,i_K));
        end
    end
end
