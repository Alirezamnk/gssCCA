function [] = Optimal_alpha_selection(X , Y , alphaus , alphavs , Omegu , Omegv , posu , posv , iterS)

global n_CV

n = size(X,2);
p = size(Y,2);

lamus = 0;
lamvs = 0;

N_alpha_u = length(alphaus); 
N_alpha_v = length(alphavs); 

Lu = eigs(speye(n) + n * max(alphaus) * Omegu , 1) + .01;%
Lv = eigs(speye(p) + p * max(alphavs) * Omegv , 1) + .01;%
% Lu = .01;%
% Lv = .01;%

% optaus = []; 
% optavs = []; 
% optlus = []; 
% optlvs = [];

% bicu = zeros(N_alpha_u , rlu); 
% bicv = zeros(N_alpha_v , rlv);

u = mean(X , 1)';       
v = mean(Y , 1)';
%         u = randn(n,1);         v = randn(p,1);
u = u/norm(u);%Orig
v = v/norm(v);%Orig

thr = 1e-6;
count_max = 5;
Delta_Corr_1 = zeros(N_alpha_u , N_alpha_v);
Delta_Corr_2 = zeros(N_alpha_u , N_alpha_v);

for i = 1 : N_alpha_u
    for j = 1 : N_alpha_v
        for i_CV = 1 : n_CV
            [X_train , X_test , Y_train , Y_test] = Do_Cross_Validation(X , Y , i_CV);
        
            [X_train , mu_X , C_X] = standardize(X_train);
            [Y_train , mu_Y , C_Y] = standardize(Y_train);

            count = 0;
            indo = 1; 
            min_indo = Inf;
            iter = 1; 
            Su = speye(n) + alphaus(i)*n*Omegu;     %Orig
            Sv = speye(p) + alphavs(j)*p*Omegv;     %Orig

            while ( (indo > thr) && (iter < iterS) )
                oldu = u; 
                oldv = v;

                indiu = 1;                
                count_indu = 0;
                min_indu = Inf;
                us = zeros(n,1);
                while indiu > thr                  
                    [us , indiu] = Bilinear_Convergence_Func(X_train , Y_train , us , v , Su , Lu , 0 , posu);
                    [min_indu , count_indu , exit] = Divergence_Detector(indiu , min_indu , count_indu , count_max);
                    if exit == 1
                        break;
                    end
                end
                u = us;

                indiv = 1;
                count_indv = 0;
                min_indv = Inf;
                vs = zeros(p,1);
                while indiv>thr
                    [vs , indiv] = Bilinear_Convergence_Func(Y_train , X_train , vs , u , Sv , Lv , 0 , posv);
                    [min_indv , count_indv , exit] = Divergence_Detector(indiv , min_indv , count_indv , count_max);                   
                    if exit == 1
                        break;
                    end
                end
                v = vs;

                iter = iter + 1;
                indo = norm(u - oldu)/norm(oldu) + norm(v- oldv)/norm(oldv);

                [min_indo , count , exit] = Divergence_Detector(indo , min_indo , count , count_max);                   
                if exit == 1
                    break;
                end
            end

            W1 = u/norm(u);%Orig
            W2 = v/norm(v);%Orig

            X_test = standardize(X_test , mu_X , C_X);
            Y_test = standardize(Y_test , mu_Y , C_Y);

            U_train = X_train * W1;
            V_train = Y_train * W2;

            U_test = X_test * W1;
            V_test = Y_test * W2;
        
            Delta_Corr_1(i , j) = Delta_Corr_1(i , j) + abs(corr(U_test , V_test));
            Delta_Corr_2(i , j) = Delta_Corr_2(i , j) + abs(corr(U_train , V_train) - corr(U_test , V_test));
        end        
        Delta_Corr_1(i , j) = Delta_Corr_1(i , j) / n_CV;
        Delta_Corr_2(i , j) = Delta_Corr_2(i , j) / n_CV;
    end
end

return

%====================================================================

% while ( (indo > thr) && (iter <iterS) )
%     oldu = u; 
%     oldv = v;
% 
%     us = zeros(n,N_alpha_u,rlu);
%     for j = 1:N_alpha_u
%         Su = speye(n) + alphaus(j)*n*Omegu;     %Orig
%         for i = 1:rlu
%             indiu = 1;                
%             count_indu = 0;
%             min_indu = Inf;
%             while indiu > thr                  
%                 [us(:,j,i) , indiu] = Bilinear_Convergence_Func(X , Y , us(:,j,i) , v , Su , Lu , 0 , posu);
%                 [min_indu , count_indu] = Divergence_Detector(indiu , min_indu , count_indu , count_max);                   
%             end
%         end
%     end
%     optlu = lamus; 
%     optau = alphaus;
%     u = us;
% 
%     vs = zeros(p,N_alpha_v,rlv);
%     for j=1:N_alpha_v
%         Sv = speye(p) + alphavs(j)*p*Omegv;   %Orig
%         for i=1:rlv
%             indiv = 1;
%             count_indv = 0;
%             min_indv = Inf;
%             while indiv>thr
%                 [vs(:,j,i) , indiv] = Bilinear_Convergence_Func(Y , X , vs(:,j,i) , u , Sv , Lv , 0 , posv);
%                 [min_indv , count_indv] = Divergence_Detector(indiv , min_indv , count_indv , count_max);                   
%             end
%         end
%     end
%     optlv = lamvs; 
%     optav = alphavs;
%     v = vs;
% 
%     iter = iter + 1;
%     indo = norm(u - oldu)/norm(oldu) + norm(v- oldv)/norm(oldv);
% 
%     [min_indo , count] = Divergence_Detector(indo , min_indo , count , count_max);                   
% end
