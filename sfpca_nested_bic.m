function [U, V, d, optaus, optavs, optlus, optlvs, bicu, bicv , Lu , Lv] = sfpca_nested_bic(X, Y, K, lamus, ...
                       lamvs, alphaus, alphavs, Omegu, Omegv, startu, startv, posu, posv, maxit, iterS , Lu , Lv, Permut_test)
                     
% [p, q] = size(X' * Y);
p = size(X,2);
q = size(Y,2);

ru = length(alphaus); 
rv = length(alphavs); 
rlu = length(lamus); 
rlv = length(lamvs); 

tic

% Lu = max(eig(speye(p) + p * max(alphaus) * Omegu)) + .01;%Orig
% Lv = max(eig(speye(q) + q * max(alphavs) * Omegv)) + .01;%Orig
if (Lu == 0 && Lv == 0)
    Lu = eigs(speye(p) + p * max(alphaus) * Omegu , 1) + .01;%    
    if (p == q && alphaus == alphavs)
        Lv = Lu;
    else
        Lv = eigs(speye(q) + q * max(alphavs) * Omegv , 1) + .01;%
    end
elseif (Lu == 0)
    if (p == q && alphaus == alphavs)
        Lu = Lv;
    else
        Lu = eigs(speye(p) + p * max(alphaus) * Omegu , 1) + .01;%    
    end
elseif (Lv == 0)
    if (p == q && alphaus == alphavs)
        Lv = Lu;
    else
        Lv = eigs(speye(q) + q * max(alphavs) * Omegv , 1) + .01;%
    end
end
% Lu = max(eig(speye(p) + max(alphaus) * Omegu)) + .01;
% Lv = max(eig(speye(q) + max(alphavs) * Omegv)) + .01;
% Lu = max(eig(X'*X + p * max(alphaus) * Omegu)) + .01;    % for CCA test
% Lv = max(eig(Y'*Y + q * max(alphavs) * Omegv)) + .01;     % for CCA test

lamus = sort(lamus,'descend');
lamvs = sort(lamvs,'descend');

optaus = []; 
optavs = []; 
optlus = []; 
optlvs = [];

% Zhat = X' * Y; 
% toc
% tic

U = zeros(p,K);
V = zeros(q,K);

for k=1:K
    bicu = zeros(ru,rlu); 
    bicv = zeros(rv,rlv);
    
    True = 1;
    while (True == 1)
    if (sum(startu) == 0)
%         [u,~,v] = svds(Zhat,1);
        [u, v] = svd_using_QR(X, Y, 1);
%         u = mean(X,1)';       v = mean(Y,1)';
%         u = median(X,1)';     v = median(Y,1)';
%         u = max(X, [], 1)';       v = max(Y, [], 1)';
%         u = min(X, [], 1)';       v = min(Y, [], 1)';
%         u = mode(X,1)';       v = mode(Y,1)';
%         u = randn(n,1);         v = randn(p,1);
        u = u/norm(u);%Orig
        v = v/norm(v);%Orig
    else
        u = startu(:,k); 
        v = startv(:,k);    
    end
    
    if abs(max(u)) < abs(min(u))
        True = 0;
    end
    end
    
    Data_DIR = '/home/mszam12/main/matlab/MATLAB_Codes/Make_Real_Data/Multi-modal_Data/IKHH/';
    Tagged_X_DIR = strcat(Data_DIR , 'Tagged_Obj_X_T1_IKHH_49_950629');
    Show_Canonical_Coeffitient_X(u , strcat('/W1_ALL_1.nii') , Tagged_X_DIR);

    indo = 1;
    iter = 1;
    thr = 1e-6;
    
    count = 0;% a counter for divergence detection
    count_max = 10;
    min_indo = Inf;
    
    while ( (indo > thr) && (iter < min([iterS maxit])) )
        oldu = u; 
        oldv = v;

        us = zeros(p,ru,rlu);
        for j = 1:ru
            Su = speye(p) + alphaus(j)*p*Omegu;     %Orig
            for i = 1:rlu
                indiu = 1;                
                count_indu = 0;
                min_indu = Inf;
                while indiu > thr                  
                    [us(:,j,i) , indiu] = Bilinear_Convergence_Func(X , Y , us(:,j,i) , v , Su , Lu , lamus(i) , posu);
                    [min_indu , count_indu , exit] = Divergence_Detector(indiu , min_indu , count_indu , count_max);                   
                    if exit == 1
                        break;
                    end
                end
                
%                 actu = us(:,j,i) ~= 0;
%                 dfu = trace(inv( speye(sum(actu)) + n*alphaus(j)*Omegu(actu,actu)));
%                 bicu(j,i) = log( norm(X'*(Y*v) - utild)^2/n ) + .5*log(n)*dfu/n;
            end
        end
%         iu = bicu==min(min(bicu));
%         [indau,indlu] = ind2sub([ru rlu],find(iu==1));
% %         optlu = lamus(indlu); 
% %         optau = alphaus(indau);
% %         u = us(:,indau,indlu);
%         optlu = lamus(min(indlu)); 
%         optau = alphaus(min(indau));
%         u = us(:,min(indau),min(indlu));
        optlu = lamus; 
        optau = alphaus;
        u = us;
    
        vs = zeros(q,rv,rlv);
        for j=1:rv
            Sv = speye(q) + alphavs(j)*q*Omegv;   %Orig
            for i=1:rlv
                indiv = 1;
                count_indv = 0;
                min_indv = Inf;
                while indiv>thr
                    [vs(:,j,i) , indiv] = Bilinear_Convergence_Func(Y , X , vs(:,j,i) , u , Sv , Lv , lamvs(i) , posv);
                    [min_indv , count_indv , exit] = Divergence_Detector(indiv , min_indv , count_indv , count_max);                   
                    if exit == 1
                        break;    
                    end
                end
%                 actv = vs(:,j,i) ~= 0;
%                 dfv = trace(inv( speye(sum(actv)) + p*alphavs(j)*Omegv(actv,actv)));
%                 bicv(j,i) = log( norm(Y'*(X*u) - vtild)^2/p ) + .5*log(p)*dfv/p;
            end
        end
%         iv = bicv==min(min(bicv));
%         [indav,indlv] = ind2sub([rv rlv],find(iv==1));
%         optlv = lamvs(min(indlv)); 
%         optav = alphavs(min(indav));
%         v = vs(:,min(indav),min(indlv));
        optlv = lamvs; 
        optav = alphavs;
        v = vs;

        iter = iter + 1;
        indo = norm(u - oldu)/norm(oldu) + norm(v- oldv)/norm(oldv);

        [min_indo , count , exit] = Divergence_Detector(indo , min_indo , count , count_max);                   
        if exit == 1
            break;
        end
    end
    
    clear vs; 
    clear us;
    
    Su = speye(p) + p*optau*Omegu;%Orig
    Sv = speye(q) + q*optav*Omegv;%Orig

    while indo>thr && iter<maxit
        oldu = u; 
        oldv = v;

        indiu = 1;        
        count_indu = 0;
        min_indu = Inf;
        while indiu > thr
            [u , indiu] = Bilinear_Convergence_Func(X , Y , u , v , Su , Lu , optlu , posu);
            [min_indu , count_indu , exit] = Divergence_Detector(indiu , min_indu , count_indu , count_max);                   
            if exit == 1
                break;
            end
        end
        
        indiv = 1;
        count_indv = 0;
        min_indv = Inf;
        while indiv > thr
            [v , indiv] = Bilinear_Convergence_Func(Y , X , v , u , Sv , Lv , optlv , posv);
            [min_indv , count_indv , exit] = Divergence_Detector(indiv , min_indv , count_indv , count_max);                   
            if exit == 1
                break;
            end            
        end

        indo = norm(oldu - u)/norm(oldu) + norm(oldv - v)/norm(oldv);
        iter = iter + 1;

        [min_indo , count , exit] = Divergence_Detector(indo , min_indo , count , count_max);                   
        if exit == 1
            break;
        end
    end
        
    u = u/norm(u);%Orig
    v = v/norm(v);%Orig

    if (Permut_test == 1)
        Permutation_test(X , Y , u , v);
    end

    tx = X * u;%test
    ty = Y * v;%test

    X = X - (tx * tx' * X) / (tx' * tx);%test
    Y = Y - (ty * ty' * Y) / (ty' * ty);%test
%     Zhat = X' * Y;%test

    dd = (u'*X')*(Y*v);%Orig
%     Zhat = Zhat - dd*u*v';%Orig
    
    U(:,k) = u; %Orig
    V(:,k) = v; %Orig
%     U(:,k) = test2; %test
%     V(:,k) = test4; %test
    d(k) = dd;
    
    toc

    optlus(k) = optlu;    
    optlvs(k) = optlv; 
    optaus(k) = optau; 
    optavs(k) = optav; 
end