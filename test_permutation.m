function [] = test_permutation(X)

M = 49;
Perm_NUM = 100000;
Perm_Mat = [];
N1 = 25;
N2 = 24;
PDF = zeros(Perm_NUM,1);

for i = 1 : Perm_NUM
    order =  randperm(M);

    Perm_Mat = [Perm_Mat ; order];
    
    [Perm_Mat] = unique(Perm_Mat , 'rows');
    
    if (size(Perm_Mat,1) ~= i)
        i = i - 1;
    else
        Perm_X = X(order);
        perm_X = [Perm_X(1:N1) ; Perm_X(N1 + (1 : N2)) * (N1 / N2)];
        
        mean_G1 = mean(Perm_X(1 : N1) , 1);
        mean_G2 = mean(Perm_X(N1 + (1 : N2)) , 1);

        tvar1 = var(Perm_X(1 : N1)) / N1;
        tvar2 = var(Perm_X(N1 + (1 : N2))) / N2;
        
        PDF(i) = abs(mean_G1 - mean_G2) ./ sqrt(tvar1 + tvar2);        
    end
end

X = [X(1:N1) ; X(N1 + (1 : N2)) * (N1 / N2)];
MU_1 = mean(X(1 : N1) , 1);
MU_2 = mean(X(N1 + (1 : N2)) , 1);

var1 = var(X(1 : N1)) / N1;
var2 = var(X(N1 + (1 : N2))) / N2;

Tvals = abs(MU_1 - MU_2) ./ sqrt(var1 + var2);
Tobs  = max(Tvals);

pval = mean(Tobs < PDF);
qval = quantile(PDF , 0.05);

figure,
hist(PDF , 100); 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w')

% cor1 = corr(X * u , Y * v);
title(sprintf('Inertia of the Original Sample:   %0.2f%',Tobs),'FontWeight','bold','FontSize',12);
xlabel('Inertia of the Permuted Samples','fontsize',12,'fontweight','bold')
ylabel('Number of Samples (out of 100,000)','fontsize',12,'fontweight','bold')
