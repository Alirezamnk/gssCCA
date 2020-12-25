function [] = Permutation_test(X , Y , u , v)

M = size(X,1);
Perm_NUM = 10000;
PDF = zeros(Perm_NUM,1);

for i = 1 : Perm_NUM
    order =  randperm(M);
    
    P_X = X(order , :);

    PDF(i) = abs(corr(P_X * u , Y * v));
end
figure,
hist(PDF , 100); 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w')

cor1 = corr(X * u , Y * v);
title(sprintf('Correlation:   %0.2f%%',cor1*100),'FontWeight','bold','FontSize',12);
xlabel('Inertia of the Permuted Sample','fontsize',12,'fontweight','bold')
ylabel('Number of Samples (out of 10,000)','fontsize',12,'fontweight','bold')

% disp('Click graph to place arrow; first tail, then head:')
% [axx axy] = ginput(2);   % Returns list of x, list of y in data space
% % Transform from data space to figure space
% [arrowx,arrowy] = dsxy2figxy(gca, axx, axy);
% har = annotation('textarrow',arrowx,arrowy);
% content = sprintf('(%4.2f,%4.2f)',axx(2), axy(2));
% % Plot anno text centered at the tail of the arrow
% set(har,'String',content,'Fontsize',8)
% 
% xlabel('Integer value');
% ylabel('Probability');
