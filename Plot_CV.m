function Plot_CV(CV, CV_test , G_train , G_test , text_index , N_train_G1 , N_test_G1)

[M1, N] = size(CV);
M2 = size(CV_test, 1);

if N == 1
    CV         = [CV             zeros(M1,1)];
    CV_test = [CV_test     zeros(M2,1)];
end

figure, hold on;

plot(CV(1 : N_train_G1 , 1) , CV(1 : N_train_G1 ,  2),'LineWidth',3,'Color','blue','Marker','square', 'LineStyle','none');
plot(CV_test(1 : N_test_G1 , 1) , CV_test(1 : N_test_G1 , 2),'LineWidth',3,'Color','black','Marker','>', 'LineStyle','none');

plot(CV(N_train_G1 + 1 : end, 1) , CV(N_train_G1 + 1 : end, 2),'LineWidth',3,'Color','red','Marker','o', 'LineStyle','none');
plot(CV_test(N_test_G1 + 1 : end, 1),CV_test(N_test_G1 + 1 : end, 2),'LineWidth',3,'Color','yellow','Marker','<', 'LineStyle','none');

labels = num2str(G_train','%d');
text(CV(:,1) , CV(:,2) , labels, 'horizontal','left', 'vertical','bottom')

labels = num2str(G_test','%d');
text(CV_test(:,1) , CV_test(:,2) , labels, 'horizontal','left', 'vertical','bottom')

if N == 1
    set(gca,'YTickLabel',{},'box','off','fontsize',12,'fontweight','bold')
end
if N == 2
    set(gca,'fontsize',12,'fontweight','bold')
end

mylegend{1} = ['HC'];
mylegend{2} = ['test sample (HC)'];
mylegend{3} = ['SZ'];
mylegend{4} = ['test sample (AD)'];
legend(mylegend);

xlabel(strcat(text_index,'_',num2str(1)),'fontsize',12,'fontweight','bold')
ylabel(strcat(text_index,'_',num2str(2)),'fontsize',12,'fontweight','bold')

if N == 1
    title('1D classification','FontWeight','bold','FontSize',12);
else
    if N == 2
        title('2D classification','FontWeight','bold','FontSize',12);
    end
end
hold off;