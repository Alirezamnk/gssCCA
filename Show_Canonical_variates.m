function [] = Show_Canonical_variates(CV_1 , CV_2 , corr , i_K , N_train_G1 , N_train_G2)

Tot = N_train_G1 + N_train_G2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                   Curve Fitting                              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xData, yData] = prepareCurveData( 1 : Tot, CV_1 );
% 
% % Set up fittype and options.
% ft = fittype( 'poly5' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Normalize = 'on';
% opts.Robust = 'Bisquare';
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure, 
plot(CV_1(1:N_train_G1),'LineWidth',4,'Color','red');
hold on,
grid on,
plot(N_train_G1 : Tot , CV_1(N_train_G1 : end),'LineWidth',4,'Color','black');
plot(CV_2(1:N_train_G1),'LineWidth',4,'Color','blue');
plot(N_train_G1 : Tot , CV_2(N_train_G1 : end),'LineWidth',4,'Color','yellow');
xlabel('Subject #','fontsize',15,'fontweight','bold')
title(sprintf('Correlation:   %0.2f%%',corr*100),'FontWeight','bold','FontSize',15);


% plot( fitresult, xData, yData );
% plot( fitresult);
% mylegend{5} = ['Curve fitted to 1st CV1, Mod #1'];

mylegend{1} = [strcat('[HC]-' , 'CV Mod #1-' , '(U' , num2str(i_K) , ')')];
mylegend{2} = [strcat('[SZ]-' , 'CV Mod #1-' , '(U' , num2str(i_K) , ')')];
mylegend{3} = [strcat('[HC]-' , 'CV Mod #2-' , '(V' , num2str(i_K) , ')')];
mylegend{4} = [strcat('[SZ]-' , 'CV Mod #2-' , '(V' , num2str(i_K) , ')')];

% mylegend{1} = ['CV Modality #1'];
% mylegend{2} = ['CV Modality #2'];
legend(mylegend,'Location','SouthWest');
% legend(mylegend);
set([gca gca],'box','off','fontsize',15,'fontweight','bold')
hold off;

print(h , 'Canonical Variates' , '-dtiff' , '-r150')

print(h , '-depsc' , 'Canonical Variates.eps')
% print –dmeta;
% print (h , '-dmeta' , '-painters' , '-r864' , 'dataFigure');
