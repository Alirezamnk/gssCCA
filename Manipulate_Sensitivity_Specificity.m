function [SEN , SPE] = Manipulate_Sensitivity_Specificity(W1, Orig, sign)

[M, N] = size(Orig);

W1_map = sign * reshape(W1,M,N);
W1_map = (W1_map ~= 0);
Orig_map = (Orig ~= 0);
res = W1_map(:,:) & Orig_map(:,:);

TP = W1_map(:,:) & Orig_map(:,:);
TP = TP * 1;

FP = W1_map(:,:) - Orig_map(:,:);
FP = FP > 0;
FP = FP * 0.75;

FN = Orig_map(:,:) - W1_map(:,:);
FN = FN > 0;
FN = FN * 0.5;

S_W = sum(W1_map(:));
S_O = sum(Orig_map(:));
S_T = M * N;
S_R = sum(res(:));

SEN = S_R / S_O;
SPE = ((S_T - S_O) - (S_W - S_R)) / (S_T - S_O);

Image = TP + FP + FN;
figure, imshow(Image, [], 'InitialMagnification','fit','Colormap',jet(255));
title(sprintf('Sensitivity:  %0.2f%%  and Specificity:  %0.2f%%',SEN*100,SPE*100),'FontWeight','bold','FontSize',12);