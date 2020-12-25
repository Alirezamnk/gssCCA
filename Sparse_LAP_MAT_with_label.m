function [b] = Sparse_LAP_MAT_with_label(X_Dim, Y_Dim, LABEL)

% clear %all;
% close all;
% clc;
% 
% global X_Dim
% global Y_Dim
% 
% X_Dim = 6;
% Y_Dim = 6;

%% ---------------------------------------------------------------------------------------------------------------------------
%                                                             Modality #1     Group #1                                                   %
%%----------------------------------------------------------------------------------------------------------------------------

% LABEL = zeros(32,32);
% dev = 1;
% LABEL(1:8 , 1:8)           = 0 * dev * ones(8,8);
% LABEL(1:8 , 9:16)         = 11 * dev * ones(8,8);
% LABEL(1:8 , 17:24)       = 9 * dev * ones(8,8);
% LABEL(1:8 , 25:32)       = 1 * dev * ones(8,8);
% 
% LABEL(9:16 , 1:8)        = 6 * dev * ones(8,8);
% LABEL(9:16 , 9:16)      = 15 * dev * ones(8,8);
% LABEL(9:16 , 17:24)    = 2 * dev * ones(8,8);
% LABEL(9:16 , 25:32)    = 13 * dev * ones(8,8);
% 
% LABEL(17:24 , 1:8)     = 7 * dev * ones(8,8);
% LABEL(17:24 , 9:16)   = 3 * dev * ones(8,8);
% LABEL(17:24 , 17:24) = 10 * dev * ones(8,8);
% LABEL(17:24 , 25:32) = 12 * dev * ones(8,8);
% 
% LABEL(25:32 , 1:8)     = 4 * dev * ones(8,8);
% LABEL(25:32 , 9:16)   = 8 * dev * ones(8,8);
% LABEL(25:32 , 17:24) = 14 * dev * ones(8,8);
% LABEL(25:32 , 25:32) = 5 * dev * ones(8,8);

temp = -1 * ones(X_Dim + 2 , Y_Dim + 2);
temp(2 : X_Dim + 1 , 2 : Y_Dim + 1) = LABEL;

% figure, imshow(temp, [], 'InitialMagnification','fit');
% 
% P = X_Dim * Y_Dim;

U = [];
V = [];
W = [];

for i = 2 : X_Dim + 1
    for j = 2 : Y_Dim + 1
        [u, v, w] = Find_neighberhood(i, j, temp(i-1 : i+1 , j-1 : j+1), X_Dim);
        U = [U u];
        V = [V v];
        W = [W w];
    end
end

b = sparse(U , V , W);
b_hat = b - diag(diag(b));
b = b + b_hat';

% spy(b);
% c = full(b);
% figure,

% b_orig = b;
% [m, n] = size(b);
% [X, Y] = meshgrid(1:m, 1:n);
% b = (X + Y) .* b;
% 
% nonzeroInd = find(b);
% [x, y] = ind2sub([m n], nonzeroInd);
% 
% figure();
% hp = patch(x, y, b_orig(nonzeroInd), ...
%            'Marker', 's', 'MarkerFaceColor', 'flat', 'MarkerSize', 4, ...
%            'EdgeColor', 'none', 'FaceColor', 'none');
% set(gca, 'XLim', [0, n + 1], 'YLim', [0, m + 1], 'YDir', 'reverse', ...
%     'PlotBoxAspectRatio', [n + 1, m + 1, 1]);
% 
% colorbar();
% colormap(flipud(hot))