function [b] = Sparse_LAP_MAT(X_Dim, Y_Dim, L1, L2, L3, connectivity)

% clear %all;
% close all;
% clc;

% X_Dim = 4;
% Y_Dim = 4;

P = X_Dim * Y_Dim;

%% Diagonal
first = [L1  L2*ones(1,X_Dim-2)   L1];
other =  [L2  L3*ones(1,X_Dim-2)   L2];

D = zeros(1,P);

D(1,1 : X_Dim) = first;
for i=2:Y_Dim-1
    D(1,(i-1)*X_Dim + 1 : i*X_Dim) = other;
end 
D(1,i*X_Dim + 1 : (i+1)*X_Dim) = first;

u  = [1:P];
v  = [1:P];
w = [D];
clear other D  P  first;

%%  Right & Left
first = [];
second = [];

for i=1:X_Dim
    first       = [first        (i-1)*Y_Dim + 1   :   i * Y_Dim - 1];
    second = [second  (i-1)*Y_Dim + 2   :   i * Y_Dim];
end

last = [-1*ones(1,X_Dim)];
R = repmat(last,1,Y_Dim - 1);

u  = [u  first         second];
v  = [v  second   first];
w = [w     R            R];
clear first second R  last;

if (connectivity == 2)
    b = sparse(u , v , w);
    return;
end

%% UP & Down
first_U       = [1               : (X_Dim - 1) * Y_Dim];
second_U = [Y_Dim + 1 : X_Dim * Y_Dim];
last_U = [-1 * ones(1,(X_Dim - 1) * Y_Dim)];

u = [u    first_U         second_U];
v = [v    second_U   first_U];
w = [w   last_U         last_U];

clear first_U second_U last_U;

%% Upper and lower right
first_UR = [];
second_UR = [];

for i=1:X_Dim - 1
    first_UR       = [first_UR         (i - 1) * Y_Dim + 1 + 1  :   i * Y_Dim];
    second_UR = [second_UR   i * Y_Dim + 1                :  (i + 1)  * Y_Dim - 1];
end

last = [-1*ones(1,Y_Dim - 1)];
UR = repmat(last,1,X_Dim - 1);

u  = [u    first_UR          second_UR];
v  = [v    second_UR    first_UR];
w = [w   UR                  UR];

clear first_UR     second_UR     UR   last;

%% Upper and lower left
first_LR = [];
second_LR = [];

for i=1:X_Dim - 1
    first_LR       = [first_LR         (i - 1) * Y_Dim + 1  :   i  * Y_Dim - 1];
    second_LR = [second_LR   i * Y_Dim + 1 + 1    :   (i + 1) * Y_Dim];
end

last = [-1*ones(1,Y_Dim - 1)];
UL = repmat(last,1,X_Dim - 1);

%%
u  = [u   first_LR         second_LR];
v  = [v   second_LR   first_LR];
w = [w  UL                 UL];

clear first_LR     second_LR    last   UL;

b = sparse(u , v , w);
% spy(b)