function [R , mu , C] = standardize(A , mu , C)

%--------------------------------------------
% To have zero mean and unit variance
%--------------------------------------------

[M, N] = size(A);
if (nargin == 1)
    sd = std(A, 0, 1);

    sd(sd == 0) = 1;
    if (min(sd) == 0) 
        error('Cannot standardize because some of the columns have std. dev. 0');
    end

% [M, N] = size(A);
% mu = mean(A,1);
% B = A - ones(M,1)*mu;

    mu = mean(A,1);
    B = A - ones(M,1)*mu;
    C = std(A, 0, 1); 
    C(C == 0) = 1;
    R = B(:,1:end) ./ (ones(M,1) * C);
else
    B = A - ones(M,1)*mu;
    R = B(:,1:end) ./ (ones(M,1) * C);
end
    
% R = B;

% [M, N] = size(A);
% mu = mean(A,1);
% B = A - ones(M,1)*mu;
% C = std(A, 0, 1); 
% R = B(:,1:end) ./ (ones(M,1) * C);

% mu = mean(A,2);
% 
% B = A - mu * ones(1,N);
% [M, N] = size(A);
% mu = mean(A,2);
% 
% B = A - mu * ones(1,N);
% C = std(A, 0, 2); 
% R = B(1:end,:) ./ (C * ones(1,N));% C = std(A, 0, 1);
 
% R = B(:,1:end) ./ (ones(M,1) * C);

