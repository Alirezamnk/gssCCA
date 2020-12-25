function [u, v, w] = Find_neighberhood(i, j, sub_LABEL, X_Dim)

% global X_Dim

T = (sub_LABEL == sub_LABEL(2 , 2));
Total = sum(T(:)) - 1;

% T = triu(T) - diag(diag(T));
% Total_Neg = sum(T(:));

T(1,1) = 0;
T(1,2) = 0;
T(1,3) = 0;
T(2,1) = 0;

[o1, o2] = ind2sub([3 3],find(T==1));
Total_Forwad = sum(T(:));

x = [o1' + (i - 3)];
y = [o2' + (j - 3)];

for k = 1 : Total_Forwad
        u(k) = (i - 2) * X_Dim + (j - 1);
        v(k) = (x(k) - 1) * X_Dim + y(k);
end

v = sort(v);

w = [Total  -1*ones(1,Total_Forwad - 1)];
