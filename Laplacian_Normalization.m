function [u] = Laplacian_Normalization(u , Su)

n = size(u,1);

if norm(u) > 0
    u  = u / sqrt(u' * Su * u);
else
    u = zeros(n,1);
end            
