function [Us , indiu] = Bilinear_Convergence_Func(X , Y , Us , v , S , L , lam , pos)

oldui = Us;
utild = Us + (X'*(Y*v) - S*Us)/L;
Us = soft_thr(utild,lam/L,pos);

[Us] = Laplacian_Normalization(Us , S);

indiu = norm(Us - oldui) / norm(oldui);
