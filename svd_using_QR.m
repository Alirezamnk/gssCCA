function [startu, startv] = svd_using_QR(X, Y, num)

[Q , R] = qr(X' , 0);
[P , T] = qr(Y' , 0);

[u , ~ , v] = svds(R * T' , num);

startu = Q * u;
startv = P * v;
