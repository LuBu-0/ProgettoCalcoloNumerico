function [cstar,Rzero,b]=cheb_tikonov(n,lambda,xsample,ysample)
%cheb_tikonov Costruisce la Matrice Rzero che corrisponde alla fattorizzazione
%di Cholesky della matrice G := R0'R0 dove R0 è la parte quadrata
%superiore di R a sua volta ottenuta dalla fattorizzazione QR della matrice
%L t.c. G = L'L ed il termine noto b del sistema lineare Gc = b in modo da
%risolvere poi tale sistema e fornire il valore cstar soluzione dello stesso.
%   
%   INPUT:
%       n := [1 x 1] grado polinomiale;
%       lambda := [1 x 1] parametro di regolarizzazione alla Tikonov;
%       xsample := [M x 1] valori del data-set;
%       ysample := [M x 1] ysample = f(xsample) + rumore.
%
%   OUTPUT
%       cstar := [n+1 x 1] vettore minimizzatore della funzione convessa;
%       Rzero := [n+1 x n+1] fattorizzazione di Choleski di G;
%       b := [n+1 x 1] termine noto del sistema lineare Gc=b.
xsample=xsample(:);
ysample=ysample(:);
M = length(xsample);
%costruzione di R0
V=cheb_vand(n,xsample);
[xquad,w] = cheb_quad(n);
Vquad=cheb_vand(n,xquad);
L=[V/sqrt(M);diag(sqrt(lambda*w))*Vquad];
[Q,R] = qr(L);
R0 = R(1:n+1,:);
%costruzione del termine noto b
b=(V'*ysample) / M;
%calcolo della soluzione cstar tramite algoritmi di sostituzione
cstar=SostituzioneIndietro(R0,SostituzioneAvanti(R0',b))';
Rzero=R0'*R0;
end