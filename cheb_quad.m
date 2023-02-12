function [xquad,w]=cheb_quad(n)
%cheb_quad calcola il vettore dei nodi di quadratura ed il vettore dei
%relativi pesi risolvendo il problema dei momenti Vquad'w = m. Il sistema 
%di equazioni generato dal problema dei momenti viene risolto mediante il
%metodo di fattorizzazione LU con pivoting parziale e gli algoritmi di
%sostituzione (sostituzioneAvanti, sostituzioneIndietro)
%
%   INPUT:
%       n := [1 x 1] grado di esattezza polinomiale.
%
%   OUTPUT:
%       xquad := [n+1 x 1] vettore dei nodi di quadratura
%       w := [1 x n+1] vettore dei pesi
%

xquad=cos((0:n)./n*pi);
Vquad = cheb_vand(n,xquad);

m=[2;zeros(n,1)];
for k = 2:2:n+1
    m(k) = 2/(1-k^2);
end

[L,U,P] = lu(Vquad');
w=SostituzioneIndietro(U,SostituzioneAvanti(L,P*m))';

end 