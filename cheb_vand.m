function V = cheb_vand(n, x)
%cheb_vand calcola la matrioce di Vandermonde V della base dei punti x.
%   
%   INPUT:
%       n := [1 x 1] grado massimo di approssimazione del polinomio;
%       x := [n+1 x 1] oppure [1 x n+1] vettore dei punti appartenente all'intervallo [-1,1].
%
%   OUTPUT
%       V := [length(x) x n+1] matrice di Vandermonde.

x=x(:);

V=cos(acos(x)*(0:n));

end