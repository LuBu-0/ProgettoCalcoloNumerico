function D=ChebDiff1D(N)
%ChebDiff1D calcola la matrice di derivazione D della base di Chebishev di grado n.
%Questa matrice può essere utilizzata per calcolare la derivata di una funzione
%nella base di Chebyshev.
%
%   INPUT:
%       N := [1 x 1] intero positivo, rappresenta il grado della base di Chebyshev.
%
%   OUTPUT:
%       D:=[N+1 x N+1] matrice di derivazione della base di Chebishev.

i=1:N;
[I,J]=meshgrid(i,i);
D0=tril(-J.*((-1).^(I+J)-1));
j=J(:,1);D1=(1-(-1).^j)./2.*j;
D=[zeros(1,N+1);D1,D0]';