function x=SostituzioneIndietro(U,b)
%SostituzioneIndietro risolve il sistema lineare Ux = b mediante l'algoritmo
%di sostituzioine indietro
%   
%   INPUT:
%       U := [n x n] matrice triangolare superiore;
%       b := [n x 1] termine noto del sistema.
%
%   OUTPUT
%       x := [n x 1] soluzione del sistema.

toll=10^-9;
if norm(U-triu(U))>toll
    error('La matrice deve essere triangolare superiore')
end
if min(abs(diag(U)))==0
    error('matrice singolare')
end
n=size(U,1);
x=zeros(n,1);
x(n)=b(end)/U(n,n);
for k=1:n-1
    x(n-k)=(b(n-k)-U(n-k,:)*x)./U(n-k,n-k);
end