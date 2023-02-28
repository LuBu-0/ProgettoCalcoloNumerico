function x=SostituzioneAvanti(L,b)
%SostituzioneAvanti risolve il sistema lineare Lx = b mediante l'algoritmo
%di sostituzioine avanti
%   
%   INPUT:
%       L := [n x n] matrice triangolare inferiore;
%       b := [n x 1] termine noto del sistema.
%
%   OUTPUT
%       x := [n x 1] soluzione del sistema.
toll=10^-12;
if norm(L-tril(L))>toll
    error('La matrice deve essere triangolare inferiore')
end
if min(abs(diag(L)))==0
    error('matrice singolare')
end
n=size(L,1);
x=zeros(n,1);
x(1)=b(1)/L(1,1);
for k=1:n-1
    x(k+1)=(b(k+1)-L(k+1,:)*x)./L(k+1,k+1);
end