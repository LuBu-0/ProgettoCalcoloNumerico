clear
close all
clc
%1
load dataset_chebtikonov.mat
%2
n=12; M=length(xsample); lambda=1/M;
cstar=cheb_tikonov(n,lambda,xsample,ysample);
%3
D=ChebDiff1D(n);
f=@(x) cstar*cos(acos(x)*(0:n))';
df=@(x) (D*cstar')'*cos(acos(x)*(0:n))';
d2f=@(x) (D*D*cstar')'*cos(acos(x)*(0:n))';
%4
[zero,res,iterates,flag]=Newton(df,d2f,0,1e-15,100,'m');
%5
fprintf("Il punto di minimo è %d con valore %d, La derivata seconda in quel punto vale %d.\nEssendo la derivata sdeconda strettamente positiva possiamo affermare che si tratta di un punto di minimo.",zero,f(zero),d2f(zero));

figure(1)
plot(xsample,ysample)
hold on
plot(xsample,f(xsample),"r-",LineWidth=2)
plot(zero,f(zero),'go',LineWidth=2)
title("Interpolazione ai minimi quadrati di una funzione affetta da rumore")
legend("f(x)","f^~_α(x)","punto di minimo")
hold off

