clear
close all
clc

load dataset_chebtikonov.mat

n=12; M=length(xsample); lambda=1/M;
cstar=cheb_tikonov(n,lambda,xsample,ysample);

D=ChebDiff1D(n);
f=@(x) cstar'*cheb_vand(n,x)';
df=@(x) (D*cstar)'*cheb_vand(n,x)';
d2f=@(x) (D*D*cstar)'*cheb_vand(n,x)';

[zero,res,iterates,flag]=Newton(df,d2f,0,1e-15,100,'m');

fprintf("Il punto di minimo è %.10f con valore %.10f, La derivata seconda in quel punto vale %.10f.\nEssendo la derivata sdeconda strettamente positiva possiamo affermare che si tratta di un punto di minimo.\n",zero,f(zero),d2f(zero));

figure(1)
plot(xsample,ysample)
hold on
plot(xsample,f(xsample),"r-",LineWidth=2)
plot(zero,f(zero),'go',LineWidth=2)
title("Interpolazione ai minimi quadrati di una funzione affetta da rumore")
legend("funzione rumorosa - f(x)","funzione interpolante - f^~_α(x)","punto di minimo")
xlabel("x")
ylabel("f(x)")
hold off