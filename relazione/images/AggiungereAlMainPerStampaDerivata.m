z=zeros(M,1);
figure(2)
plot (xsample,z,"k-")
hold on
plot(xsample,df(xsample),Color="#0072BD")
plot(zero,df(zero),"r.",LineWidth=2)
plot(0.6,df(0.6),".",LineWidth=2,Color="#77AC30")
title("Grafico della derivata f^~_α(x)")
legend("asse x","df^~_α(x)","ε","x_0",Location="northwest")
xlabel("x")
ylabel("y")
hold off