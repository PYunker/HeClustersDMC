load data1d
nb = 500;
x=linspace(-5,5,nb);
y = sqrt(0.5/pi) * exp(-0.5*x.^2);

figure 1;
plot(x,Ef,"-r;presne;",'linewidth',2)
hold on
plot(x,y,"--b;akde;",'linewidth',2)