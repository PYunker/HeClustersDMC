load data
nb = 50;
x=linspace(-5,5,nb);
[X,Y]=meshgrid(x,x);
grid = [X(:),Y(:)]';
Z = reshape(Ef,[nb,nb]);
y = sqrt(0.5/pi) * exp(-0.5*x.^2);
Z2 = y' * y;

figure 1;
contour(X,Y,Z2,0.02:0.02:0.1,"-r;presne;",'linewidth',2)
hold on
contour(X,Y,Z,0.02:0.02:0.1,"--b;akde;",'linewidth',2)

figure 2;
bias = (Z - Z2).^2;
variance = reshape(Ef2,[nb,nb])-Z.^2;
subplot(1,2,1)
surf(X,Y,bias)
title('Bias')
subplot(1,2,2)
surf(X,Y,variance)
title('Variance')

# ===== monte carlo integrace ==========
#{
subplot(1,2,1)
surf(X,Y,Z)
subplot(1,2,2)
surf(X,Y,Z2)
#}
