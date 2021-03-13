# load the estimate
fid = fopen('estimate','r');
grid_n = fread(fid,1,'int32');
Ef = fread(fid,grid_n,'float64');
Ef2 = fread(fid,grid_n,'float64');
Es = fread(fid,grid_n,'float64');
Es2 = fread(fid,grid_n,'float64');
GRID = fread(fid,grid_n,'float64');
exact = fread(fid,grid_n,'float64');
fclose(fid)
# load the estimator
fid = fopen('estimator','r');
d = fread(fid,1,'int32');
gam = fread(fid,1,'int32');
runs = fread(fid,1,'int32');
mu = reshape(fread(fid,d*gam*runs,'float64'),[d,gam*runs]);
w = fread(fid,gam*runs,'float64');
sig = reshape(fread(fid,d*d*gam*runs,'float64'),[d,d,gam*runs]);
det = fread(fid,gam*runs,'float64');
fclose(fid)
# plot AKDE and KDE mean estimate and true PDF
figure 1
hold on
plot(GRID,Ef,"-b;AKDE;",'linewidth',2)
plot(GRID,Es,"-r;KDE;",'linewidth',2)
plot(GRID,exact,"--k;Presne;",'linewidth',2)
# PLOT AKDE and KDE bias and variance
figure 2
hold on
plot(GRID, (Ef - exact).^2, "-b;AKDE bias;",'linewidth',2)
plot(GRID, (Es - exact).^2, "-r;KDE bias;",'linewidth',2)
plot(GRID, Ef2 - Ef.^2, "--b;AKDE variance;",'linewidth',2)
plot(GRID, Es2 - Es.^2, "--r;KDE variance;",'linewidth',2)
