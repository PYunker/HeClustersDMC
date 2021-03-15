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

# plot AKDE and KDE mean estimate and true PDF
fig1 = figure(1);
hold on
plot(GRID,Ef,"-b;AKDE;",'linewidth',2)
plot(GRID,Es,"-r;KDE;",'linewidth',2)
plot(GRID,exact,"--k;exact ;",'linewidth',2)

# PLOT AKDE and KDE bias and variance
fig2 = figure(2);
hold on
plot(GRID, (Ef - exact).^2, "-b;AKDE bias;",'linewidth',2)
plot(GRID, (Es - exact).^2, "-r;KDE bias;",'linewidth',2)
plot(GRID, Ef2 - Ef.^2, "--b;AKDE variance;",'linewidth',2)
plot(GRID, Es2 - Es.^2, "--r;KDE variance;",'linewidth',2)

# hang until figures are closed, so that script doesn't terminate prematurely
# when invoked directly from command line. when invoking thrugh Octave, it's not
# needed
waitfor(fig1); waitfor(fig2)
