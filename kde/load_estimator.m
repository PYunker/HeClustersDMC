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
