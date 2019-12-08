function hubbardU(R)
[rx,ry]=meshgrid(linspace(-2*sqrt(3)*parameters.aM,2*sqrt(3)*parameters.aM,101));
[wbgrid,wtgrid]=w([0,0],rx,ry,parameters);
alpha=0.00729735;
