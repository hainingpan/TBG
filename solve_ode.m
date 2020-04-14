param=mainTri('W',27*3,'L',ceil(27/sqrt(3)),'left',8*sqrt(3)-5,'right',8*sqrt(3)+5,'phi',5/4*pi);
param=construct_triangular(param);
param=initialize(1,param);
eq=@(t,s) ode(t,s,param);
tlist=0:0.1:10;
[t,s]=ode45(eq,tlist,param.s0);
param.t=t;
param.s=s;