function re=Htb_tri(parameters)
vf=parameters.vf;
NN=parameters.NN;
a0=parameters.a0;
sigma0=parameters.sigma0;
sigmax=parameters.sigmax;
sigmay=parameters.sigmay;
sigmaz=parameters.sigmaz;
bM1=parameters.bM1;
bM2=parameters.bM2;
T0=parameters.T0;
T1=parameters.T1;
Tm1=parameters.Tm1;
theta=parameters.theta;
kb=parameters.kb;
kt=parameters.kt;
aM=parameters.a/parameters.theta;
n=parameters.n;
parameters.a0=aM/sqrt(3)/n;




ilist=1:NN;
jlist=1:NN;
[imat,jmat]=meshgrid(ilist,jlist);
banddiff=spdiags(ones(NN),-1,NN,NN)-spdiags(ones(NN),1,NN,NN);
kxmat=kron(banddiff,speye(NN));
kymat=kron(speye(NN),banddiff);
Hb=vf*(-1i/(2*a0)*kron((cos(theta/2)*sigma0-1i*sigmaz*sin(theta/2))*sigmax,kxmat)...
       -1i/(2*a0)*kron((cos(theta/2)*sigma0-1i*sigmaz*sin(theta/2))*sigmay,kymat)...
       -kron((cos(theta/2)*sigma0-1i*sigmaz*sin(theta/2))*(kb(1)*sigmax+kb(2)*sigmay),speye(NN^2)));
Ht=vf*(-1i/(2*a0)*kron((cos(theta/2)*sigma0+1i*sigmaz*sin(theta/2))*sigmax,kxmat)...
       -1i/(2*a0)*kron((cos(theta/2)*sigma0+1i*sigmaz*sin(theta/2))*sigmay,kymat)...
       -kron((cos(theta/2)*sigma0+1i*sigmaz*sin(theta/2))*(kt(1)*sigmax+kt(2)*sigmay),speye(NN^2)));
expmat1=exp(-1i*(bM1(1)*a0*imat+bM1(2)*a0*jmat)).';
expmat2=exp(-1i*(bM2(1)*a0*imat+bM2(2)*a0*jmat)).';
T=kron(T0,speye(NN^2))+kron(T1,spdiags(expmat1(:),0,NN^2,NN^2))...
    +kron(Tm1,spdiags(expmat2(:),0,NN^2,NN^2));
re=[Hb,T;T',Ht];
% re=Hb;
end