function [val,vec]=energyTBG(kx,ky,parameters)
%energy at (kx,ky)
Tmat=parameters.Tmat;
TTmat=parameters.TTmat;

theta=parameters.theta;
vf=parameters.vf;
sigma0=parameters.sigma0;
sigmax=parameters.sigmax;
sigmay=parameters.sigmay;
sigmaz=parameters.sigmaz;

bM1=parameters.bM1;
bM2=parameters.bM2;
kb=parameters.kb;
kt=parameters.kt;
h1index=parameters.h1index;
h2index=parameters.h2index;

klist=[kx,ky]+h1index(:)*bM1+h2index(:)*bM2;

hb=arrayfun(@(x,y) vf*(cos(theta/2)*sigma0-1i*sigmaz*sin(theta/2))*((x-kb(1))*sigmax+(y-kb(2))*sigmay),klist(:,1),klist(:,2),'UniformOutput',false) ;    %bottom takes value of +1
Hb=blkdiag(hb{:});
ht=arrayfun(@(x,y) vf*(cos(theta/2)*sigma0+1i*sigmaz*sin(theta/2))*((x-kt(1))*sigmax+(y-kt(2))*sigmay),klist(:,1),klist(:,2),'UniformOutput',false) ;    %top takes value of -1
Ht=blkdiag(ht{:});
H=[Hb,Tmat;TTmat,Ht];
[vec,val]=eig(H);
val=diag(val);
val=val(end:-1:1);
vec=vec(:,end:-1:1);
end


