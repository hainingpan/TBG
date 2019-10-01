function parameters=mainTBG()
parameters=struct('a',0.246e-3*5.076,'theta',1.2/360*2*pi,...
    'w0',0e-3,'w1',0e-3,'vf',1e6/3e8,...
    'Nmax',5,...
    'a0',0.246e-3*5.076*1,'NN',300);

%Pauli matrix
parameters.sigma0=eye(2);
parameters.sigmax=[0,1;1,0];
parameters.sigmay=[0,-1i;1i,0];
parameters.sigmaz=[1,0;0,-1];

%Reciprocal lattice
parameters.bM1= (4*pi/(sqrt(3)*parameters.a/parameters.theta))*[1/2,sqrt(3)/2];  %reciprocal unit vector of small lattice  b1=b+ in Wu,PRB, 2019
parameters.bM2= (4*pi/(sqrt(3)*parameters.a/parameters.theta))*[-1/2,sqrt(3)/2];   %reciprocal unit vector of small lattice b2=b- in Wu,PRB, 2019
parameters.kb= (4*pi/(3*parameters.a/parameters.theta))*[-sqrt(3)/2,-1/2];  %bottom takes l=1
parameters.kt= (4*pi/(3*parameters.a/parameters.theta))*[-sqrt(3)/2,1/2];  %top takes l=-1
%T_j=w0 sigma0+ w1 *(cos(2pi j/3)sigmax+sin(2pi/3)sigmay)
parameters.T0=parameters.w0*parameters.sigma0+parameters.w1*(parameters.sigmax);
parameters.T1=parameters.w0*parameters.sigma0+parameters.w1*(cos(2*pi/3)*parameters.sigmax+sin(2*pi/3)*parameters.sigmay);
parameters.Tm1=parameters.w0*parameters.sigma0+parameters.w1*(cos(2*pi/3)*parameters.sigmax-sin(2*pi/3)*parameters.sigmay);
parameters.kp=(parameters.bM1+parameters.bM2)/3;
parameters.kn=parameters.kp*[cos(pi/3),-sin(pi/3);sin(pi/3),cos(pi/3)];


%Potential term
Nrange=-parameters.Nmax:parameters.Nmax;
[h1index,h2index]=meshgrid(Nrange,Nrange);
parameters.h1index=h1index;
parameters.h2index=h2index;
[h1matX,h1matY]=meshgrid(h1index(:));
[h2matX,h2matY]=meshgrid(h2index(:));
h1mat=h1matX-h1matY;
h2mat=h2matX-h2matY;
parameters.Tmat=cell2mat(reshape(arrayfun(@(h1,h2) T(h1,h2,parameters),h1mat(:),h2mat(:),'UniformOutput',false),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2));
parameters.TTmat=cell2mat(reshape(arrayfun(@(h1,h2) TT(h1,h2,parameters),h1mat(:),h2mat(:),'UniformOutput',false),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2));
end



