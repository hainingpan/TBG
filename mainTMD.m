function parameters=mainTMD()
parameters=struct('a',3.472e-4*5.076,'m',0.62*0.511e6,'theta',1.2/360*2*pi,'V',8e-3,'psi',-89.6/360*2*pi,'w',-8.5e-3,...
    'Nmax',5);
%Reciprocal lattice
parameters.G1=[0,4*pi/(sqrt(3)*parameters.a)];
parameters.G2=[-2*pi/parameters.a,(2*pi)/(sqrt(3)*parameters.a)];
parameters.G3=[-2*pi/parameters.a,-(2*pi)/(sqrt(3)*parameters.a)];
parameters.G4=[0,-4*pi/(sqrt(3)*parameters.a)];
parameters.G5=[2*pi/parameters.a,-(2*pi)/(sqrt(3)*parameters.a)];
parameters.G6=[2*pi/parameters.a,(2*pi)/(sqrt(3)*parameters.a)];
parameters.b1=parameters.G1;    %reciprocal unit vector of small lattice 
parameters.b2=parameters.G5;    %reciprocal unit vector of small lattice
parameters.bM1=parameters.theta.*[parameters.b1(2),-parameters.b1(1)];  %bM1=theta* cross(b1,z); reciprocal unit vector of Moire lattice
parameters.bM2=parameters.theta.*[parameters.b2(2),-parameters.b2(1)];  %bM2=theta* cross(b2,z); reciprocal unit vector of Moire lattice
parameters.kp=[-parameters.bM1(1)/2,-parameters.bM1(1)/2/sqrt(3)];  %k+=[-bM1.x/2,-bM1.x/2 tan(pi/6)]
parameters.kn=[-parameters.bM1(1)/2,parameters.bM1(1)/2/sqrt(3)]; %k-=[-bM1.x/2,bM1.x/2 tan(pi/6)]

%Potential term
Nrange=-parameters.Nmax:parameters.Nmax;
[h1index,h2index]=meshgrid(Nrange,Nrange);
parameters.h1index=h1index;
parameters.h2index=h2index;
[h1matX,h1matY]=meshgrid(h1index(:));
[h2matX,h2matY]=meshgrid(h2index(:));
h1mat=h1matX-h1matY;
h2mat=h2matX-h2matY;
parameters.DeltaTmat=reshape(arrayfun(@(h1,h2) DeltaT(h1,h2,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
parameters.DeltaTTmat=reshape(arrayfun(@(h1,h2) DeltaTT(h1,h2,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
parameters.Deltatmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,-1,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
parameters.Deltabmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,1,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
end


