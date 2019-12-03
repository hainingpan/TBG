function [ldosAA,ldosAB,enlist]=LDOS_TMD_r(parameters)
%AA:r=(0,0)
%AB:r=(1/sqrt(3) aM,0) aM=a/theta
rAA=[0,0];
rAB=[1/sqrt(3),0]*parameters.a/parameters.theta;

eta=1e-3;

n=5;
xrange=-n:n;
yrange=-n:n;
bM1=parameters.bM1;
kp=parameters.kp;
kn=parameters.kn;
a1=[bM1(1)/2,3*kp(2)]/(2*n);
a2=[bM1(1)/2,3*kn(2)]/(2*n);
kxmap=zeros(2*n+1,2*n+1);
kymap=zeros(2*n+1,2*n+1);
enmap=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
psiAA2=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
psiAB2=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);


Nx=length(xrange);
Ny=length(yrange);

parfor xindex=1:Nx
    kx=xrange(xindex);
    for yindex=1:Ny
        ky=yrange(yindex);
        k=kx*a1+ky*a2;
        [val,vec]=energyTMD(k(1),k(2),parameters);
        enmap(xindex,yindex,:)=val;
        psiAA2(xindex,yindex,:)=u2(vec,rAA,parameters);      
        psiAB2(xindex,yindex,:)=u2(vec,rAB,parameters); 
%         
%         kxmap(xindex,yindex)=k(1);
%         kymap(xindex,yindex)=k(2);
    end
end

enlist=linspace(min(enmap(:,:,4),[],'all'),max(enmap(:,:,1),[],'all'),100);
ldosAA=zeros(1,length(enlist));
parfor i=1:length(enlist)
%     fprintf("i_r=%d of %d\n",i,length(enlist));
    deltaf=eta./((enlist(i)-enmap).^2+eta^2);
    ldosprod=psiAA2.*deltaf;
    ldosAA(i)=sum(ldosprod(:));
end 
% 
ldosAB=zeros(1,length(enlist));
parfor i=1:length(enlist)
%     fprintf("i_r=%d of %d\n",i,length(enlist));
    deltaf=eta./((enlist(i)-enmap).^2+eta^2);
    ldosprod=psiAB2.*deltaf;
    ldosAB(i)=sum(ldosprod(:));
end

end
