function [ldos,enlist]=LDOS_TMD_rx(parameters)
%AA:r=(0,0)
%AB:r=(1/sqrt(3) aM,0) aM=a/theta
rlistx=linspace(0,sqrt(3),100)*parameters.a/parameters.theta;
eta=1e-3;

n=20;
xrange=-n:n;
yrange=-n:n;
bM1=parameters.bM1;
bM2=parameters.bM2;
kp=parameters.kp;
kn=parameters.kn;
a1=-bM1/(2*n);
a2=(bM1+bM2)/(2*n);
enmap=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);


Nx=length(xrange);
Ny=length(yrange);

for xindex=1:Nx
    kx=xrange(xindex);
    for yindex=1:Ny
        ky=yrange(yindex);
        k=kx*a1+ky*a2;
        [val,vec]=energyTMD(k(1),k(2),parameters);
        enmap(xindex,yindex,:)=val;
        for i=1:length(rlistx)
            psir{i}(xindex,yindex,:)=u2(vec,[rlistx(i),0],parameters);
        end
    end
end

enlist=linspace(min(enmap(:,:,4),[],'all'),max(enmap(:,:,1),[],'all'),100);

ldos=zeros(length(rlistx),length(enlist));
for i=1:length(rlistx)
    for j=1:length(enlist)
        deltaf=eta./((enlist(j)-enmap).^2+eta^2);
        ldosprod=psir{i}.*deltaf;
        ldos(i,j)=sum(ldosprod(:));
    end
end
end
        

% dos=zeros(1,length(enlist));
% parfor i=1:length(enlist)
%     deltaf=eta./((enlist(i)-enmap).^2+eta^2);
%     dos(i)=sum(deltaf(:));
% end

