function [dos1,dos2,ldosAA1,ldosAA2,ldosAB1,ldosAB2,enlist1,enlist2,int1,int2]=LDOS_TMD_r(parameters)
%AA:r=(0,0)
%AB:r=(1/sqrt(3) aM,0) aM=a/theta
rAA=[0,0];
rAB=[1/sqrt(3),0]*parameters.a/parameters.theta;
eta=1e-3;

n=10;
xrange=-n:n;
yrange=-n:n;
bM1=parameters.bM1;
bM2=parameters.bM2;
% kp=parameters.kp;
% kn=parameters.kn;
a1=bM2/(2*n);
a2=(-2*bM1-bM2)/2/(2*n);
enmap=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
psiAA2=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
psiAB2=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);


Nx=length(xrange);
Ny=length(yrange);

for xindex=1:Nx
    kx=xrange(xindex);
    for yindex=1:Ny
        ky=yrange(yindex);
        k=kx*a1+ky*a2;
        [val,vec]=energyTMD(k(1),k(2),parameters);
        enmap(xindex,yindex,:)=val;
        psiAA2(xindex,yindex,:)=u2(vec,rAA,parameters);      
        psiAB2(xindex,yindex,:)=u2(vec,rAB,parameters); 
    end
end

enlist1=linspace(min(enmap(:,:,1),[],'all'),max(enmap(:,:,1),[],'all'),50);
enlist2=linspace(min(enmap(:,:,2),[],'all'),max(enmap(:,:,2),[],'all'),50);

ldosAA1=zeros(1,length(enlist1));
ldosAA2=zeros(1,length(enlist2));
for i=1:length(enlist1)
    deltaf=eta./((enlist1(i)-enmap).^2+eta^2);
    ldosprod=psiAA2.*deltaf;
    ldosAA1(i)=sum(ldosprod(:));
end 
for i=1:length(enlist2)
    deltaf=eta./((enlist2(i)-enmap).^2+eta^2);
    ldosprod=psiAA2.*deltaf;
    ldosAA2(i)=sum(ldosprod(:));
end 

% intAA=sum(ldosAA1(:))*(enlist1(2)-enlist1(1))/(sum(ldosAA2(:))*(enlist2(2)-enlist2(1)));
ldosAB1=zeros(1,length(enlist1));
ldosAB2=zeros(1,length(enlist2));

for i=1:length(enlist1)
    deltaf=eta./((enlist1(i)-enmap).^2+eta^2);
    ldosprod=psiAB2.*deltaf;
    ldosAB1(i)=sum(ldosprod(:));
end
for i=1:length(enlist2)
    deltaf=eta./((enlist2(i)-enmap).^2+eta^2);
    ldosprod=psiAB2.*deltaf;
    ldosAB2(i)=sum(ldosprod(:));
end
% intAB=sum(ldosAB1(:))*(enlist1(2)-enlist1(1))/(sum(ldosAB2(:))*(enlist2(2)-enlist2(1)));
int1=sum(ldosAA1(:))*(enlist1(2)-enlist1(1))/(sum(ldosAB1(:))*(enlist1(2)-enlist1(1)));
int2=sum(ldosAA2(:))*(enlist2(2)-enlist2(1))/(sum(ldosAB2(:))*(enlist2(2)-enlist2(1)));

dos1=zeros(1,length(enlist1));
for i=1:length(enlist1)
    deltaf=eta./((enlist1(i)-enmap).^2+eta^2);
    dos1(i)=sum(deltaf(:));
end
dos2=zeros(1,length(enlist2));
for i=1:length(enlist2)
    deltaf=eta./((enlist2(i)-enmap).^2+eta^2);
    dos2(i)=sum(deltaf(:));
end

end
