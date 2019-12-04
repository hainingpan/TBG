function [int1,int2]=LDOS_TMD_r(parameters)
%AA:r=(0,0)
%AB:r=(1/sqrt(3) aM,0) aM=a/theta
rAA=[0,0];
rAB=[1/sqrt(3),0]*parameters.a/parameters.theta;
eta=1e-3;

n=10;
xrange=-n:n;
yrange=-n:n;
bM1=parameters.bM1;
kp=parameters.kp;
kn=parameters.kn;
a1=[bM1(1)/2,3*kp(2)]/(2*n);
a2=[bM1(1)/2,3*kn(2)]/(2*n);
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
    end
end

enlist1=linspace(min(enmap(:,:,1),[],'all'),max(enmap(:,:,1),[],'all'),50);
enlist2=linspace(min(enmap(:,:,2),[],'all'),max(enmap(:,:,2),[],'all'),50);

ldosAA1=zeros(1,length(enlist1));
ldosAA2=zeros(1,length(enlist2));

parfor i=1:length(enlist1)
%     fprintf("i_r=%d of %d\n",i,length(enlist));
    deltaf=eta./((enlist1(i)-enmap).^2+eta^2);
    ldosprod=psiAA2.*deltaf;
    ldosAA1(i)=sum(ldosprod(:));
end 


parfor i=1:length(enlist2)
%     fprintf("i_r=%d of %d\n",i,length(enlist));
    deltaf=eta./((enlist2(i)-enmap).^2+eta^2);
    ldosprod=psiAA2.*deltaf;
    ldosAA2(i)=sum(ldosprod(:));
end 

% intAA=sum(ldosAA1(:))*(enlist1(2)-enlist1(1))/(sum(ldosAA2(:))*(enlist2(2)-enlist2(1)));
%
ldosAB1=zeros(1,length(enlist1));
ldosAB2=zeros(1,length(enlist2));

parfor i=1:length(enlist1)
%     fprintf("i_r=%d of %d\n",i,length(enlist));
    deltaf=eta./((enlist1(i)-enmap).^2+eta^2);
    ldosprod=psiAB2.*deltaf;
    ldosAB1(i)=sum(ldosprod(:));
end
parfor i=1:length(enlist2)
%     fprintf("i_r=%d of %d\n",i,length(enlist));
    deltaf=eta./((enlist2(i)-enmap).^2+eta^2);
    ldosprod=psiAB2.*deltaf;
    ldosAB2(i)=sum(ldosprod(:));
end

% intAB=sum(ldosAB1(:))*(enlist1(2)-enlist1(1))/(sum(ldosAB2(:))*(enlist2(2)-enlist2(1)));


int1=sum(ldosAA1(:))*(enlist1(2)-enlist1(1))/(sum(ldosAB1(:))*(enlist1(2)-enlist1(1)));

int2=sum(ldosAA2(:))*(enlist2(2)-enlist2(1))/(sum(ldosAB2(:))*(enlist2(2)-enlist2(1)));

% dos=zeros(1,length(enlist));
% parfor i=1:length(enlist)
%     deltaf=eta./((enlist(i)-enmap).^2+eta^2);
%     dos(i)=sum(deltaf(:));
% end

% intAA=sum(ldosAA(min(enmap(:,:,1),[],'all')<=enlist & enlist<=max(enmap(:,:,1),[],'all')))...
%     /sum(ldosAA(min(enmap(:,:,2),[],'all')<=enlist & enlist<=max(enmap(:,:,2),[],'all')));
% intAB=sum(ldosAB(min(enmap(:,:,1),[],'all')<=enlist & enlist<=max(enmap(:,:,1),[],'all')))...
%     /sum(ldosAB(min(enmap(:,:,2),[],'all')<=enlist & enlist<=max(enmap(:,:,2),[],'all')));

end
