function gap_Zeeman(epsilon,theta)

NVz=51;
Vzlist=linspace(0,50,NVz);
NB=41;
Blist=linspace(0,15,NB);

gapmap=zeros(NVz,NB);
szmap=zeros(NVz,NB,3);
isconverge=zeros(NVz,NB);
parfor Vzindex=1:NVz    
parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'Vz',Vzlist(Vzindex),'B',0);
[t,neighborlist]=t_calc_func(3,parameters);
U=U_calc_func(0,parameters);
Ux=U{1}/epsilon;

n=10;
counter=1;
N=3*n^2+3*n+1;
kx3list=zeros(N,1);
ky3list=zeros(N,1);
a1=-(2*parameters.bM1+parameters.bM2)/3/n;
a2=(parameters.bM1+2*parameters.bM2)/3/n;
a31=a1/sqrt(3)*[0,1;-1,0]; %rotate 90 deg clockwisely
a32=a2/sqrt(3)*[0,1;-1,0];  %rotate 90 deg clockwisely
for yindex=-n:n
    for xindex=max(-n,-n+yindex):min(n+yindex,n)
        k2=xindex*a31+yindex*a32;
        kx3list(counter)=k2(1);
        ky3list(counter)=k2(2);
        counter=counter+1;
    end
end

k=3;
m=[1,0,0.;cos(-2.*pi/3),sin(-2.*pi/3),0.;cos(-4.*pi/3),sin(-4.*pi/3),0.]';

Blist=linspace(0,10,NB);
    for Bindex=1:NB

    parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'Vz',Vzlist(Vzindex),'B',Blist(Bindex));
    S0=init(m,parameters);
    S=S0;
    itermax=2000;
    phi=zeros(itermax,3);

    i=0;
    gap=0;
    gap_old=1;
    while (abs(gap-gap_old)>1e-5 && i<itermax)
        i=i+1;
        [energyall,wfall]=energyMF(kx3list,ky3list,S,[neighborlist{1:k+1}],[t{1:k+1}],Ux,parameters);
        S=S_calc(energyall,wfall);
        m=spin(energyall,wfall,parameters);
        [~,phi]=orientation(m);
        gap_old=gap;
        gap=min(1e3*energyall(:,4))-max(1e3*energyall(:,3));        
    end
    szmap(Vzindex,Bindex,:)=cos(phi*pi/180);
    gapmap(Vzindex,Bindex)=gap;
    isconverge(Vzindex,Bindex)=i;
    end
end
save(strcat('gap_ep',num2str(epsilon),'theta',num2str(theta),'.mat'),'gapmap','szmap','Vzlist','Blist','isconverge');
end


