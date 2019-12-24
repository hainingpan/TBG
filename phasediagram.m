function [re,isconverge]=phasediagram(theta,Vz,epsilon)
parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'Vz',Vz);
[t,neighborlist]=t_calc_func(3,parameters);
U=U_calc_func(0,parameters);
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
re=zeros(1,length(epsilon));
isconverge=zeros(1,length(epsilon));
for epi=1:length(epsilon)
    Ux=epsilon(epi)*U{1};
    if Vz>=0
        m0=[1,0,0;cos(-2*pi/3),sin(-2*pi/3),0;cos(-4*pi/3),sin(-4*pi/3),0]';
    else
        m0=[1,0,0;cos(2*pi/3),sin(2*pi/3),0;cos(4*pi/3),sin(4*pi/3),0]';
    end
    S0=init(m0,parameters);
    S=S0;
    i=0;
    gap=0;
    while (length(gap)<=2 || abs(gap(i)-gap(i-1))>1e-5 || i>1e3)
        i=i+1;
        [energyall,wfall]=energyMF(kx3list,ky3list,S,[neighborlist{1:k+1}],[t{1:k+1}],Ux,parameters);
        [S,~]=S_calc(energyall,wfall);
    %     SS{i}=S;
    %     m(:,:,i)=spin(energyall,wfall,parameters);
    %     theta(i,:)=orientation(m(:,:,i));
        gap(i)=min(1e3*energyall(:,4))-max(1e3*energyall(:,3));
    %     fprintf("Iter:%d Gap:%f Spin:(%f,%f,%f)\n",i,gap(i),theta(i,:))
    end
    re(epi)=min(1e3*energyall(:,4))-max(1e3*energyall(:,3));
    isconverge(epi)=1-(i>1e4);
end
end


