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
% S=zeros(2,2,3);
% S=ones(2,2,3);
m0=[1,0,0;cos(-2*pi/3),sin(-2*pi/3),0;cos(-4*pi/3),sin(-4*pi/3),0]';
%  m0=[1,0,0;0,1,0;0,-1,0]';
S0=init(m0,parameters);
S=S0;
clear mu SS m total
theta=zeros(100,3);
sm=zeros(100,3);

i=0;
gap=0;
Ux=U{1}/20;
while ((length(gap)<=2 || abs(gap(i)-gap(i-1))>1e-5) && i<1000)
    i=i+1;
    [energyall,wfall]=energyMF(kx3list,ky3list,S,[neighborlist{1:k+1}],[t{1:k+1}],Ux,parameters);
    [S,mu(i)]=S_calc(energyall,wfall);
    SS{i}=S;
    m(:,:,i)=spin(energyall,wfall,parameters);
    theta(i,:)=orientation(m(:,:,i));
    sm(i,:)=sum(abs(m(:,:,i)).^2,1);
    gap(i)=min(1e3*energyall(:,4))-max(1e3*energyall(:,3));
    if i>=2
%         total(i)=1000*real(total_energy(energyall,wfall,energyall2,wfall2,Ux,parameters));
        total(i)=1000*real(total_energy_def(kx3list,ky3list,energyall,wfall,[neighborlist{1:k+1}],[t{1:k+1}],Ux,parameters));
        fprintf("Iter:%d Gap:%f E:%f Spin:(%f,%f,%f) |S|:(%f,%f,%f)\n",i,gap(i),total(i),theta(i,:),sm(i,:));
    end
%     energyall=energyall;
%     wfall=wfall;
end
