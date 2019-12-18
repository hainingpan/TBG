function m=number(energyall,wfall,parameters)
%m=<..>=[SAx,SBx,SCx;SAy,SBy,SCy;SAz,SBz,SCz]
energyall_sort=sort(energyall(:));
mu=energyall_sort(end/2);
occupied=energyall<=mu;
wf_occupied=[wfall{occupied}];
N=length(wf_occupied);
sigma(:,:,1)=eye(2);
% sigma(:,:,2)=[0,-1i;1i,0];
% sigma(:,:,3)=[1,0;0,-1];
m=zeros(1,3);
Q=[0,0;parameters.kn-parameters.kp;parameters.kp-parameters.kn];
r=[0,0;sqrt(3)/2*parameters.aM,1/2*parameters.aM;sqrt(3)*parameters.aM,parameters.aM];
for i=1:1   %i=x,y,z
    for j=1:3   %j=A,B,C
        re=0;
        for k=1:N
            wf=reshape(wf_occupied(:,k),3,2);
            for q=1:3   %q for bra
                for q2=1:3  %q2 for momentum of bra-ket
                    bra=wf(q,:);
                    bra=bra.';
                    ket=wf(mod(q+q2-1-1,3)+1,:);
                    ket=ket.';
                    re=re+(bra'*sigma(:,:,i)*ket)*exp(1i*Q(q2,:)*r(j,:)');
                end
            end
        end
        m(i,j)=re/N;
    end
end
end