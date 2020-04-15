function re=total_energy_def(kxlist,kylist,energyall,wfall,bond,t,Ux,parameters)
% ek{i,j} i=1 q=0, i=2 q=Q, i=3 q=-Q, j=1 up, j=2 down
energyall_sort=sort(energyall(:));
mu=energyall_sort(end/2);
occupied=energyall<=mu;
wf_occupied=[wfall{occupied}];
N=length(wf_occupied);
Q=parameters.kn-parameters.kp;
Ez=parameters.Ez;

kx1list=kxlist; %q=0
ky1list=kylist;
kx2list=kxlist+Q(1);
ky2list=kylist+Q(2);    %q=Q
kx3list=kxlist-Q(1);
ky3list=kylist-Q(2);    %q=-Q
energylist=real(tb(bond,t,[kx1list,kx2list,kx3list,-kx1list,-kx2list,-kx3list],[ky1list,ky2list,ky3list,-ky1list,-ky2list,-ky3list],parameters));
ek{1,1}=repmat(energylist(:,1),[1,6])+Ez;
ek{2,1}=repmat(energylist(:,2),[1,6])+Ez;
ek{3,1}=repmat(energylist(:,3),[1,6])+Ez;
ek{1,2}=repmat(energylist(:,4),[1,6])-Ez;
ek{2,2}=repmat(energylist(:,5),[1,6])-Ez;
ek{3,2}=repmat(energylist(:,6),[1,6])-Ez;

for i=1:3
    for j=1:2
        eko{i,j}=ek{i,j}(occupied);
    end
end
K=0;
for i=1:N
    wf=reshape(wf_occupied(:,i),3,2);
    for j=1:3
        for k=1:2
            K=K+eko{j,k}(i)*abs(wf(j,k))^2;
        end
    end
end
K=K/N;

c=cc(energyall,wfall,parameters);
re=K+Ux*(c(1,1)*c(2,2)-c(1,2)*c(2,1));
end


    
    
