function [energyall,wfall]=energyMF(kxlist,kylist,S,bond,t,U,parameters)
%S(a,b,c): a,b=1(up) or 2(down), c=1,2,3 (0, Q, -Q)
N=length(kxlist);
Q=parameters.kn-parameters.kp;
Ez=parameters.Ez;
kx1list=kxlist; %q=0
ky1list=kylist;
kx2list=kxlist+Q(1);
ky2list=kylist+Q(2);    %q=Q
kx3list=kxlist-Q(1);
ky3list=kylist-Q(2);    %q=-Q
energylist=real(tb(bond,t,[kx1list,kx2list,kx3list,-kx1list,-kx2list,-kx3list],[ky1list,ky2list,ky3list,-ky1list,-ky2list,-ky3list],parameters));
% energylist=reshape(energylist,6,[]);
off_H=U*[S(1,1,1),S(1,1,3),S(1,1,2),S(1,2,1),S(1,2,3),S(1,2,2);
       S(1,1,2),S(1,1,1),S(1,1,3),S(1,2,2),S(1,2,1),S(1,2,3);
       S(1,1,3),S(1,1,2),S(1,1,1),S(1,2,3),S(1,2,2),S(1,2,1);
       S(2,1,1),S(2,1,3),S(2,1,2),S(2,2,1),S(2,2,3),S(2,2,2);
       S(2,1,2),S(2,1,1),S(2,1,3),S(2,2,2),S(2,2,1),S(2,2,3);
       S(2,1,3),S(2,1,2),S(2,1,1),S(2,2,3),S(2,2,2),S(2,2,1)];
% Zeeman=Ez*[ones(3),zeros(3);zeros(3),-ones(3)];
Zeeman=Ez*diag([1,1,1,-1,-1,-1]);
energyall=zeros(N,6);
wfall=cell(N,6);
for i=1:N
    H=diag(energylist(i,:))+off_H+Zeeman;
    [vec,val]=eig(H);
    val=real(diag(val));
    [val,I]=sort(val);
    vec=vec(:,I);
    energyall(i,:)=val;
    for ii=1:6
        wfall{i,ii}=vec(:,ii);
    end
end

end
    
    


