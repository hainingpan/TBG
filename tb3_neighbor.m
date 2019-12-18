function [Neighborlist,tmat]=tb3_neighbor(neighborlist,t,k,parameters)
Neighbor0list={[0,0]};
Neighbor1list={[1,0],[0,1],[-1,0],[0,-1]};
Neighbor2list={[1,1],[-1,1],[-1,-1],[1,-1]};
% Neighbor3list={[3,0],[2,1],[1,2],[0,3],[-1,2],[-2,1],[-3,0],[-2,-1],[-1,-2],[0,-3],[1,-2],[2,-1]};

% tmat0=[t0,t1(6),t3(6);
%        t1(3),t0,t1(6);
%        t3(3),t1(3),t0];
% tmat1{1}=[t2(5),0,0;t1(5),t2(5),0;t1(4),t1(5),t2(5)];
% tmat1{2}=[t2(6),0,0;0,t2(6),0;0,0,t2(6)];
% tmat1{3}=[t2(2),t1(2),t1(1);
%           0,t2(2),t1(2);
%           0,0,t2(2)];
% tmat1{4}=[t2(3),t1(4),t1(5);
%           0,t2(3),t1(4);
%           0,0,t2(3)];
% tmat2{1}=[0,0,0;t3(6),0,0,;t1(6),t3(6),]
Neighborlist=[Neighbor0list,Neighbor1list,Neighbor2list];
r={[0,0],[1/3,1/3],[2/3,2/3]};
% k=3;
neighbor=[neighborlist{1:k+1}];
cor=cell2mat(neighbor')*[parameters.aM1;parameters.aM2];
t_all=[t{1:k+1}];
tmat={};
for nn=1:length(Neighborlist)
    for i=1:3
        for j=1:3
            coor=Neighborlist{nn}+r{j}-r{i};
            pos=coor(1)*parameters.A1+coor(2)*parameters.A2;
            dev=sum(abs(pos-cor),2);
            [minval,ind]=min(dev);
            if minval<1e-4
                tmat{nn}(i,j)=t_all(ind);
            else
                tmat{nn}(i,j)=0;
            end
        end
    end
end
end

         
         
 