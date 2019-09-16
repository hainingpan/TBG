function [kcxmap,kcymap,kcx2map,kcy2map,bcmap]=berrycurvature(level,parameters)
n=90;
bM1=parameters.bM1;
kp=parameters.kp;
kn=parameters.kn;
a1=[bM1(1)/2,3*kp(2)]/(2*n);
a2=[bM1(1)/2,3*kn(2)]/(2*n);
xrange=-n:n;
yrange=-n:n;
kxmap=zeros(2*n+1,2*n+1);
kymap=zeros(2*n+1,2*n+1);
kx2map=zeros(2*n+1,2*n+1); %after shifting to Hexagon
ky2map=zeros(2*n+1,2*n+1);  %after shifting to Hexagon
bcmap=zeros(2*n+1,2*n+1);
umap=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
omega=abs(cross([a1,0],[a2,0]));
omega=omega(3);
parfor xindex=1:length(xrange)
    kx=xrange(xindex);
    kxlist=zeros(1,2*n+1);
    kylist=zeros(1,2*n+1);
    kx2list=zeros(1,2*n+1);
    ky2list=zeros(1,2*n+1);
    ulist=zeros(2*n+1,2*(2*parameters.Nmax+1)^2);
    for yindex=1:length(yrange)        
        ky=yrange(yindex);
        k=kx*a1+ky*a2;
        kxlist(yindex)=k(1);
        kylist(yindex)=k(2);
        shift=[0,0];
        if (k(1)<=0) && (k(2)>=2*kn(2)/bM1(1)*k(1)+2*kn(2))
            shift=a1*2*n;
        end
        if (k(1)>=0) && (k(2)>=-2*kn(2)/bM1(1)*k(1)+2*kn(2))
            shift=-a2*2*n;
        end
        if (k(1)<=0) && (k(2)<=-2*kn(2)/bM1(1)*k(1)-2*kn(2))
            shift=a2*2*n;
        end
        if (k(1)>=0) && (k(2)<=2*kn(2)/bM1(1)*k(1)-2*kn(2))
            shift=-a1*2*n;
        end        
        kx2list(yindex)=k(1)+shift(1);
        ky2list(yindex)=k(2)+shift(2);
        [~,vec]=energy(k(1),k(2),parameters);
        ulist(yindex,:)=vec(:,level);
    end
    kxmap(xindex,:)=kxlist;
    kymap(xindex,:)=kylist;
    kx2map(xindex,:)=kx2list;
    ky2map(xindex,:)=ky2list;
    umap(xindex,:,:)=ulist;
end


kcxmap=(kxmap(1:end-1,1:end-1)+kxmap(1:end-1,2:end)+kxmap(2:end,1:end-1)+kxmap(2:end,2:end))/4;
kcymap=(kymap(1:end-1,1:end-1)+kymap(1:end-1,2:end)+kymap(2:end,1:end-1)+kymap(2:end,2:end))/4;
kcx2map=(kx2map(1:end-1,1:end-1)+kx2map(1:end-1,2:end)+kx2map(2:end,1:end-1)+kx2map(2:end,2:end))/4;
kcy2map=(ky2map(1:end-1,1:end-1)+ky2map(1:end-1,2:end)+ky2map(2:end,1:end-1)+ky2map(2:end,2:end))/4;


bcmap=angle(dot(umap(1:end-1,1:end-1,:),umap(1:end-1,2:end,:),3)...
    .*dot(umap(1:end-1,2:end,:),umap(2:end,2:end,:),3)...
    .*dot(umap(2:end,2:end,:),umap(2:end,1:end-1,:),3)...
    .*dot(umap(2:end,1:end-1,:),umap(1:end-1,1:end-1,:),3))/omega;
    
% parfor xindex=1:length(xrange)
%     x=xrange(xindex);
%     kxlist=zeros(1,2*n+1);
%     kylist=zeros(1,2*n+1);
%     kx2list=zeros(1,2*n+1);
%     ky2list=zeros(1,2*n+1);
%     bclist=zeros(1,2*n+1);
%     for yindex=1:length(yrange)        
%         y=xrange(yindex);
%         k=x*a1+y*a2;
%         right=k+[bM1(1)/(4*n),0];
%         left=k+[-bM1(1)/(4*n),0];
%         up=k+[0,3*kn(2)/(2*n)];
%         down=k+[0,3*kp(2)/(2*n)];
%         kxlist(yindex)=k(1);
%         kylist(yindex)=k(2);
%         shift=[0,0];
%         if (k(1)<=0) && (k(2)>=2*kn(2)/bM1(1)*k(1)+2*kn(2))
%             shift=a1*2*n;
%         end
%         if (k(1)>=0) && (k(2)>=-2*kn(2)/bM1(1)*k(1)+2*kn(2))
%             shift=-a2*2*n;
%         end
%         if (k(1)<=0) && (k(2)<=-2*kn(2)/bM1(1)*k(1)-2*kn(2))
%             shift=a2*2*n;
%         end
%         if (k(1)>=0) && (k(2)<=2*kn(2)/bM1(1)*k(1)-2*kn(2))
%             shift=-a1*2*n;
%         end        
%         kx2list(yindex)=k(1)+shift(1);
%         ky2list(yindex)=k(2)+shift(2);
%         bclist(yindex)=wilson(level,right,up,left,down,parameters);      
%     end
%     kxmap(xindex,:)=kxlist;
%     kymap(xindex,:)=kylist;
%     kx2map(xindex,:)=kx2list;
%     ky2map(xindex,:)=ky2list;
%     bcmap(xindex,:)=bclist;
% end
% kcxmap=kxmap;
% kcymap=kymap;
% kcx2map=kx2map;
% kcy2map=ky2map;
end