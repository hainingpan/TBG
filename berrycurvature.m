function [kxmap,kymap,kx2map,ky2map,bcmap]=berrycurvature(level,parameters)
n=40;
bM1=parameters.bM1;
kp=parameters.kp;
kn=parameters.kn;
a1=[bM1(1)/2,3*kp(2)]/(2*n);
a2=[bM1(1)/2,3*kn(2)]/(2*n);
xrange=-n:n;
yrange=-n:n;
kxmap=zeros(2*n+1,2*n+1);
kymap=zeros(2*n+1,2*n+1);
kx2map=zeros(2*n+1,2*n+1);
ky2map=zeros(2*n+1,2*n+1);
bcmap=zeros(2*n+1,2*n+1);

parfor xindex=1:length(xrange)
    x=xrange(xindex);
    kxlist=zeros(1,2*n+1);
    kylist=zeros(1,2*n+1);
    kx2list=zeros(1,2*n+1);
    ky2list=zeros(1,2*n+1);
    bclist=zeros(1,2*n+1);
    for yindex=1:length(yrange)        
        y=xrange(yindex);
        k=x*a1+y*a2;
        right=k+[bM1(1)/(2*n),0];
        left=k+[-bM1(1)/(2*n),0];
        up=k+[0,3*kn(2)/n];
        down=k+[0,3*kp(2)/n];
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
        bclist(yindex)=wilson(level,right,up,left,down,parameters);      
    end
    kxmap(xindex,:)=kxlist;
    kymap(xindex,:)=kylist;
    kx2map(xindex,:)=kx2list;
    ky2map(xindex,:)=ky2list;
    bcmap(xindex,:)=bclist;
end
end