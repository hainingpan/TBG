function [wbgrid,wtgrid]=w(R,rx,ry,parameters)
%R: center of wannier state
%r: scalar|array real space position

n=10;
state=1;
xrange=-n:n;
yrange=-n:n;
Nmax=parameters.Nmax;
% aM=parameters.aM;
bM1=parameters.bM1;
bM2=parameters.bM2;
%shift to diamond
a1=-bM1/(2*n);
a2=(bM1+bM2)/(2*n);
%shift to rectangular
% a1=bM2/(2*n);
% a2=(-2*bM1-bM2)/2/(2*n);
% [rx,ry]=meshgrid(linspace(-sqrt(3).aM,sqrt(3).aM,Nrx));
[Nrx,Nry]=size(rx);
kxmap=zeros(2*n+1,2*n+1);
kymap=zeros(2*n+1,2*n+1);
gauge=zeros(2*n+1,2*n+1);
expR=zeros(2*n+1,2*n+1);
psibline=zeros(2*n+1,2*n+1,Nrx*Nry);
psitline=zeros(2*n+1,2*n+1,Nrx*Nry);


Nkx=length(xrange);
Nky=length(yrange);
omega=abs(cross([bM1,0],[bM2,0]));
omega=omega(3);


parfor xindex=1:Nkx
    kx=xrange(xindex);
    for yindex=1:Nky
        ky=yrange(yindex);
        k=kx*a1+ky*a2;
        [~,vec]=energyTMD(k(1),k(2),parameters);
        psib0=sum(vec(1:(2*Nmax+1)^2,state)); %bloch wf at r=(0,0) on the bottom layer
        gauge(xindex,yindex)=conj(abs(psib0)/psib0);
        [ubgrid,utgrid]=u(vec(:,state),rx,ry,parameters);
        expR(xindex,yindex)=exp(-1i*dot(k,R));
        psibgrid=ubgrid.*exp(1i*(k(1)*rx+k(2)*ry));
        psitgrid=utgrid.*exp(1i*(k(1)*rx+k(2)*ry));
        psibline(xindex,yindex,:)=psibgrid(:);
        psitline(xindex,yindex,:)=psitgrid(:);        
        kxmap(xindex,yindex)=k(1);
        kymap(xindex,yindex)=k(2);
    end
end

expRgauge=(gauge.*expR);
%use integral
Fb=expRgauge.*psibline;
Ft=expRgauge.*psitline;
wbline=trapz(kxmap(:,1),trapz(kymap(1,:),Fb,2))/omega;
wtline=trapz(kxmap(:,1),trapz(kymap(1,:),Ft,2))/omega;

%use summation
psib=reshape(psibline,[Nkx*Nky,Nrx*Nry]);
psit=reshape(psitline,[Nkx*Nky,Nrx*Nry]);
wbline=expRgauge(:).'*psib/(Nkx*Nky);
wtline=expRgauge(:).'*psit/(Nkx*Nky);

wbgrid=reshape(wbline,[Nrx,Nry]);
wtgrid=reshape(wtline,[Nrx,Nry]);

end


