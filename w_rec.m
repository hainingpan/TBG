function [wbgrid,wtgrid]=w_rec(R,rx,ry,parameters)
%R: center of wannier state
%r: scalar|array real space position
%use rectangular (skew) grid
n=9;
state=1;
xrange=-n:n-1;
yrange=-n:n-1;
Nkx=length(xrange);
Nky=length(yrange);
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
% rectangular(skew grid)
kxmap=zeros(Nkx,Nky);
kymap=zeros(Nkx,Nky);
gauge=zeros(Nkx,Nky);
expR=zeros(Nkx,Nky);
psibline=zeros(Nkx,Nky,Nrx*Nry);
psitline=zeros(Nkx,Nky,Nrx*Nry);



omega=abs(cross([bM1,0],[bM2,0]));
omega=omega(3);

% rectangular(skew) grid, for diamond mesh and rectangular mesh
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
% 
%use integral, rectangular(skew) grid
% Fb=expRgauge.*psibline;
% Ft=expRgauge.*psitline;
% wbline=trapz(kxmap(:,1),trapz(kymap(1,:),Fb,2))/omega;
% wtline=trapz(kxmap(:,1),trapz(kymap(1,:),Ft,2))/omega;

%use summation, rectangular(skew) grid,
psib=reshape(psibline,[Nkx*Nky,Nrx*Nry]);
psit=reshape(psitline,[Nkx*Nky,Nrx*Nry]);
wbline=expRgauge(:).'*psib/(Nkx*Nky);
wtline=expRgauge(:).'*psit/(Nkx*Nky);

wbgrid=reshape(wbline,[Nrx,Nry]);
wtgrid=reshape(wtline,[Nrx,Nry]);

end


